"""The campaign, stage by stage.

Design shape, and why:

  target    1N8Z gives the HER2 domain IV sequence, a known binder
            (trastuzumab) as a positive control, and the epitope measured off
            the crystal structure rather than quoted from memory.
  msa       one alignment for the target, reused by every later co-fold.
  baseline  trastuzumab + HER2 co-folded first. Without a yardstick a
            confidence of 0.74 means nothing; against a known binder scored
            the same way, it means something.
  generate  new CDRs on a VHH framework (AntiFold, antibody-aware), or a
            general binder backbone (RFD4-Proteina), or de novo from the
            target sequence alone (Mosaic). Default is AntiFold: it is the
            only one of the three with an antibody prior, and it is 110 s
            against 5 and 22 minutes.
  triage    local, free, instant. Drops candidates no team would carry
            forward before they consume a GPU slot.
  build     NanoBodyBuilder2 per survivor: an antibody-aware model in 130 s.
  cofold    Boltz2, candidate against target. Submitted with --no-wait so the
            whole panel queues at once instead of serially.
  score     Boltz2 reports no ipTM, so the interface is measured from the
            returned coordinates: contacts, interface size, and whether the
            candidate sits on the trastuzumab epitope.
  confirm   OpenFold3 on the top few. Two models agreeing is worth more than
            one high number.

Every stage is resumable. Run ids live in forge.py's receipts, campaign state
in results/campaign.json, so the laptop can close mid-run.
"""
import json
import os
import re

from . import bio, config, triage
from .forge_client import ForgeError

STATE_FILE = "campaign.json"


def slug(text):
    return re.sub(r"[^a-z0-9]+", "-", str(text).lower()).strip("-")[:48]


class Campaign:
    def __init__(self, forge, work_dir, target_chain=None, heavy=None, light=None):
        self.f = forge
        self.work = work_dir
        self.inputs = os.path.join(work_dir, config.INPUTS_DIR)
        self.results = os.path.join(work_dir, config.RESULTS_DIR)
        os.makedirs(self.inputs, exist_ok=True)
        os.makedirs(self.results, exist_ok=True)
        self.state = self._load()
        self.state.setdefault("chains", {})
        if target_chain:
            self.state["chains"]["target"] = target_chain
        if heavy:
            self.state["chains"]["heavy"] = heavy
        if light:
            self.state["chains"]["light"] = light

    # ------------------------------------------------------------------ state

    def _path(self):
        return os.path.join(self.results, STATE_FILE)

    def _load(self):
        try:
            with open(self._path()) as fh:
                return json.load(fh)
        except (OSError, json.JSONDecodeError):
            return {}

    def save(self):
        with open(self._path(), "w") as fh:
            json.dump(self.state, fh, indent=1, sort_keys=True)

    def rel(self, path):
        return os.path.relpath(path, self.work)

    # ----------------------------------------------------------- 1. the target

    def stage_target(self):
        print("[target] %s from %s" % (config.TARGET_NAME, config.TARGET_PDB.upper()))
        pdb = bio.fetch_pdb(config.TARGET_PDB, self.inputs)
        atoms = bio.load_structure(pdb)
        seqs = {c: "".join(o for _, o in res) for c, res in bio.chains(atoms).items()}
        for cid, seq in sorted(seqs.items()):
            print("  chain %s  %3d aa  %s" % (
                cid, len(seq), "antibody-like" if bio.looks_like_antibody(seq) else "target-like"))

        chosen = self.state["chains"]
        ab = [c for c, s in seqs.items() if bio.looks_like_antibody(s)]
        non_ab = [c for c in seqs if c not in ab]
        if "target" not in chosen:
            if not non_ab:
                raise ForgeError("no non-antibody chain found in %s; pass --target-chain"
                                 % config.TARGET_PDB)
            chosen["target"] = max(non_ab, key=lambda c: len(seqs[c]))
        if "heavy" not in chosen and len(ab) >= 2:
            # The heavy chain of a Fab is the longer of the two V-domain chains
            # in 1N8Z; override with --heavy-chain if that is not true here.
            pair = sorted(ab, key=lambda c: -len(seqs[c]))
            chosen["heavy"], chosen["light"] = pair[0], pair[1]

        t = chosen["target"]
        self.state["target_seq"] = seqs[t]
        print("  target chain %s, %d aa" % (t, len(seqs[t])))

        bio.write_fasta(os.path.join(self.inputs, "target.fasta"),
                        [(config.TARGET_NAME, seqs[t])])

        # The reference epitope: which target residues the known binder touches.
        if chosen.get("heavy") and chosen.get("light"):
            ep = set()
            for ab_chain in (chosen["heavy"], chosen["light"]):
                _, _, res_t, _ = bio.interchain_contacts(atoms, ab_chain, t)
                ep |= res_t
            self.state["reference_epitope"] = sorted(r[0] for r in ep)
            self.state["reference_seqs"] = {
                "H": seqs[chosen["heavy"]], "L": seqs[chosen["light"]]}
            print("  trastuzumab epitope on chain %s: %d residues (measured from the structure)"
                  % (t, len(ep)))
        self.save()
        return self.state["target_seq"]

    # -------------------------------------------------------------- 2. the MSA

    def stage_msa(self):
        model, _lane, est = config.MODELS["msa"]
        label = "g06-msa-target"
        print("[msa] alignment for the target")
        op = self.f.run(model, {
            "sequence": self.state["target_seq"],
            "databases": [config.MSA_DB],
            "output_alignment_formats": ["a3m"],
        }, label=label, est=est)
        self.f.wait(op)
        dest = self.f.fetch(op, os.path.join(config.RESULTS_DIR, label))
        a3m = None
        if dest:
            hits = bio.find_files(dest, (".a3m",))
            if hits:
                a3m = self.rel(hits[0])
        if not a3m:
            print("  no a3m came back; falling back to a single-sequence alignment")
        self.state["target_a3m"] = a3m
        self.state["ops"] = dict(self.state.get("ops", {}), msa=op)
        self.save()
        return a3m

    def _msa_block(self, seq, a3m_rel=None):
        """Boltz2 wants an alignment per chain. A real a3m when we have one,
        otherwise the one-line form, which is enough for a first look."""
        if a3m_rel:
            return {"msa_search": {"a3m": {"alignment": {"$text": a3m_rel}}}}
        return {"msa_search": {"a3m": {"alignment": ">query\n%s\n" % seq}}}

    # ---------------------------------------------- 3. baseline and target fold

    def _cofold(self, label, candidate_seq, est=None, no_wait=True, extra=None):
        model, _lane, default_est = config.MODELS["fold"]
        polymers = [{"id": "A", "molecule_type": "protein", "sequence": candidate_seq,
                     "msa": self._msa_block(candidate_seq)}]
        for cid, seq in (extra or []):
            polymers.append({"id": cid, "molecule_type": "protein", "sequence": seq,
                             "msa": self._msa_block(seq)})
        polymers.append({"id": "T", "molecule_type": "protein",
                         "sequence": self.state["target_seq"],
                         "msa": self._msa_block(self.state["target_seq"],
                                                self.state.get("target_a3m"))})
        return self.f.run(model, {"polymers": polymers}, label=label,
                          no_wait=no_wait, est=est or default_est)

    def stage_baseline(self):
        """Fold the target alone, and co-fold the known binder against it."""
        model, _lane, est = config.MODELS["fold"]
        print("[baseline] target alone, and the known binder against it")
        ops = {}
        ops["target_fold"] = self.f.run(model, {"polymers": [{
            "id": "T", "molecule_type": "protein", "sequence": self.state["target_seq"],
            "msa": self._msa_block(self.state["target_seq"], self.state.get("target_a3m"))}]},
            label="g06-fold-target", no_wait=True, est=est)
        ref = self.state.get("reference_seqs")
        if ref:
            ops["control"] = self._cofold("g06-control-trastuzumab", ref["H"],
                                          extra=[("B", ref["L"])])
        for name, op in ops.items():
            self.f.wait(op)
            self.f.fetch(op, os.path.join(config.RESULTS_DIR, "g06-" + slug(name)))
        self.state["ops"] = dict(self.state.get("ops", {}), **ops)
        self.save()
        return ops

    # ---------------------------------------------------------- 4. generation

    def stage_generate(self, route="antifold", n=config.N_GENERATE):
        print("[generate] route=%s, asking for %d sequences" % (route, n))
        if route == "antifold":
            return self._generate_antifold(n)
        if route == "rfd4":
            return self._generate_rfd4()
        if route == "mosaic":
            return self._generate_mosaic()
        raise ForgeError("unknown route %r (antifold, rfd4, mosaic)" % route)

    def _generate_antifold(self, n):
        """New CDR loops on a VHH framework. Antibody-aware, 110 s."""
        model, _lane, est = config.MODELS["redesign"]
        scaffold = bio.fetch_pdb(config.SCAFFOLD_PDB, self.inputs)
        tar = bio.make_tar(os.path.join(self.inputs, "antifold.tar"), [scaffold])
        label = "g06-antifold-cdrs"
        op = self.f.submit(model,
                           [("workload", "inverse-folding-tar/v1", self.rel(tar))],
                           {"pdb_file": os.path.basename(scaffold),
                            "nanobody_chain": config.SCAFFOLD_CHAIN,
                            "nanobody_mode": True,
                            "regions": config.CDR_REGIONS,
                            "num_seq_per_target": n,
                            "sampling_temp": config.SAMPLING_TEMP},
                           label=label, est=est)
        self.f.wait(op)
        dest = self.f.fetch(op, os.path.join(config.RESULTS_DIR, label))
        return self._read_designs(dest, label, op)

    def _generate_rfd4(self):
        """A general binder backbone guided to the target, then a sequence for
        it. RFD4-Proteina has no antibody prior: what comes out is a protein
        binder, not a VHH, and it needs ProteinMPNN after it."""
        model, _lane, est = config.MODELS["backbone"]
        cif = os.path.join(self.inputs, "%s.cif" % config.TARGET_PDB)
        if not os.path.exists(cif):
            raise ForgeError(
                "RFD4-Proteina needs an annotated mmCIF of the target at %s "
                "(the annotation marks the site to design against). Produce it "
                "first, or use --route antifold." % self.rel(cif))
        label = "g06-rfd4-backbone"
        op = self.f.submit(model,
                           [("annotated-cif", "rfd4-proteina-annotated-cif/v1", self.rel(cif))],
                           {"num_replicates": 1, "batch_size": 1, "seed": 42,
                            "nsteps": 60, "recycle_strategy": "thinking_medium",
                            "compile": True, "save_trajectory_every": 0,
                            "save_chain_of_interest": False},
                           label=label, est=est)
        self.f.wait(op)
        dest = self.f.fetch(op, os.path.join(config.RESULTS_DIR, label))
        backbones = bio.find_files(dest or "", (".pdb", ".cif")) if dest else []
        if not backbones:
            raise ForgeError("no backbone came back from %s" % model)
        mpnn, _lane, mest = config.MODELS["seq_design"]
        designs = []
        for i, bb in enumerate(backbones[:4], 1):
            dlabel = "g06-mpnn-%d" % i
            dop = self.f.run(mpnn, {
                "input_pdb": {"$upload": self.rel(bb), "media_type": "chemical/x-pdb"},
                "input_pdb_chains": ["A"], "num_seq_per_target": 4,
                "sampling_temp": 0.1, "random_seed": 42,
            }, label=dlabel, est=mest)
            self.f.wait(dop)
            d = self.f.fetch(dop, os.path.join(config.RESULTS_DIR, dlabel))
            designs += self._read_designs(d, dlabel, dop)
        return designs

    def _generate_mosaic(self):
        """De novo binder from the target sequence alone. 22 minutes, and it has
        no tested recipe here, so the entry comes from the model's own
        contract rather than from a guess."""
        model, _lane, est = config.MODELS["binder"]
        entries = self.f.contract(model, "design-binder")
        if not entries or not entries[0].get("name"):
            raise ForgeError(
                "Mosaic publishes no input contract I can read. Run\n"
                "  python3 forge.py schema mosaic --full\n"
                "and pass the entry explicitly, or use --route antifold.")
        fasta = bio.write_fasta(os.path.join(self.inputs, "target.fasta"),
                                [(config.TARGET_NAME, self.state["target_seq"])])
        e = entries[0]
        label = "g06-mosaic-binder"
        op = self.f.submit(model,
                           [(e["name"], e["semantic_type"], self.rel(fasta))],
                           {}, label=label, est=est)
        self.f.wait(op)
        dest = self.f.fetch(op, os.path.join(config.RESULTS_DIR, label))
        return self._read_designs(dest, label, op)

    def _read_designs(self, dest, label, op):
        """Pull sequences out of whatever the generator returned, and keep the
        per-residue region CSV when there is one."""
        if not dest:
            return []
        records = []
        for path in bio.find_files(dest, (".fasta", ".fa", ".mfasta", ".fas")):
            records += bio.read_fasta(path)
        csvs = bio.find_files(dest, (".csv",))
        regions = {}
        if csvs:
            labels = triage.regions_from_antifold_csv(csvs[0], config.SCAFFOLD_CHAIN)
            if labels:
                regions["_by_position"] = labels
        out = []
        for i, (header, seq) in enumerate(records, 1):
            seq = "".join(seq.split()).upper()
            if len(seq) < 60:
                continue
            out.append(("%s-%02d" % (label.replace("g06-", ""), i), seq, header))
        self.state["designs"] = [{"name": n, "sequence": s, "source_header": h,
                                 "generator": label, "op": op} for n, s, h in out]
        self.state["design_regions"] = regions
        self.save()
        print("  %d sequences returned" % len(out))
        return out

    # ------------------------------------------------------------- 5. triage

    def stage_triage(self, keep=config.N_COFOLD):
        designs = self.state.get("designs") or []
        if not designs:
            raise ForgeError("no designs to triage; run the generate stage first")
        print("[triage] %d designs, keeping the %d cleanest" % (len(designs), keep))
        kept, dropped = triage.triage([(d["name"], d["sequence"]) for d in designs], keep)
        for a in kept:
            print("  keep %-22s penalty %-5s CDR3 %s" % (a["name"], a["penalty"], a["cdr3"]))
        for a in dropped:
            why = ", ".join(a["blocking"]) or "penalty %s" % a["penalty"]
            print("  drop %-22s %s" % (a["name"], why))
        self.state["kept"] = kept
        self.state["dropped"] = dropped
        self.save()
        return kept

    # ------------------------------------------------- 6/7. build and co-fold

    def stage_build(self):
        """An antibody-aware model per survivor. Useful on its own, and the
        input if you later want structure-aware scoring."""
        model, _lane, est = config.MODELS["build_vhh"]
        kept = self.state.get("kept") or []
        print("[build] %s on %d candidates" % (model, len(kept)))
        ops = {}
        for a in kept:
            fasta = bio.write_fasta(
                os.path.join(self.inputs, "%s.fasta" % slug(a["name"])),
                [("H", a["sequence"])])
            label = "g06-nbb2-%s" % slug(a["name"])
            ops[a["name"]] = self.f.submit(
                model, [("sequence", "nanobody-fasta/v1", self.rel(fasta))],
                {"numbering_scheme": "imgt"}, label=label, no_wait=True, est=est)
        for name, op in ops.items():
            self.f.wait(op)
            self.f.fetch(op, os.path.join(config.RESULTS_DIR, "g06-nbb2-" + slug(name)))
        self.state["build_ops"] = ops
        self.save()
        return ops

    def stage_cofold(self):
        """The expensive stage. Everything is queued first, then collected."""
        kept = self.state.get("kept") or []
        print("[cofold] %d candidates against %s, ~75 s each once a slot opens"
              % (len(kept), config.TARGET_NAME))
        ops = {}
        for a in kept:
            ops[a["name"]] = self._cofold("g06-cofold-%s" % slug(a["name"]), a["sequence"])
        dests = {}
        for name, op in ops.items():
            self.f.wait(op)
            dests[name] = self.f.fetch(
                op, os.path.join(config.RESULTS_DIR, "g06-cofold-" + slug(name)))
        self.state["cofold_ops"] = ops
        self.save()
        return dests

    # -------------------------------------------------------------- 8. scoring

    def stage_score(self):
        kept = {a["name"]: a for a in (self.state.get("kept") or [])}
        ref_ep = set(self.state.get("reference_epitope") or [])
        rows = []
        print("[score] measuring each interface from the returned coordinates")
        for name, op in (self.state.get("cofold_ops") or {}).items():
            dest = os.path.join(self.results, "g06-cofold-" + slug(name))
            row = {"name": name, "op": op,
                   "sequence": kept.get(name, {}).get("sequence", ""),
                   "cdr3": kept.get(name, {}).get("cdr3", ""),
                   "penalty": kept.get(name, {}).get("penalty"),
                   "pI": kept.get(name, {}).get("pI"),
                   "gravy": kept.get(name, {}).get("gravy")}
            row.update(self._score_one(dest, ref_ep))
            rows.append(row)
            print("  %-22s conf %-6s pTM %-6s contacts %-5s epitope overlap %s"
                  % (name, row.get("confidence"), row.get("ptm"),
                     row.get("contacts"), row.get("epitope_overlap")))
        control = self._score_one(os.path.join(self.results, "g06-control-trastuzumab"), ref_ep)
        self.state["scores"] = rows
        self.state["control_score"] = control
        if control.get("confidence") is not None:
            print("  %-22s conf %-6s pTM %-6s contacts %-5s  <- known binder, same scoring"
                  % ("trastuzumab (control)", control.get("confidence"),
                     control.get("ptm"), control.get("contacts")))
        self.save()
        return rows

    def _score_one(self, dest, ref_ep):
        out = {"confidence": None, "ptm": None, "contacts": None,
               "interface_residues": None, "epitope_overlap": None,
               "epitope_shared": None, "structure": None}
        if not dest or not os.path.isdir(dest):
            return out
        for path in bio.find_files(dest, (".json",)):
            with open(path) as fh:
                try:
                    doc = json.load(fh)
                except json.JSONDecodeError:
                    continue
            conf = _first_number(doc, ("confidence_score", "confidence"))
            ptm = _first_number(doc, ("ptm_score", "ptm"))
            if conf is not None and out["confidence"] is None:
                out["confidence"] = round(conf, 4)
            if ptm is not None and out["ptm"] is None:
                out["ptm"] = round(ptm, 4)
        structures = bio.find_files(dest, (".cif", ".mmcif", ".pdb"))
        if structures:
            out["structure"] = os.path.basename(structures[0])
            atoms = bio.load_structure(structures[0])
            present = {a.chain for a in atoms}
            target_chain = "T" if "T" in present else sorted(present)[-1]
            binder_chains = [c for c in sorted(present) if c != target_chain]
            pairs, res_b, res_t, respairs = 0, set(), set(), 0
            for c in binder_chains:
                p, rb, rt, rp = bio.interchain_contacts(atoms, c, target_chain)
                pairs += p
                res_b |= rb
                res_t |= rt
                respairs += rp
            out["contacts"] = pairs
            out["interface_residues"] = len(res_b) + len(res_t)
            if ref_ep:
                frac, shared = bio.epitope_overlap(res_t, {(r, "") for r in ref_ep})
                out["epitope_overlap"] = round(frac, 3)
                out["epitope_shared"] = shared
        return out

    # ------------------------------------------------------------- 9. confirm

    def stage_confirm(self, k=config.N_CONFIRM):
        model, _lane, est = config.MODELS["confirm"]
        rows = sorted((self.state.get("scores") or []),
                      key=lambda r: -(r.get("contacts") or 0))[:k]
        print("[confirm] %s on the top %d" % (model, len(rows)))
        ops = {}
        for r in rows:
            label = "g06-of3-%s" % slug(r["name"])
            ops[r["name"]] = self.f.run(model, {
                "request_id": slug(r["name"]),
                "inputs": [{"input_id": slug(r["name"]), "output_format": "cif",
                            "molecules": [
                                {"id": "A", "type": "protein", "sequence": r["sequence"]},
                                {"id": "T", "type": "protein",
                                 "sequence": self.state["target_seq"]}]}],
            }, label=label, no_wait=True, est=est)
        for name, op in ops.items():
            self.f.wait(op)
            self.f.fetch(op, os.path.join(config.RESULTS_DIR, "g06-of3-" + slug(name)))
        self.state["confirm_ops"] = ops
        self.save()
        return ops


def _first_number(doc, keys):
    """Find the first numeric value under any key containing one of `keys`.

    The result payloads nest differently per model and the field names move
    between versions, so matching on substring is more durable than a path.
    """
    stack = [doc]
    while stack:
        node = stack.pop(0)
        if isinstance(node, dict):
            for k, v in node.items():
                kl = str(k).lower()
                if any(key in kl for key in keys):
                    if isinstance(v, (int, float)):
                        return float(v)
                    if isinstance(v, list) and v and isinstance(v[0], (int, float)):
                        return float(sum(v)) / len(v)
                stack.append(v)
        elif isinstance(node, list):
            stack.extend(node)
    return None
