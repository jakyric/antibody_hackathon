"""Developability triage, run locally before anything is sent to a GPU.

Why this sits *between* generation and scoring rather than at the end: a
co-fold costs ~75 s and the jobs queue, so with sixteen teams sharing the
platform the scarce resource is wall-clock, not ideas. A generated sequence
carrying a free cysteine or a new N-glycosylation sequon in CDR3 is not a
candidate a drug discovery team would take forward whatever its predicted
interface looks like, so it should never consume a slot.

Liabilities are weighted by *where* they fall. The same chemistry is a
different problem in CDR3, where it can cost binding and is hard to engineer
out, than in framework, where it usually can be left alone.

The CDR3 locator here uses the same anchoring rule as the abagent repo
(conserved [YF]C of framework 3, J motif at the end, re-anchor when the loop
comes out implausibly long). The measured accuracy quoted in that repo is for
*its* implementation against SAbDab2; this is a compact restatement of the
rule, not a re-measured engine. When AntiFold returns its per-residue region
CSV we use that instead, because it is the model's own assignment.
"""
import re

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}

PKA = {"Cterm": 2.34, "Nterm": 9.69, "C": 8.18, "D": 3.65, "E": 4.25,
       "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07}

NEGATIVE = ("D", "E", "C", "Y")
POSITIVE = ("K", "R", "H")

J_MOTIF = re.compile(r"[WF]G[A-Z]G")
FR3_ANCHOR = re.compile(r"[YF]C")

# kind -> (regex, weight, what it costs)
MOTIFS = {
    "n_glycosylation": (re.compile(r"N[^P][ST]"), 3,
                        "new glycan site: heterogeneity and a changed interface"),
    "deamidation_NG": (re.compile(r"N[GS]"), 2,
                       "Asn deamidation, fastest in NG"),
    "isomerisation_DG": (re.compile(r"D[GS]"), 1,
                         "Asp isomerisation"),
    "fragmentation_DP": (re.compile(r"DP"), 2,
                         "Asp-Pro acid-labile bond"),
    "oxidation_MW": (re.compile(r"[MW]"), 1,
                     "Met/Trp oxidation if solvent exposed"),
}

REGION_WEIGHT = {"CDR3": 3.0, "CDR1": 2.0, "CDR2": 2.0, "framework": 1.0}

# A liability a team would not carry forward at all, when it is in a loop.
HARD_FAIL = {"n_glycosylation", "free_cysteine"}


PLAUSIBLE_CDR3 = (3, 38)


def find_cdr3(seq):
    """(start, end) of heavy-chain CDR3, end exclusive, or None.

    End: the J motif that closes the loop. Start: just after the conserved
    cysteine (IMGT 104) that precedes it.

    Which cysteine, when there is more than one, is the whole difficulty. The
    conserved one is scanned for in ascending order and the first that leaves a
    plausible loop length wins, which handles both failure modes in one rule:

      * the intradomain partner (IMGT 23) sits ~80 residues upstream and leaves
        an implausibly long loop, so it is skipped;
      * an extra cysteine *inside* CDR3 sits downstream of the conserved one,
        so it is never reached and cannot truncate the loop.

    A regex anchor on the framework-3 `[YF]C` motif was tried first and
    rejected: real VHHs carry `TC` there often enough to matter -- the
    nanobody in the platform's own RECIPES.md ends framework 3 with `...TAIYTC`
    and would have been missed entirely.
    """
    js = list(J_MOTIF.finditer(seq))
    if not js:
        return None
    end = js[-1].start()
    lo, hi = PLAUSIBLE_CDR3
    for i, c in enumerate(seq[:end]):
        if c != "C":
            continue
        start = i + 1
        if lo <= end - start <= hi:
            return start, end
    return None


def regions_from_sequence(seq):
    """Position -> region label. CDR3 only; CDR1/CDR2 need real IMGT
    numbering and are reported as framework rather than guessed at."""
    labels = ["framework"] * len(seq)
    cdr3 = find_cdr3(seq)
    if cdr3:
        for i in range(cdr3[0], cdr3[1]):
            labels[i] = "CDR3"
    return labels


def regions_from_antifold_csv(path, chain=None):
    """AntiFold returns per-residue perplexity and the region it assumed.
    That assignment is better than ours, so prefer it when present."""
    labels = {}
    try:
        with open(path) as fh:
            header = fh.readline().strip().split(",")
            cols = {name.strip().lower(): n for n, name in enumerate(header)}
            ri = next((cols[k] for k in ("region", "imgt_region", "loop") if k in cols), None)
            pi = next((cols[k] for k in ("pdb_pos", "position", "pos", "resi") if k in cols), None)
            ci = next((cols[k] for k in ("chain", "pdb_chain") if k in cols), None)
            if ri is None or pi is None:
                return {}
            for line in fh:
                parts = line.rstrip("\n").split(",")
                if len(parts) <= max(ri, pi):
                    continue
                if chain and ci is not None and parts[ci].strip() != chain:
                    continue
                labels[parts[pi].strip()] = (parts[ri].strip() or "framework")
    except OSError:
        return {}
    return labels


def liabilities(seq, labels=None):
    labels = labels or regions_from_sequence(seq)
    found = []
    for kind, (rx, weight, cost) in MOTIFS.items():
        for m in rx.finditer(seq):
            i = m.start()
            region = labels[i] if i < len(labels) else "framework"
            if kind == "oxidation_MW" and region == "framework":
                continue  # buried Met/Trp in framework is normal, not a finding
            found.append({
                "kind": kind, "motif": m.group(0), "pos": i + 1,
                "region": region, "cost": cost,
                "penalty": weight * REGION_WEIGHT.get(region, 1.0),
            })
    # Cysteines beyond the conserved intradomain pair are free thiols: they
    # drive dimerisation and mispairing, and they are the one liability that is
    # almost always disqualifying.
    cys = [i for i, c in enumerate(seq) if c == "C"]
    if len(cys) > 2:
        for i in cys[2:]:
            region = labels[i] if i < len(labels) else "framework"
            found.append({
                "kind": "free_cysteine", "motif": "C", "pos": i + 1,
                "region": region, "cost": "unpaired thiol: dimerisation, mispairing",
                "penalty": 3.0 * REGION_WEIGHT.get(region, 1.0),
            })
    return sorted(found, key=lambda d: -d["penalty"])


def net_charge(seq, ph=7.0):
    q = 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"])) - 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
    for aa in POSITIVE:
        q += seq.count(aa) / (1.0 + 10 ** (ph - PKA[aa]))
    for aa in NEGATIVE:
        q -= seq.count(aa) / (1.0 + 10 ** (PKA[aa] - ph))
    return q


def isoelectric_point(seq):
    lo, hi = 2.0, 13.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def gravy(seq):
    vals = [KD[c] for c in seq if c in KD]
    return round(sum(vals) / len(vals), 3) if vals else 0.0


def assess(name, seq, labels=None):
    items = liabilities(seq, labels)
    in_loop = [d for d in items if d["region"].startswith("CDR")]
    blocking = [d for d in in_loop if d["kind"] in HARD_FAIL]
    cdr3 = find_cdr3(seq)
    return {
        "name": name,
        "sequence": seq,
        "length": len(seq),
        "cdr3": seq[cdr3[0]:cdr3[1]] if cdr3 else "",
        "pI": isoelectric_point(seq),
        "charge_pH7": round(net_charge(seq), 2),
        "gravy": gravy(seq),
        "liabilities": items,
        "n_liabilities": len(items),
        "n_in_cdr": len(in_loop),
        "penalty": round(sum(d["penalty"] for d in items), 1),
        "blocking": [d["kind"] + "@" + str(d["pos"]) for d in blocking],
        "pass": not blocking,
    }


def triage(records, keep, region_labels=None):
    """Assess, drop duplicates and hard failures, keep the cleanest `keep`.

    Returns (kept, rejected). Both carry their full assessment so the report
    can say why something was dropped -- a candidate list without the
    rejections is not auditable.
    """
    seen, assessed = set(), []
    for name, seq in records:
        seq = "".join(seq.split()).upper()
        if not seq or seq in seen:
            continue
        seen.add(seq)
        assessed.append(assess(name, seq, (region_labels or {}).get(name)))
    passing = [a for a in assessed if a["pass"]]
    failing = [a for a in assessed if not a["pass"]]
    passing.sort(key=lambda a: (a["penalty"], a["n_in_cdr"]))
    return passing[:keep], failing + passing[keep:]
