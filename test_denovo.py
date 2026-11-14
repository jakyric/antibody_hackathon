"""Checks for the parts that can be checked without the platform.

Geometry is tested against hand-placed coordinates, so a wrong contact count
fails here rather than silently producing a wrong ranking.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from denovo import bio, report, triage  # noqa: E402

# The VHH from the platform's RECIPES.md (1ZVY chain A plus a His tag).
VHH = ("DVQLVESGGGSVQAGGSLRLSCAASGSTDSIEYMTWFRQAPGKAREGVAALYTHTGNTYYTDSVKGRFT"
       "ISQDKAKNMAYLRMDSVKSEDTAIYTCGATRKYVPVRFALDQSSYDYWGQGTQVTVSSRGRHHHHHH")


def pdb_line(serial, name, resname, chain, resseq, x, y, z, element):
    return ("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, name, resname, chain, resseq, x, y, z, element))


class TestCdr3(unittest.TestCase):
    def test_real_vhh(self):
        span = triage.find_cdr3(VHH)
        self.assertIsNotNone(span, "CDR3 not found in a real VHH")
        cdr3 = VHH[span[0]:span[1]]
        self.assertEqual(cdr3, "GATRKYVPVRFALDQSSYDY")
        self.assertTrue(VHH[span[0] - 1] == "C")
        self.assertTrue(VHH[span[1]:span[1] + 4].startswith("WGQG"))

    def test_extra_cysteine_inside_cdr3_does_not_truncate(self):
        # Put a cysteine inside the loop; the conserved one must still win.
        mutated = VHH.replace("GATRKYVPVRF", "GATCKYVPVRF")
        span = triage.find_cdr3(mutated)
        self.assertIsNotNone(span)
        self.assertEqual(len(mutated[span[0]:span[1]]), 20)

    def test_no_j_motif(self):
        self.assertIsNone(triage.find_cdr3("ACDEFGHIKLMNPQRSTVWY" * 5))


class TestLiabilities(unittest.TestCase):
    def test_glycosylation_in_cdr3_is_blocking(self):
        # NIS is an N-glycosylation sequon; place it inside the loop.
        bad = VHH.replace("GATRKYVPVRF", "GANISYVPVRF")
        a = triage.assess("bad", bad)
        self.assertFalse(a["pass"])
        self.assertTrue(any(b.startswith("n_glycosylation") for b in a["blocking"]))

    def test_clean_vhh_passes(self):
        a = triage.assess("vhh", VHH)
        self.assertTrue(a["pass"], a["blocking"])
        self.assertEqual(a["cdr3"], "GATRKYVPVRFALDQSSYDY")

    def test_free_cysteine_flagged(self):
        seq = VHH[:60] + "C" + VHH[61:]
        kinds = {d["kind"] for d in triage.liabilities(seq)}
        self.assertIn("free_cysteine", kinds)

    def test_framework_oxidation_not_reported(self):
        # Met in framework is normal and must not inflate the penalty.
        for d in triage.liabilities(VHH):
            if d["kind"] == "oxidation_MW":
                self.assertTrue(d["region"].startswith("CDR"))

    def test_biophysics_sane(self):
        self.assertTrue(3.0 < triage.isoelectric_point(VHH) < 11.0)
        self.assertLess(triage.gravy(VHH), 1.0)
        self.assertGreater(triage.net_charge("KKKKK"), 3.0)
        self.assertLess(triage.net_charge("DDDDD"), -3.0)

    def test_triage_drops_and_keeps(self):
        bad = VHH.replace("GATRKYVPVRF", "GANISYVPVRF")
        kept, dropped = triage.triage([("good", VHH), ("bad", bad)], keep=5)
        self.assertEqual([a["name"] for a in kept], ["good"])
        self.assertEqual([a["name"] for a in dropped], ["bad"])

    def test_triage_dedupes(self):
        kept, _ = triage.triage([("a", VHH), ("b", VHH)], keep=5)
        self.assertEqual(len(kept), 1)


class TestGeometry(unittest.TestCase):
    """Two chains with atoms placed by hand, so the expected counts are known."""

    def _structure(self, gap):
        lines = []
        # chain A: one residue, 2 heavy atoms at x=0 and x=1
        lines.append(pdb_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"))
        lines.append(pdb_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C"))
        # chain T: one residue, 2 heavy atoms offset along x by `gap`
        lines.append(pdb_line(3, "CA", "GLY", "T", 10, gap, 0.0, 0.0, "C"))
        lines.append(pdb_line(4, "CB", "GLY", "T", 10, gap + 1.0, 0.0, 0.0, "C"))
        path = os.path.join(self.tmp, "s%s.pdb" % gap)
        with open(path, "w") as fh:
            fh.writelines(lines)
        return path

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = self.tmpdir.name

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_contacts_counted_exactly(self):
        atoms = bio.load_structure(self._structure(2.0))
        # pairs within 4.5 A: (0,2.0)=2.0 yes, (0,3.0)=3.0 yes,
        #                     (1,2.0)=1.0 yes, (1,3.0)=2.0 yes  -> 4
        pairs, ra, rt, respairs = bio.interchain_contacts(atoms, "A", "T")
        self.assertEqual(pairs, 4)
        self.assertEqual(len(ra), 1)
        self.assertEqual(len(rt), 1)
        self.assertEqual(respairs, 1)

    def test_no_contacts_when_far(self):
        atoms = bio.load_structure(self._structure(40.0))
        pairs, ra, rt, _ = bio.interchain_contacts(atoms, "A", "T")
        self.assertEqual(pairs, 0)
        self.assertEqual(ra, set())
        self.assertEqual(rt, set())

    def test_cutoff_boundary(self):
        atoms = bio.load_structure(self._structure(5.0))
        # (1, 5.0) = 4.0 yes ; (1,6.0)=5.0 no ; (0,5.0)=5.0 no ; (0,6.0)=6.0 no
        pairs, _, _, _ = bio.interchain_contacts(atoms, "A", "T")
        self.assertEqual(pairs, 1)

    def test_water_and_hydrogen_ignored(self):
        path = self._structure(2.0)
        with open(path, "a") as fh:
            fh.write(pdb_line(5, "O", "HOH", "T", 11, 1.0, 0.0, 0.0, "O"))
            fh.write(pdb_line(6, "H1", "ALA", "A", 1, 1.2, 0.0, 0.0, "H"))
        atoms = bio.load_structure(path)
        self.assertFalse(any(a.resname == "HOH" for a in atoms))
        pairs, _, _, _ = bio.interchain_contacts(atoms, "A", "T")
        self.assertEqual(pairs, 4)  # the hydrogen must not add a contact

    def test_cif_round_trip(self):
        cif = os.path.join(self.tmp, "x.cif")
        with open(cif, "w") as fh:
            fh.write("data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n"
                     "_atom_site.type_symbol\n_atom_site.label_atom_id\n"
                     "_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
                     "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n"
                     "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
                     "ATOM 1 C CA . ALA A 1 0.000 0.000 0.000\n"
                     "ATOM 2 C CA . GLY T 9 2.000 0.000 0.000\n"
                     "ATOM 3 O O . HOH T 10 2.500 0.000 0.000\n")
        atoms = bio.parse_cif(cif)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(bio.sequence_of(atoms, "A"), "A")
        self.assertEqual(bio.sequence_of(atoms, "T"), "G")
        pairs, _, _, _ = bio.interchain_contacts(atoms, "A", "T")
        self.assertEqual(pairs, 1)

    def test_epitope_overlap(self):
        a = {("1", "ALA"), ("2", "GLY"), ("3", "SER")}
        b = {("2", ""), ("3", ""), ("4", "")}
        frac, shared = bio.epitope_overlap(a, b)
        self.assertEqual(shared, 2)
        self.assertAlmostEqual(frac, 0.5)
        self.assertEqual(bio.epitope_overlap(set(), b), (0.0, 0))

    def test_antibody_heuristic(self):
        self.assertTrue(bio.looks_like_antibody(VHH))
        self.assertFalse(bio.looks_like_antibody("A" * 200))
        self.assertFalse(bio.looks_like_antibody(VHH[:40]))

    def test_tar_has_flat_members(self):
        import tarfile
        src = self._structure(2.0)
        tar = bio.make_tar(os.path.join(self.tmp, "t.tar"), [src])
        with tarfile.open(tar) as t:
            names = t.getnames()
        self.assertEqual(names, [os.path.basename(src)])
        self.assertFalse(any(n.startswith("./") for n in names))

    def test_fasta_round_trip(self):
        p = bio.write_fasta(os.path.join(self.tmp, "a.fasta"), [("H", VHH)])
        self.assertEqual(bio.read_fasta(p), [("H", VHH)])


class TestReport(unittest.TestCase):
    def test_writes_and_ranks(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = {
                "reference_epitope": ["1", "2", "3"],
                "control_score": {"confidence": 0.77, "ptm": 0.57, "contacts": 200},
                "scores": [
                    {"name": "weak", "op": "op-1", "sequence": VHH, "confidence": 0.60,
                     "ptm": 0.40, "contacts": 20, "interface_residues": 8,
                     "epitope_overlap": 0.1, "penalty": 6.0, "pI": 8.1},
                    {"name": "strong", "op": "op-2", "sequence": VHH, "confidence": 0.85,
                     "ptm": 0.65, "contacts": 180, "interface_residues": 40,
                     "epitope_overlap": 0.6, "penalty": 1.0, "pI": 7.4},
                ],
                "dropped": [{"name": "nope", "blocking": ["n_glycosylation@101"],
                             "penalty": 9.0}],
            }
            csv_path, md_path = report.write(state, tmp)
            self.assertTrue(os.path.exists(csv_path))
            body = open(md_path).read()
            self.assertIn("strong", body)
            self.assertIn("op-2", body)
            self.assertIn("Not carried forward", body)
            # the better candidate must rank first
            self.assertEqual([r["name"] for r in sorted(
                state["scores"], key=lambda r: r["rank"])], ["strong", "weak"])

    def test_report_survives_missing_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            state = {"scores": [{"name": "a", "op": "op-1", "sequence": "AC",
                                 "contacts": None, "confidence": None}]}
            _, md = report.write(state, tmp)
            self.assertIn("No control was scored", open(md).read())


class TestCli(unittest.TestCase):
    """The CLI must plan a whole campaign without credentials or network."""

    def test_dry_run_plans_without_creds(self):
        with tempfile.TemporaryDirectory() as tmp:
            stub = os.path.join(tmp, "forge.py")
            open(stub, "w").write("#!/usr/bin/env python3\n")
            os.makedirs(os.path.join(tmp, "results"), exist_ok=True)
            # Seed the state so the stages that need a target can run offline.
            with open(os.path.join(tmp, "results", "campaign.json"), "w") as fh:
                json.dump({"target_seq": "GGGG" * 20,
                           "designs": [{"name": "d-01", "sequence": VHH,
                                        "source_header": "h", "generator": "g",
                                        "op": "op-x"}]}, fh)
            proc = subprocess.run(
                [sys.executable, "-m", "denovo", "--forge", stub, "--dry-run",
                 "--work-dir", tmp, "--stages", "triage,cofold,report"],
                cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
            self.assertEqual(proc.returncode, 0, proc.stdout)
            self.assertIn("[dry-run]", proc.stdout)
            self.assertIn("boltz2", proc.stdout)
            self.assertIn("keep", proc.stdout)

    def test_requires_creds_when_not_dry_run(self):
        env = dict(os.environ)
        env.pop("NEBIUS_CREDS", None)
        proc = subprocess.run(
            [sys.executable, "-m", "denovo", "--forge", "/nope/forge.py"],
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, timeout=60)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("no credentials", proc.stdout)

    def test_rejects_unknown_stage(self):
        proc = subprocess.run(
            [sys.executable, "-m", "denovo", "--dry-run", "--stages", "wat"],
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("unknown stage", proc.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
