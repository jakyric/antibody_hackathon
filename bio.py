"""Structure and sequence handling. Standard library only.

Deliberately no Biopython: the whole point of the hackathon brief is that a
demo should not fail because a dependency would not install, and everything
needed here is a few hundred lines of parsing.

The one piece of real science in this file is `interchain_contacts`. Boltz2's
output has no ipTM, so the interface has to be measured from the returned
coordinates rather than read off a score. That is also what lets us ask the
question that actually matters for HER2: does this candidate sit on the
trastuzumab epitope, or somewhere else on the domain?
"""
import gzip
import os
import tarfile
import urllib.request

from . import config

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

NON_POLYMER = {"HOH", "WAT", "NAG", "SO4", "GOL", "EDO", "ZN", "MG", "CA",
               "NA", "CL", "PO4", "ACT", "MAN", "BMA", "FUC"}


# --------------------------------------------------------------------- download

def fetch_pdb(pdb_id, dest_dir):
    """Download a PDB entry from RCSB once and keep it."""
    os.makedirs(dest_dir, exist_ok=True)
    path = os.path.join(dest_dir, "%s.pdb" % pdb_id.lower())
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    url = config.RCSB_URL % pdb_id.upper()
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    with open(path, "wb") as fh:
        fh.write(data)
    return path


# ---------------------------------------------------------------------- parsing

class Atom(object):
    __slots__ = ("chain", "resseq", "resname", "name", "x", "y", "z", "element")

    def __init__(self, chain, resseq, resname, name, x, y, z, element):
        self.chain, self.resseq, self.resname, self.name = chain, resseq, resname, name
        self.x, self.y, self.z, self.element = x, y, z, element

    @property
    def residue(self):
        return (self.chain, self.resseq, self.resname)


def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname in NON_POLYMER:
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            try:
                atoms.append(Atom(
                    line[21].strip() or "_",
                    line[22:27].strip(),
                    resname,
                    line[12:16].strip(),
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                    (line[76:78].strip() or line[12:16].strip()[:1]).upper(),
                ))
            except ValueError:
                continue
    return atoms


def parse_cif(path):
    """Minimal mmCIF reader: the _atom_site loop and nothing else."""
    atoms = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i, cols = 0, None
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("loop_"):
            cols, j = [], i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if cols and cols[0].startswith("_atom_site."):
                idx = {name.split(".", 1)[1]: n for n, name in enumerate(cols)}
                i = j
                while i < len(lines):
                    row = lines[i].strip()
                    if not row or row.startswith(("#", "loop_", "_", "data_")):
                        break
                    parts = row.split()
                    if len(parts) >= len(cols):
                        atoms.append(_cif_atom(parts, idx))
                    i += 1
                return [a for a in atoms if a is not None]
            i = j
            continue
        i += 1
    return atoms


def _cif_atom(parts, idx):
    def get(*names):
        for n in names:
            if n in idx:
                return parts[idx[n]]
        return ""
    resname = get("label_comp_id", "auth_comp_id").strip()
    if resname in NON_POLYMER:
        return None
    alt = get("label_alt_id")
    if alt not in (".", "?", "", "A"):
        return None
    try:
        return Atom(
            (get("auth_asym_id", "label_asym_id") or "_").strip(),
            (get("auth_seq_id", "label_seq_id") or "").strip(),
            resname,
            get("label_atom_id", "auth_atom_id").strip().strip('"'),
            float(get("Cartn_x")), float(get("Cartn_y")), float(get("Cartn_z")),
            (get("type_symbol") or "C").upper(),
        )
    except ValueError:
        return None


def load_structure(path):
    return parse_cif(path) if path.lower().endswith((".cif", ".mmcif")) else parse_pdb(path)


def chains(atoms):
    """Chain id -> ordered list of (resseq, one-letter). Built from CA atoms,
    so it is what the model actually resolved, not what was in the construct."""
    out = {}
    for a in atoms:
        if a.name != "CA":
            continue
        one = THREE_TO_ONE.get(a.resname)
        if one is None:
            continue
        out.setdefault(a.chain, []).append((a.resseq, one))
    return out


def sequence_of(atoms, chain):
    return "".join(one for _, one in chains(atoms).get(chain, []))


# ------------------------------------------------------------------- interfaces

def _cell_index(atoms, size):
    grid = {}
    for n, a in enumerate(atoms):
        key = (int(a.x // size), int(a.y // size), int(a.z // size))
        grid.setdefault(key, []).append(n)
    return grid


def interchain_contacts(atoms, chain_a, chain_b, cutoff=4.5):
    """Heavy-atom contacts across two chains.

    Returns (n_atom_pairs, residues_a, residues_b, n_residue_pairs) where the
    residue sets are the interface on each side. Hydrogens are ignored because
    predicted structures do not place them consistently.
    """
    a_atoms = [a for a in atoms if a.chain == chain_a and a.element != "H"]
    b_atoms = [a for a in atoms if a.chain == chain_b and a.element != "H"]
    if not a_atoms or not b_atoms:
        return 0, set(), set(), 0
    grid = _cell_index(b_atoms, cutoff)
    cut2 = cutoff * cutoff
    pairs = 0
    res_a, res_b, res_pairs = set(), set(), set()
    for a in a_atoms:
        cx, cy, cz = int(a.x // cutoff), int(a.y // cutoff), int(a.z // cutoff)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for n in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        b = b_atoms[n]
                        if (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2 <= cut2:
                            pairs += 1
                            res_a.add((a.resseq, a.resname))
                            res_b.add((b.resseq, b.resname))
                            res_pairs.add((a.resseq, b.resseq))
    return pairs, res_a, res_b, len(res_pairs)


def epitope_overlap(candidate_epitope, reference_epitope):
    """Jaccard overlap of two residue sets on the same target chain.

    Used to ask whether a designed binder lands on the trastuzumab epitope.
    Both sets must come from the same numbering, which they do when the target
    chain is the one we supplied.
    """
    a = {r[0] for r in candidate_epitope}
    b = {r[0] for r in reference_epitope}
    if not a or not b:
        return 0.0, 0
    shared = a & b
    return len(shared) / float(len(a | b)), len(shared)


# ----------------------------------------------------------------------- files

def write_fasta(path, records, width=0):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(">%s\n" % header)
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")
    return path


def read_fasta(path):
    records, header, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def make_tar(tar_path, file_paths):
    """Tar with the members at the archive root and no ./ prefix, which is what
    the inverse-folding models expect (the shell equivalent is
    `COPYFILE_DISABLE=1 tar -cf out.tar -C dir .`)."""
    os.makedirs(os.path.dirname(tar_path) or ".", exist_ok=True)
    with tarfile.open(tar_path, "w") as tar:
        for path in file_paths:
            tar.add(path, arcname=os.path.basename(path))
    return tar_path


def find_files(root, suffixes):
    hits = []
    for base, _dirs, names in os.walk(root):
        for name in names:
            if name.lower().endswith(tuple(suffixes)):
                hits.append(os.path.join(base, name))
    return sorted(hits)


# -------------------------------------------------------- antibody chain guess

J_MOTIFS = ("WGQG", "WGKG", "WGRG", "WGAG", "FGQG", "FGGG")


def looks_like_antibody(seq):
    """Rough test for an immunoglobulin V domain.

    Two conserved features: the J motif that closes CDR3, and a cysteine pair
    roughly 60-80 residues apart that forms the domain's internal disulfide.
    This is a heuristic for labelling chains in a downloaded structure, not a
    classifier -- the chain assignment can always be set explicitly.
    """
    if not 90 <= len(seq) <= 260:
        return False
    if not any(m in seq for m in J_MOTIFS):
        return False
    cys = [i for i, c in enumerate(seq) if c == "C"]
    return any(60 <= b - a <= 85 for a in cys for b in cys if b > a)
