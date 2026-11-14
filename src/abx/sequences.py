"""Reference sequences, taken from an experimental structure, not from memory.

Source: PDB 1N8Z, 2.52 Å — trastuzumab (Herceptin) Fab bound to HER2.
    chain A = light chain, chain B = heavy chain, chain C = HER2 extracellular domain.
The sequences below are the SEQRES records, not the ATOM records, so they have no gaps
where the crystal structure is missing residues.

Fetch them again with:  make data
"""

from __future__ import annotations

# ── Trastuzumab, variable domains ───────────────────────────────────────────
# The heavy chain in 1N8Z is the entire Fd fragment; VH is the first part up to
# "...VTVSS", where CH1 begins with "ASTKGPSVFPLAP".
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
    "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)

# The light chain is entirely kappa; VL ends at "...KVEIK", where CL begins with "RTVAAPSVFIFPP".
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
    "RFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

# ── The CDRs, Kabat boundaries ──────────────────────────────────────────────
# Read out of the sequence above, not reproduced from memory: each of them is
# verified with a substring check in tests/test_sequences.py.
TRASTUZUMAB_CDRS = {
    "H1": "DTYIH",
    "H2": "RIYPTNGYTRYADSVKG",
    "H3": "WGGDGFYAMDY",
    "L1": "RASQDVNTAVA",
    "L2": "SASFLYS",
    "L3": "QQHYTTPPT",
}

# CDR-H3 is what the task measures novelty against. 11 residues.
TRASTUZUMAB_CDRH3 = TRASTUZUMAB_CDRS["H3"]

# Anchors around CDR-H3 in VH. They are conserved and let us find the loop in a
# variant without numbering the whole domain.
CDRH3_PREFIX = "YYCSR"     # end of framework 3
CDRH3_SUFFIX = "WGQGT"     # start of framework 4

# ── The target ──────────────────────────────────────────────────────────────
HER2_ECD = (
    "TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQG"
    "YVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPGGLRELQL"
    "RSLTEILKGGVLIQRNPQLCYQDTILWKDIFHKNNQLALTLIDTNRSRACHPCSPMCKGS"
    "RCWGESSEDCQSLTRTVCAGGCARCKGPLPTDCCHEQCAAGCTGPKHSDCLACLHFNHSG"
    "ICELHCPALVTYNTDTFESMPNPEGRYTFGASCVTACPYNYLSTDVGSCTLVCPLHNQEV"
    "TAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAVTSANIQEFAGCKKIFGSLAFLPESFD"
    "GDPASNTAPLQPEQLQVFETLEEITGYLYISAWPDSLPDLSVFQNLQVIRGRILHNGAYS"
    "LTLQGLGISWLGLRSLRELGSGLALIHHNTHLCFVHTVPWDQLFRNPHQALLHTANRPED"
    "ECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLP"
    "CHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGA"
    "CQPCPIN"
)

# The scFv linker everyone uses: (G4S)x3.
SCFV_LINKER = "GGGGSGGGGSGGGGS"


def trastuzumab_scfv(order: str = "VH-VL") -> str:
    """Trastuzumab as an scFv. 'VH-VL' or 'VL-VH'."""
    if order == "VH-VL":
        return TRASTUZUMAB_VH + SCFV_LINKER + TRASTUZUMAB_VL
    if order == "VL-VH":
        return TRASTUZUMAB_VL + SCFV_LINKER + TRASTUZUMAB_VH
    raise ValueError(f"unknown order: {order!r}")

