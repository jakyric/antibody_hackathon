"""Everything that is a choice rather than logic.

Model ids, lanes and first-call times are copied from the CATALOG in
`forge.py` (measured 12 November 2026). If forge.py's catalog changes, this
file is the only place that has to follow.
"""
import os

# ----------------------------------------------------------------- the target
#
# The hackathon target is HER2 (ERBB2). 1N8Z is the trastuzumab Fab bound to
# HER2 extracellular domain IV, so one download gives us three things:
#   * the target sequence,
#   * a positive control (an antibody that is known to bind it),
#   * the epitope, which we compute from the structure instead of quoting
#     residue numbers from memory.
TARGET_PDB = "1n8z"
TARGET_NAME = "HER2"

# 1ZVY chain A is a VHH (nanobody) bound to lysozyme. We use it only as a
# framework donor: AntiFold keeps the framework and redesigns the loops, so the
# CDRs it proposes are new, and the scaffold's original target is irrelevant.
SCAFFOLD_PDB = "1zvy"
SCAFFOLD_CHAIN = "A"

RCSB_URL = "https://files.rcsb.org/download/%s.pdb"

# --------------------------------------------------------------- the platform
#
# lane matters: "inline" models go through `forge.py run --args`, "file" models
# through `forge.py submit --entry`. Getting this wrong is an immediate error,
# so we read it from here rather than remembering it per call site.
MODELS = {
    "msa":        ("msa-search-pdb70", "inline", 40),
    "fold":       ("boltz2", "inline", 80),
    "confirm":    ("openfold3", "inline", 165),
    "redesign":   ("antifold", "file", 110),
    "build_vhh":  ("nanobodybuilder2", "file", 130),
    "build_pair": ("abodybuilder3", "file", 120),
    "backbone":   ("rfd4-proteina", "file", 290),
    "seq_design": ("proteinmpnn", "inline", 100),
    "binder":     ("mosaic", "file", 1310),
    "plausible":  ("esm2-650m", "file", 230),
    "ddg":        ("stab-ddg", "file", 730),
}

# OpenFold2 never got an endpoint up on 12 November and burns 55 minutes before
# it gives up. Named here so nothing in this package reaches for it.
AVOID = {"openfold2"}

# Models with a tested recipe in RECIPES.md. For anything else we read
# `schema <model> --full` and build the call from input_artifact_contract
# rather than guessing an input shape.
RECIPE_KNOWN = {
    "msa-search-pdb70", "boltz2", "openfold3", "proteinmpnn",
    "antifold", "nanobodybuilder2", "abodybuilder3", "rfd4-proteina",
}

MSA_DB = "pdb70_220313"

# ------------------------------------------------------------------ where things go
#
# forge.py writes its receipts to ./nebius-runs relative to its working
# directory, and `resume` finds a run's key by reading them back. So every
# forge call in this package runs with the same cwd, and that cwd is WORK_DIR.
WORK_DIR = os.environ.get("DENOVO_WORK_DIR", os.getcwd())
RUNS_DIR = "nebius-runs"          # forge.py's own, relative to WORK_DIR
INPUTS_DIR = os.path.join(RUNS_DIR, "inputs")
RESULTS_DIR = "results"

# ------------------------------------------------------------------ campaign size
#
# Defaults sized for a hackathon afternoon rather than for completeness:
# generation is cheap, co-folding is what costs time (75 s per candidate, and
# they queue), so we triage hard before spending any of it.
N_GENERATE = 12          # sequences AntiFold is asked for
N_COFOLD = 6             # survivors we co-fold against the target
N_CONFIRM = 2            # top candidates that get a second opinion
SAMPLING_TEMP = "0.20"   # AntiFold: higher is more diverse and less conservative
CDR_REGIONS = "CDR1 CDR2 CDR3"
