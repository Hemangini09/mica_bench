"""Local MICA benchmark configuration."""

from pathlib import Path

ROOT = Path("/home/hemangini/Documents/temp/mica_bench")
DATA = ROOT / "data"
MODELS = ROOT / "models"
LOGS = ROOT / "logs"
RESULTS = ROOT / "results"
USALIGN = ROOT / "USalign"
WORK = ROOT / "work"
PHENIX_ACT = ROOT / "phenix" / "phenix-2.2.1-6174" / "phenix_env.sh"
PULCHRA = ROOT / "src" / "MICA" / "modules" / "pulchra304" / "src" / "pulchra"
MICA_SRC = ROOT / "src" / "MICA"

# pdb_id -> emdb numeric id (no EMD- prefix), resolution Å
# Source: 2026 panel in mica_colab_gpu.ipynb plus README giant-map examples.
PANEL = {
    "9VF5": {"emdb": "53095", "resolution": 1.52},
    "9XDE": {"emdb": "66764", "resolution": 1.64},
    "9TBI": {"emdb": "55772", "resolution": 1.70},
    "9WP6": {"emdb": "66135", "resolution": 1.98},
    "9S3E": {"emdb": "54531", "resolution": 2.10},
    "9RZV": {"emdb": "54414", "resolution": 2.23},
    "9P9L": {"emdb": "71417", "resolution": 2.28},
    "9RVN": {"emdb": "54289", "resolution": 2.30},
    "9S3M": {"emdb": "54540", "resolution": 2.31},
    "9N5N": {"emdb": "48929", "resolution": 2.35},
    "9XON": {"emdb": "67084", "resolution": 2.50},
    "9UIT": {"emdb": "64204", "resolution": 3.14},
    "9W22": {"emdb": "65563", "resolution": 3.12},
    "9UPC": {"emdb": "64391", "resolution": 3.14},
    "33AW": {"emdb": "59318", "resolution": 3.84},
    "9WSD": {"emdb": "66197", "resolution": 4.00},
}

# B2 includes high-res (<2 Å) as a slice, not a separate challenge.
BENCHMARKS = {
    "B2_resolution_decay": {
        "high": ["9VF5", "9XDE", "9TBI", "9WP6"],
        "medium": ["9S3M", "9RVN", "9XON", "9N5N"],
        "hard": ["9W22", "9UIT", "33AW", "9WSD"],
    },
    # EMDB ids from README; confirm maps exist before treating as 2026 held-out.
    "B3_giant": [
        {"pdb": None, "emdb": "26841", "note": "README example ~4834 res"},
        {"pdb": None, "emdb": "33528", "note": "README example ~4770 res"},
    ],
    "B4_bad_af_prior": [],  # filled after AF3 vs native TM-score
    "B5_flexible": [],  # filled after B-factor / local-resolution screen
}

# Full 12-map B2 list kept for ModelAngelo/AF3 scores already on disk.
RUN_PDB_IDS = (
    BENCHMARKS["B2_resolution_decay"]["high"]
    + BENCHMARKS["B2_resolution_decay"]["medium"]
    + BENCHMARKS["B2_resolution_decay"]["hard"]
)

# Fast wave: four maps spanning 1.7–4.0 Å. Giant copies dropped. EModelX(+AF)
# and CryoAtom are queued on GPU behind Cryo2Struct.
FAST_FIRST = (
    "9XON",  # 2.50 Å, 111 residues — shortest
    "9W22",  # 3.12 Å
    "9WSD",  # 4.00 Å
    "9TBI",  # 1.70 Å, 742 residues — run last
)

# Cheng email: MICA vs ModelAngelo + EModelX(+AF) plus newer methods.
# EMProt dropped: Huang lab weight server returned 503. Replacement is
# phenix.map_to_model (already installed, CPU, real de novo baseline).
METHODS = [
    "mica",
    "modelangelo",
    "emodelx",
    "cryoatom",
    "cryo2struct",
    "phenix_map_to_model",
    "af3_alone",
]

CA_MATCH_CUTOFF = 3.0
AUC_DMIN = 0.5
AUC_DMAX = 5.0
AUC_DSTEP = 0.1
