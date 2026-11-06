







# MICA cryo-EM model-building benchmark (2026)

Local benchmark of **MICA** vs ModelAngelo, EModelX(+AF), CryoAtom, Cryo2Struct, **phenix.map_to_model**, and an AF3-alone control on **12 × 2026 EMDB/PDB maps** (resolution 1.52–4.0 Å).

**Read [`DOCUMENTATION.md`](DOCUMENTATION.md) for the full decision log** (why EMProt was dropped, why this panel, AF3 Server, Phenix, GPU/CPU split, failures kept in the table).

## Methods

| Method              | Role                                           | How it runs here                                   |
| ------------------- | ---------------------------------------------- | -------------------------------------------------- |
| MICA                | Method under test                              | AF3 Server → Phenix `dock_in_map` → `run.py` |
| ModelAngelo         | Paper baseline                                 | GPU, 12/12 finished                                |
| EModelX(+AF)        | Paper baseline                                 | GPU, after MICA                                    |
| CryoAtom            | Newer DL                                       | GPU, after EModelX                                 |
| Cryo2Struct         | Cheng de novo control                          | GPU, overlapping CPU dock                          |
| phenix.map_to_model | **Replaces EMProt** (weights server 503) | CPU, `7_run_phenix_map_to_model.py`              |
| AF3-alone           | Prior-only control                             | Server CIFs; F1=0 until docked                     |

EMProt is cloned under `src/EMProt` but **not run**.

## Fast wave (what is actually running)

Four small maps only: **9XON (2.50 Å), 9TBI (1.70 Å), 9W22 (3.12 Å), 9WSD (4.00 Å)**.

Dropped as too slow: giant B3, B4/B5, 9VF5/9XDE/9UIT and other large copies, Phenix-docked MICA, EModelX, CryoAtom. ModelAngelo + AF3-alone scores for the original 12 maps stay on disk.

## Detached jobs

```bash
bash run_detached.sh
```

Logs: `logs/*_detached.log`. Scores: `python3 8_evaluate_metrics.py` → `results/scores_master.csv`.

Hardware: RTX 6000 Ada 48 GB, 128 cores. One GPU — do not start two DL builders at once.
