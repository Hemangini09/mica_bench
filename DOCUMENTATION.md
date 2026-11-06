# MICA cryo-EM model-building benchmark (2026)

## 1. What the assignment actually is

Read the MICA paper (*Commun. Chem.*).

1. Benchmark MICA vs other recent DL cryo-EM model builders on **2026 EMDB/PDB** maps.
3. **Must include** the two methods from the MICA paper: **ModelAngelo** and **EModelX(+AF)**.
4. Add newer methods.
5. Write a report: methods, datasets, procedure, metrics, results, conclusions.

It did **not** require five named “scientific challenges” as separate products. An early local README invented B1–B5 (near-atomic / resolution decay / giant maps / bad AF prior / flexible loops) as a proposal. I kept **resolution-stratified 2026 maps (B2)** as the first wave because that is what can be grade without waiting for giant-map downloads and AF-prior tagging.

---

## 2. Panel (wave 1)

Twelve deposited 2026 entries, unique protein sequences only in FASTA (copies collapsed so AF3 Server stayed under token/job limits).

| Bin    | PDB                    | EMDB                       | Å         |
| ------ | ---------------------- | -------------------------- | ---------- |
| High   | 9VF5, 9XDE, 9TBI, 9WP6 | 53095, 66764, 55772, 66135 | 1.52–1.98 |
| Medium | 9S3M, 9RVN, 9XON, 9N5N | 54540, 54289, 67084, 48929 | 2.30–2.50 |
| Hard   | 9W22, 9UIT, 33AW, 9WSD | 65563, 64204, 59318, 66197 | 3.12–4.00 |

Also listed in `config.py` but **not run** this wave: 9S3E, 9RZV, 9P9L, 9UPC.

**Deferred**

- **B3 giant complexes** (README examples EMD-26841 / 33528): maps are huge; one GPU-night is not enough until B2 MICA exists.
- **B4 bad AF prior:** fill after AF3-vs-native TM; several AF3-alone TMs are already low (e.g. 9XON 0.04, 9WP6 0.10) and can seed B4 later.
- **B5 flexible loops:** needs local-resolution / B-factor screen; not blocking the  table.

Maps, FASTA, and deposited CIFs: `1_download_datasets.py` → `data/<PDB>/`.

---

## 3. Methods: chosen, dropped, replaced

| Method                        | Role                                | Status                                                                  | Why this decision                                                                                                                                                                                                             |
| ----------------------------- | ------------------------------------ | ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **MICA**                | The method under test                | Env + weights + Pulchra ready; official run waits on Phenix dock of AF3 | Required. Official protocol is AF3 → `phenix.dock_in_map` → FPN. We did **not** use the paper’s ancient conda (PyTorch 1.8); that stack is incompatible with this CUDA 13.2 box. Env is Python 3.10 + torch cu128. |
| **ModelAngelo**         | Paper baseline                       | **12/12 done** (~9.5 h)                                           | Required. ESM-1b via `TORCH_HOME`.                                                                                                                                                                                          |
| **EModelX(+AF)**        | Paper baseline                       | Code cloned, not run                                                    | Required. Needs docked AF / similar map-frame models; queued after MICA on the GPU.                                                                                                                                           |
| **CryoAtom**            | Newer DL (NSMB)                      | Cloned, not run                                                         | Newer than the paper. GPU queue after MICA/EModelX.                                                                                                                                                                           |
| **Cryo2Struct**         | -lab de novo (not in the MICA paper) | Running detached; 9VF5/9XDE previously failed                           | Included because it is this lab’s method and a true de novo control vs AF3-conditioned MICA. Failures stay in the table.                                                                                                     |
| **AF3-alone**           | Control                              | 16 Server jobs unpacked and scored                                      | Isolates “AF3 prior vs map-conditioned building”. Server models are **not in the map frame**, so Cα F1 is 0 until dock; US-align TM still measures fold.                                                             |
| **phenix.map_to_model** | **Replaces EMProt**            | Runner added; detached CPU job                                          | See §4.                                                                                                                                                                                                                      |
| **EMProt**              | Was a “newer DL” candidate         | **Dropped**                                                       | Huang lab weight HTTP returned **503 / timeout**. Cloning the repo without weights is theatre. We do not pretend it ran.                                                                                                |

**Not added:** more GitHub cryo-EM DL builders (DeepMainmast, CryoREAD, …). One GPU, and  asked for MICA vs the paper pair plus a few extras, not a zoo.

**Hugging Face:** no drop-in cryo-EM builders; only ESM-1b (already used inside ModelAngelo).

---

## 4. Why EMProt → `phenix.map_to_model`

EMProt could not download weights. Options were: wait on Huang lab, clone another DL repo, or use software already on this machine.

**`phenix.map_to_model`** (Phenix 2.2.1-6174, Linux x86-64 installer, no Coot, no CUDA build):

- Already installed for MICA docking; academic non-profit license; perpetual for this version.
- Real, citable de novo map interpretation (Terwilliger et al.), not a stub.
- **CPU only** (Phenix GPU is for `aquaref`, which we do not use), so it runs **in parallel** with Cryo2Struct/MICA on the GPU.
- Uses the same FASTA + MRC + deposited resolution as everyone else.
- Gives  a classical baseline next to DL, which the email did not forbid and the report should have.

We did **not** pick the CUDA Phenix installer or Coot: unused for `dock_in_map` / `map_to_model` / `real_space_refine`.

---

## 5. AlphaFold3

Local AF3 weights/AlphaFold-3 install were not available. **AlphaFold Server** was used (jobs JSON in `af3_server/`). Unique chains only, Copies=1, 16 jobs, zip `folds_2026_10_25_08_09.zip` unpacked to `data/<PDB>/AF3/`.

Staging: `2_stage_af3_and_score.py` writes CA-only combined PDBs for AF3-alone scoring. MICA needs Server-style `AF3_results/*/*_model_0.cif`; we symlink unique jobs to FASTA ids (`9VF5_1`, …) under `work/<PDB>/AF3_results/`.

**TOS / science:** Server models must be docked before they are a map-space model. Docking is `utils/dock_in_map.py` (Merizo domains, then `phenix.dock_in_map` nproc=48). Contours from EMDB API (`work/contours.json`), else map 90th percentile.

---

## 6. Metrics

Same engine for every method: `8_evaluate_metrics.py` → `results/scores_master.csv`.

- Map-frame Cα match at 3 Å: precision, recall, F1, Cα RMSD, sequence match.
- AUC-PR sweeping 0.5–5.0 Å.
- US-align (`mica_bench/USalign`, compiled locally): TM-score, aligned length, seq identity.

Failed jobs are **rows with `status=missing` or fail logs**, not deleted.

Known scoring quirks (keep, then footnote in the report):

- ModelAngelo **9TBI** high TM, very low F1: raw CIF has huge extra CA count (`9TBI_raw.cif`). Prefer pruned CIF for the table.
- **9RVN** ModelAngelo TM high, F1 0: likely frame/which-CIF issue; inspect before the report.
- AF3-alone F1 = 0 everywhere until dock: expected.

---

## 7. Engineering decisions

| Decision                                            | Why                                                                                                     |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| One method per GPU at a time                        | 48 GB; MICA/ModelAngelo/Cryo2Struct each want tens of GB.                                               |
| CPU Phenix dock ∥ GPU Cryo2Struct                  | Dock does not use the GPU.                                                                              |
| Detached `setsid`/`nohup` (`run_detached.sh`) | Chat-tied shells died overnight; zero docked PDBs.                                                      |
| `bash -c` not `bash -lc`                        | Login shells source `.profile`, which errors on a missing snap `env` and can abort with `set -e`. |
| Cryo2Struct resample origin `.x/.y/.z`            | Assigning a (3,) vector into MRC `header.origin` crashed the job.                                     |
| MICA env not paper `env.yml`                      | PyTorch 1.8 + this driver is a dead end.                                                                |
| Merizo on CPU (`-d cpu`)                          | Do not steal the GPU from Cryo2Struct. Needed `matplotlib`.                                           |
| Domainwise dock, chainwise fallback                 | Paper default is domains; empty `AF3_domains` used to fail dock.                                      |
| Keep EMProt git clone                               | Provenance; we still do not run it.                                                                     |
| Phenix batch prefix `.../phenix-2.2.1-6174`       | Vendor requirement so the env is isolated.                                                              |

---

## 8. Pipeline files (what exists now)

| Script                           | Function                                   |
| -------------------------------- | ------------------------------------------ |
| `config.py`                    | Panel, methods, paths, cutoffs             |
| `1_download_datasets.py`       | MRC / FASTA / CIF                          |
| `2_stage_af3_and_score.py`     | Combine AF3 CA PDBs + score                |
| `3_run_modelangelo.py`         | ModelAngelo GPU                            |
| `4_run_cryo2struct.py`         | 1 Å resample + Cryo2Struct                |
| `5_prepare_and_dock.py`        | Stage MICA inputs, Merizo, Phenix dock     |
| `6_run_mica.py`                | MICA `AF3_struct` after docks + GPU free |
| `7_run_phenix_map_to_model.py` | EMProt replacement                         |
| `8_evaluate_metrics.py`        | Master CSV                                 |
| `run_detached.sh`              | Survive hangup                             |
| `af3_server/`                  | Server JSON/FASTAs                         |

Code/weights live under `src/` and `weights/`. MICA input copies under `work/<PDB>/`. Outputs under `models/<method>/<PDB>/`. Logs under `logs/`.

---

## 8b. Fast-wave cut (25 Oct evening)

Dropped anything that would take overnight-plus: B3 giant maps, B4/B5, 9VF5/9XDE/9UIT/9RVN and other large targets, official `phenix.dock_in_map` + MICA refine, EModelX, CryoAtom (no env).

**Keep now:** 9XON, 9TBI, 9W22, 9WSD. Methods: ModelAngelo (already 12/12), AF3-alone (already 12/12), Cryo2Struct, MICA `AF3_struct_free`, `map_to_model` on this short list only.

## 9. Live status at last documentation write (25 Oct 2026 evening)

- ModelAngelo 12/12 scored.
- AF3-alone 12/12 scored (fold TM only).
- Phenix installed; 9VF5 `dock_in_map` running (~15+ min, no PDB yet).
- Cryo2Struct on 9VF5 (EMD-53095); earlier 9VF5/9XDE failed after GPU nets (no final PDB).
- MICA waiting on first dock + GPU.
- EModelX / CryoAtom not started.
- `phenix.map_to_model` added and launched on CPU (`nproc=16` so it shares cores with dock’s 48).

Wall clock if the machine stays on: docks 8–16 h; MICA 8–16 h after GPU free; map_to_model overlapping on CPU (hours per map); EModelX then CryoAtom another GPU day. Honest full table: **~1.5–2 days**, not the old README’s 2.5 hours.

---

## 10. Report caveats to keep

- AF3 Server models vs deposited coordinates: TM can be high with F1 = 0.
- Cryo2Struct failures are data, not to be silently dropped.
- MICA official numbers need docked AF3; `AF3_struct_free` is a weaker protocol and should be labelled if used.
- Phenix academic license is non-profit R&D; do not redistribute the installer tarball.
- 9TBI/9RVN ModelAngelo F1 anomalies need a one-line methods note after we lock which CIF is scored.
