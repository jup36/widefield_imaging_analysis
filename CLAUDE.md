# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MATLAB analysis code for the Buschman Lab's mesoscale widefield calcium imaging pipeline (mouse cortex-wide imaging, often paired with fiber photometry / dual-channel red-calcium + green-dopamine recordings, and behavior in go/no-go or passive tasks). There is no build system, package manager, or test runner — this is a research codebase run interactively in MATLAB (locally) and via batch jobs on Princeton's HPC clusters (Spock, Scotty).

There is no top-level entry point. Work happens by opening a specific script under `Pipelines/` (or a subfolder script) in MATLAB and running/adapting it, or by calling individual analysis functions directly from the command window.

## Running things

- No CLI build/lint/test commands exist. To "run" something, open the relevant `.m` file in MATLAB and execute it (or call the function from the MATLAB command window / another script).
- No automated test suite (no `matlab.unittest` usage in this repo). `Utilities/scotty_smoke_test.m` is a manual, interactive smoke test for the SSH2 connection to the Scotty cluster.
- Path setup is not centralized in a `startup.m`. Scripts add what they need at the top via `addpath(genpath(...))`, e.g. `Pipelines/Example_Pipeline_Local.m` does:
  ```matlab
  addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
  addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
  ```
  Note the dependency on a **sibling repo, `fpCNMF`**, for the convolutional/seqNMF motif-discovery step — it is not part of this repo and must be on the path separately.

## Cross-platform / cross-machine paths

Raw and processed data live on a networked drive mounted differently depending on where MATLAB is running: `Z:\...` on Windows, `/Volumes/buschman/...` on Mac, `/jukebox/buschman/...` on the Spock/Scotty clusters. Path-translation helpers in `Utilities/` are used pervasively any time a path crosses one of these boundaries — expect to see them called in almost every pipeline/cluster-submission script:
- `compatiblepath.m` — auto-converts Mac↔Windows based on current OS (`ispc`/`ismac`).
- `ConvertToBucketPath.m`, `ConvertMacToBucketPath.m`, `ConvertBucketToMacPath.m`, `ConvertBucketToWinPath.m`, `ConvertMacToWinPath.m` — explicit conversions to/from the cluster ("bucket") path form.

When editing or adding pipeline code, follow this existing pattern rather than hardcoding one platform's path.

## Cluster job submission (Spock / Scotty)

Long-running steps (preprocessing a stack, motif fitting, refitting) are dispatched as SLURM jobs on Princeton's clusters rather than run locally:
- `Spock/` and cluster-facing scripts elsewhere (e.g. `Preprocessing/Spock_CombineStacksBVcorrectTrial.m`) hold the code that actually runs *on* the cluster.
- `Spock/WriteBashScript*.m` (Win/Mac/Scotty variants) generate `sbatch` scripts from a template (`Spock/spock_base.sh`) plus a target MATLAB function name and arguments.
- Job submission from a local/Mac MATLAB session to **Spock** is done over SSH using the vendored `Toolboxes/ssh2_v2_m1_r7` toolbox (Ganymed SSH2-over-Java): open a connection with `ssh2_config`, submit with `ssh2_command(conn, 'sbatch ...')`, chain dependent jobs with `sbatch --dependency=afterok:<jobid>`, and always close with `ssh2_close(conn)` — leaving connections open is called out repeatedly in comments as something that annoys PNI IT. **Scotty** submission increasingly goes through `Spock/ssh2_command_scotty.m` instead, a lighter wrapper shelling out to the OS `ssh` executable rather than the Java toolbox.
- "Spock" and "Scotty" are two different Princeton clusters; many pipeline/job scripts exist in parallel `_Spock` and `_Scotty` variants (and further `_MotifVer*`, `_GreenDA`, `_RedCA`, `_L#K#` variants for different hyperparameter sweeps — see below).

## Parameter objects

Pipeline configuration is passed around as MATLAB `classdef` objects rather than plain structs/config files, defined in `ParameterClasses/` (e.g. `general_params.m`, `general_params_dual*.m`, `general_params_mac.m`, `general_params_win.m`, `general_params_asdmodels.m`, `general_params_parietalephys.m`). Convention for using them, seen throughout `Pipelines/`:
```matlab
parameter_class = 'general_params_example';   % pass the CLASS NAME as a string
gp = loadobj(feval(parameter_class));          % instantiate it via feval
```
Functions frequently accept `parameter_class` (a string) rather than the object itself, especially when the function is being shipped off to run on the cluster (where `feval` re-instantiates it remotely). The many `general_params_dual_L{1,3,5,10,15}_K{1,3,5,10,15,20}` variants are per-hyperparameter-sweep configs for CNMF/seqNMF motif fitting (`L` = factor length, `K` = number of factors) — new sweep points get added as new subclass-like variant files rather than by parameterizing one file.

## Directory map

- **`Pipelines/`** — top-level orchestration scripts, one per experiment type/cluster/variant (e.g. `Data_DualPipeline_GNG_Scotty_func_MotifVerGreenDA_L10.m`). `Example_Pipeline_Local.m` and `Example_Data_Pipeline.m` are the annotated, didactic versions to read first — they walk through the full preprocess → hemodynamic correction → dFF → deconvolve/split → motif fit → cluster → refit sequence and explain the cluster job-dependency chaining. Filename suffixes encode the variant: cluster (`Spock`/`Scotty`), task (`GNG`, `passive`), channel (`GreenDA`/`RedCA`), and sweep params (`L#`, `K#`).
- **`Preprocessing/`** — turns raw `.tif` stacks into ΔF/F: alignment/registration (`ManualAlignment.m`, `RegisterReferenceImages.m`), vasculature masking (`MaskVasculature*.m`), hemodynamic correction (`HemodynamicCorrection*.m` — several numbered/lettered variants, e.g. `FF0`–`FF3`), `makeDFF.m`, then `ProcessAndSplitData*.m` (deconvolution + train/test chunking, dispatched to the cluster).
- **`Deconvolution/`** — holds only `out/` (accumulated Slurm `.out` job logs, e.g. `DeconvMethods*.out`); no source lives here. Actual deconvolution code (`lucric.m`, `lucyrichardson.m`) lives in `Preprocessing/`.
- **`Postprocess/`** — post-hoc analysis of processed dF/F, operating on the `tbytDat` struct (see "Core data structures" below): trial-aligned averaging (`dffPostprocess*_averaging*.m`, split heavily by brain region — `_m1`, `_m2`, `_rs`, `_ss`, `_v1`/`_vq`), cross-session/animal collection (`*_collect.m`), SVM/classifier analyses (`Postprocess/SVM/`, `trainDffClassifier.m`).
- **`TDR/`** — Targeted Dimensionality Reduction / GLM-based trial-type decoding on motif ("H") time series. `TDR/TDR_workflow.m` and the `glm_predict_H*.m` / `glm_predict_greenDA_H*.m` / `glm_predict_redCalcium_H*.m` family (design matrix from `tbytDat`/`tbytDat_hAligned` via `events_to_design.m` + raised-cosine bases from `make_rcos_basis_ortho.m`, then ridge/GLM fit of stacked motif activity `H`) are the fit stage — their output is saved as a **`glmRez`** struct (see below).
  - **`TDR/glmTDR/`** — the active area of development (per current git status), consuming `glmRez` structs produced by the fit stage above. Builds "anchor axes" from GLM coefficients per predictor group and projects motif trajectories onto them, including leave-one-out (LOO) per-mouse variants to guard against circularity. Key chain: `find_latest_glmRez_file.m` (locate saved fit results by `<header>_<keyword>_<MMDDYY>.mat` naming) → `targetedDimRed_from_glmRez.m` (build axes from a `glmRez` struct: predictor-group SVD → optional Gram-Schmidt orthogonalization via `gs_orth_rows.m` → sign-fixing via `sign_fix_axis.m` → projection) → `buildAnchorAxes_from_glmRezList.m` / `buildPerMouseLOOAnchors_fromExperts.m` (aggregate axes across sessions/mice) → `projectGlmRezC_to*Anchors.m` (project new data onto fixed anchors) → `computePerMouseLOOProjectionAndStats.m` / `validatePerMouseLOO_noCircularity.m` (stats + LOO sanity checks) → `visualize_glmTDR_rez*.m` (plotting). The `glmRez` struct convention (documented in `targetedDimRed_from_glmRez.m`'s header comment) is: `beta`, `X_design`, `muX`/`sdX`, `Yz`, `group` (struct array with `.name`/`.cols`), `decBins.time` — read that docstring before touching this subfolder. The `_redCal_L10K10` suffixed files are the current focus of active development (red-calcium channel, CNMF `L=10,K=10`).
  - **`TDR/pca/`, `TDR/motifTransition/`, `TDR/tensor_decomposition/`** — related dimensionality-reduction/analysis variants on motif data.
- **`Analysis_functions/`, `BehavioralAnalysisFunctions/`** — general analysis (entropy, residuals) and behavior-side analysis (DeepLabCut import/parsing via `parse_dlc.m`/`CompileDLCData.m`, GMM-based behavioral state classification).
- **`Plotting/`** — shared plotting helpers (motif visualization/GIFs, montages, custom colormaps).
- **`GutcheckFigures/`** — sanity-check scripts comparing pipeline choices (deconvolution, PCA denoising, hemodynamic correction, smoothing) against alternatives — useful reference for *why* a default parameter was chosen.
- **`Spock/`** — cluster-side scripts and SSH/job-submission helpers (see above).
- **`Local/`** — small local-machine-only variants of cluster scripts (e.g. `CombineStacksBVcorrect_Local.m`).
- **`Sync/`** — stimulus/timing sync code (currently only an unsaved `.asv` is present, no committed `.m`).
- **`ParameterClasses/`** — pipeline configuration classes (see above).
- **`Utilities/`** — grab-bag of shared helpers: path conversion (see above), file discovery (`GrabFiles*.m`, `findFilePattern.m`), statistics (`PermutationStatistic.m`, `SlidingWindowPEV.m`), video/GIF generation for trial-aligned face/widefield data (`make_dff*_videos*.m`), motif statistics (`motif_statistics*.m`).
- **`Toolboxes/`** — vendored third-party code (ExportFig, distributionPlot, natsortfiles, an SSH2-over-Java client for cluster job submission). Treat as read-only external dependencies; don't refactor into house style.
- **`Sandbox/`** — exploratory/one-off scripts, not part of any maintained pipeline.
- **`CommunicationSubspace/`, `Network Graphs/`** — currently contain only MATLAB autosave (`.asv`) files, not committed `.m` source; treat as inactive/abandoned unless the user says otherwise.

## Core data structures

These names recur across `Postprocess/` and `TDR/` and are worth recognizing on sight:
- **`tbytDat`** — the central trial-by-trial struct array for a session (loaded from `Matfiles/tbytDat.mat`), with fields like `.evtType`, `.evtOn`/`.evtOff`, `.periEvtWin`, `.rewardTrI`/`.punishTrI`.
- **`tbytDat_hAligned`** — the same trial structure with motif-activity (`H`) aligned per trial; this is what feeds the GLM fit stage in `TDR/`.
- **`trI`** — a trial-index/trial-type struct (logical masks per condition, e.g. `.goI`), produced by `trialTypeInfoAuditoryGng.m` / `trialTypeInfoAuditoryGngTbytDat.m`.
- **`header`** — a session identifier of the form `<mouseID>_<sessionDate>` (e.g. `'m1045_122424'`), produced by `extract_date_animalID_header.m` and used as the filename prefix for nearly all saved results (see the `<header>_<descriptor>_<MMDDYY>.mat` convention below).
- **`glmRez`** — see `TDR/glmTDR/` above.

## Conventions worth knowing

- **Naming encodes experiment variant, not just function.** Many files are near-duplicates differing only by cluster (`_Spock`/`_Scotty`), channel (`_GreenDA`/`_RedCA`), brain region (`_m1`/`_m2`/`_rs`/`_ss`/`_v1`), or hyperparameters (`_L{n}K{n}`). When fixing a bug, check whether sibling variants need the same fix — this codebase does not centralize shared logic across variants.
- **Saved result files follow a `<header>_<descriptor>_<MMDDYY>.mat` naming pattern** (e.g. `m1045_122424_glmRez_redCal_L10K10_070126.mat`), where `header` is typically `<mouseID>_<sessionDate>`. `find_latest_glmRez_file.m` is the canonical example of parsing this convention to pick the most recent fit.
- Motif discovery uses CNMF/seqNMF terminology throughout: `K` (number of motifs/factors), `L` (motif length in timebins), `W`/`H` (motif shapes / temporal loadings). `SetAnalysisOptions.m` documents the default CNMF options.
- `.asv` files (MATLAB autosave) are committed in several places in this repo — don't treat their presence as accidental, but also don't treat them as authoritative source (they're often stale copies alongside a real `.m` file, or the only remnant of an abandoned script).
