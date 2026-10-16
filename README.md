# Branch-aware lookup tables for selective harmonic mitigation: code and data

This package reproduces the results of

> M. S. Sukkar and K. K. Abdalla, "Why Interpolated Lookup Tables Fail for Selective Harmonic Mitigation: Branch Decomposition and Branch-Aware Tables in Three-Level Neutral-Point-Clamped Converters," *International Journal of Intelligent Engineering and Systems* (under review).

It contains the MATLAB code that computes every result of the paper, a Python script that draws its figures, the Simulink models of the circuit simulations, and the complete results of the run the paper reports. Every switching-angle solution is computed by MATLAB's `fmincon` with analytic gradients. The neural networks, Gaussian processes, ANFIS models and the workbook writer are written in plain MATLAB, so no other toolbox is needed.

## Requirements

| Purpose | Requirement |
|---|---|
| All computations | MATLAB R2016b or newer (tested with R2024b), with the Optimization Toolbox |
| Circuit simulations (Sections 6.4 to 6.6) | Simulink and Simscape Electrical (Specialized Power Systems), R2024b or newer: the models are saved in R2024b |
| Figures 2 to 7 of the paper | Python 3 with numpy, scipy and matplotlib |

The published results were computed with MATLAB R2024b on Windows (64-bit).

## Quick start

In MATLAB, from the package folder:

```matlab
check_setup           % one-minute check of MATLAB, fmincon, the model and file output
reproduce_paper       % figures, tables and numbers from the stored stage results (a few minutes)
verify_reproduction   % compares output/ with results/
```

Every run writes to `output/` (or the folder given); `results/` is never changed.

To recompute everything from scratch:

```matlab
reproduce_paper('all')              % every stage, into output_all/ (about 8 hours)
verify_reproduction('output_all')
```

The figures of the paper are drawn from a results file by

```
python make_figures.py                            # from results/paper_results.mat
python make_figures.py output/paper_results.mat   # from your own run
```

and written to `output/figures/manuscript/`.

## Contents

```
check_setup.m             pre-flight check and run-time estimate
reproduce_paper.m         reproduces the MATLAB part (stages, figures, tables, numbers)
reproduce_circuit.m       reruns the Simulink simulations
verify_reproduction.m     compares a run with the published results, value by value
make_figures.py           draws Figs. 2-7 of the paper from a results file
setup_paths.m             adds code/ and models/ to the MATLAB path
code/
  run_paper.m             the stage runner: cache, random streams, outputs
  core/                   waveform model and solver, multistart, continuation, branch set
  converter/              exact R-L currents, gate-pattern check, device data,
                          switching-function model with non-ideal switching
  learners/               MLP, tanh network (Levenberg-Marquardt), GP/RBF, ANFIS, classifier
  experiments/            one function per stage (shm_x_*)
  circuit/                Simulink drivers: cases, runs, THD evaluation
  output/                 figures, tables, number files, CSV/JSON export
models/
  test.slx                three-level NPC converter, ideal switches, 64 ohm + 180 mH per phase
  test_80ohm_25mH.slx     the same with 80 ohm + 25 mH per phase (second load of Section 6.5)
  NPC_SHM.slx             the same with datasheet IGBT and diode models (Section 6.6)
results/                  the results the paper reports (see below)
```

Every function has a help text: `help shm_branchset`, `help run_paper`.

## How the computation is organised

The study is split into stages. Each stage is one function of `code/experiments` (or `code/core`), draws its random numbers from its own fixed stream, `rng(1000 + k)`, and stores its result in `<output>/cache/<stage>.mat`. A stage whose result is in the cache is loaded instead of recomputed, so an interrupted run continues where it stopped.

| Stage | Computes | Paper | Time* |
|---|---|---|---|
| D | branch set: multistart at 11 anchors, continuation over 201 grid points, probes, fold-end searches, stored subset | Section 3, Figs. 2-3 | 29 min |
| S, V | stored branches extended over (Ma, delta) and (Ma, delta, theta_gap) | Sections 4.4, 5 | 9 min |
| E1, E1b | multiplicity of optima; interior ends of the owning branches | Section 3 | 1 min |
| E2 | lookup tables on the grid, continuity-constrained path | Fig. 2 | <1 min |
| E2b | table comparison at the 200 midpoints | Table 1, Figs. 4-5 | 20 min |
| E3, E3b, E3c | regression inside a branch, d = 1, 2, 3 | Table 2 | 14 min |
| E4 | selection of the owning branch, d = 1, 2, 3 | Section 4.4 | 18 min |
| E5 | storage, scaling and on-line time of the deployment options | (supplementary) | <1 min |
| E6 | dense multistart, spectrum above the band, per-order limits, upper angle bound | Sections 2.2, 6.2 | 11 min |
| E7 | gate pattern and dead time | Section 6.3 | <1 min |
| E9 | tables under non-ideal switching (switching-function model) | Table 4 | 18 min |
| IND | independent search: interior-point from Halton starts | Section 6.1 | 124 min |
| CMP | per-order compliance | Section 6.2 | <1 min |
| STO | Table 1 at equal storage, with the published baselines | Table 1, Section 4.3 | 47 min |
| RGS | regression statistics: seeds, model selection, ranking | Table 2, Section 5 | 174 min |
| BSW | branch switching under R-L loads, exact | Table 3, Section 6.5 | <1 min |
| TIM | operations and time of a table read | Section 6.5 | <1 min |

\* On the reference PC (MATLAB R2024b, Windows). The whole set takes about 8 hours.

`reproduce_paper('stages', {'IND'})` recomputes one stage, and every stage that depends on it, and loads the others from `results/cache`. `reproduce_paper('quick')` runs the whole chain on small grids (a test of the installation, not the paper's numbers).

## Circuit simulations

```matlab
reproduce_paper       % first: the angle sets come from output/paper_results.mat
reproduce_circuit     % about one hour; then verify_reproduction
```

`reproduce_circuit` runs three sets, each of which can be run alone with `reproduce_circuit('runs', {...})`:

| Set | Model | Runs | Paper |
|---|---|---|---|
| `'ideal'` | `test` | 9 steady-state cases, the branch switches of Fig. 7, the switch at all 11 switch points | Section 6.4, Fig. 6, Section 6.5, Fig. 7 |
| `'datasheet'` | `NPC_SHM` | the 9 cases under the 6 conditions of Table 4 | Section 6.6, Table 4 |
| `'second_load'` | `test_80ohm_25mH` | as `'ideal'`, with a load of 80 ohm + 25 mH | Section 6.5 |

The models are three-phase three-level NPC converters with an 800-V DC link split by two 470-uF capacitors. The angles are written to the base workspace by `shm_sim_run`; the line-to-line voltage `vab`, the phase current `ia` and the neutral-point voltage `dVnp` are logged. `NPC_SHM` reads its device, dead-time and load parameters from the base workspace, where `shm_sim_nonideal` sets them from `shm_device` and `shm_nonideal_conditions`.

## Results

`results/` holds the results the paper reports, in the layout every run writes:

| File | Content |
|---|---|
| `numbers.tex`, `numbers_ijies.txt`, `numbers_revision.txt`, `numbers_supplement.txt` | every number of the paper, one `key = value % description` per line (`numbers.tex`: LaTeX macros) |
| `summary.txt`, `summary_revision.txt` | the same results as readable lists and tables |
| `tables/*.tex` | result tables in LaTeX |
| `figures/` | an extended set of figures drawn by MATLAB from the same results |
| `figures/manuscript/` | Figs. 2-7 of the paper, drawn by `make_figures.py` |
| `artifact/` | the main results as CSV, and every setting as JSON (see below) |
| `paper_results.mat`, `paper_results.xlsx` | all results as one struct `R`, and as a workbook |
| `cache/` | the result of every stage |
| `simulink/` | the circuit results |

`artifact/` can be read without MATLAB:

| File | Content |
|---|---|
| `branches.csv` | every traced branch: Ma, THD49, the 15 angles |
| `branch_info.csv` | per branch: range, end types, owned share of the grid, stored or not |
| `envelope.csv` | per grid point: the lowest optimum found, its owner, the stored-subset envelope |
| `owner_changes.csv` | every owner change and the switch point the branch-aware table uses |
| `midpoints.csv` | the 200 evaluation midpoints of Table 1 |
| `table1_R100.csv` | angles of the optimum, the plain table and the branch-aware table at the midpoints |
| `config.json` | the settings of every stage, the random seeds, the solver tolerances and acceptance rule, the subset-selection rule, the per-order limits |

### Where each result of the paper is

| Paper | Stages | Keys and files |
|---|---|---|
| Section 3, Figs. 2-3 | D, E1, E1b | `branch.*`, `distinct.*`, `ownerchanges`, `stored.*` in `numbers.tex`; `d.*`, `e1b.*` in `numbers_ijies.txt`; `artifact/` |
| Section 4, Table 1, Figs. 4-5 | E2, E2b, STO | `lut.*`, `split.*`, `cont.*`, `nearest.*` in `numbers.tex`; `e2b.*`; `rev.sto.*` |
| Section 4.4 | E4 | `rev.sel.*` |
| Section 5, Table 2 | E3, E3b, E3c, RGS | `reg1.*`, `reg2.*`, `reg3.*`; `rev.rgs.*` |
| Section 2.2 | E6 | `bound.*` |
| Section 6.1 | IND | `rev.ind.*` |
| Section 6.2 | CMP, E6 | `rev.cmp.*`; `spec.*`, `comply.*`; `numbers_supplement.txt` |
| Section 6.3 | E7 | `sw.*` |
| Section 6.4, Fig. 6 | circuit `'ideal'` | `sim.*` |
| Section 6.5, Table 3, Fig. 7 | BSW, TIM, circuit | `rev.bsw.*`, `rev.tim.*`, `rev.simsw.*`; `simulink/switches_two_loads.txt` |
| Section 6.6, Table 4 | E9, circuit `'datasheet'` | `ni.*`, `sim.ni.*`; `tables/tab9_nonideal.tex` |

The operation count of a pchip read with on-line slopes (Section 6.5, about 290 multiplications and 190 additions) is counted by hand; stage TIM counts the read with precomputed coefficients. Fig. 1 is a circuit diagram.

## Reproducibility

- Every stage draws its random numbers from a fixed stream, so a rerun with the same MATLAB release on the same computer reproduces the published results bit for bit. `verify_reproduction` checks this for every stage result, circuit result and text output.
- Wall-clock times (the run times, the on-line times of stages E5 and TIM, the time per corrector solve) are measurements of the computer and differ from run to run, as do the dates and the platform recorded with each run. `verify_reproduction` lists them separately and does not count them as differences.
- With another MATLAB release or processor, `fmincon` may differ in the last bits of its iterates. Where many local optima lie close together this can change which optimum a random start reaches, so individual numbers can differ slightly; the stored stage results in `results/cache` are the reference.
- Workbooks (`.xlsx`) and figures (`.pdf`, `.png`) carry creation dates and are not compared; they are drawn from the compared data.

## Contact

Muntadher Saeed Sukkar, Department of Electrical Engineering, College of Engineering, University of Babylon, Hilla, Iraq (muntazer.moaawin459@student.uobabylon.edu.iq).
