# pgenfit3 — total decay-curve fit (no neutron gates)

Unbinned, extended maximum-likelihood fit of the beta-decay curve of an
implanted parent and its whole daughter network, done **simultaneously** in
forward and backward time:

| category | data | model |
|---|---|---|
| `pos` | `startTime <= x <= T` | Σ<sub>m</sub> C<sub>m</sub> e<sup>−λ<sub>m</sub>t</sup> + `bkg` |
| `neg` | `−T <= x < 0` | `bkg` |

`x = t_beta − t_ion`. A beta cannot precede its own implant, so the
backward region holds only accidental (wrong-ion) correlations. Its rate
`bkg` (counts/s) is the same parameter as the forward-time background, so the
background is measured and propagated within one likelihood. pgenfit's
`mainT12fit` (unbinfit `fitOpt=2`) instead fixed `nbkg` from the backward
count before the fit. The neutron multiplicity `y` is not used.

C<sub>m</sub> is the sum of the Bateman terms of every decay path containing
e<sup>−λ<sub>m</sub>t</sup>. The parent term is `N0raw*be`. P1n/P2n and isomer ratios
enter only through how the daughters are fed.

## Files

| file | role |
|---|---|
| `main.cc` | reads a parmsex file, writes `path.txt` and `decayModel_cal.cc` (the model) |
| `decayModel.hh` | interface of the generated model |
| `fit_decay_curve.cxx` | fit driver: data, fit, goodness of fit, plots, `PULLDATA` line |
| `decaypath.cc/.hh`, `common.hh` | decay network (from pgenfit2; isomer error columns fixed) |
| `simulation.cc/.hh`, `mainsimulation.cc` | toy generator (from pgenfit2, with the fixes below) |
| `build.sh` | builds everything for one decay network into a directory |
| `run_single_toy.sh`, `run_pull_study.sh`, `run_real_data.sh` | drivers |
| `plot_pulls.py` | pull summary and figure from a pull study |
| `mem_guard.sh`, `env.sh` | memory limits; ROOT version check |
| `examples/` | isomer and alpha test chains, alpha-chain input example, beam settings |

## Requirements

ROOT ≥ 6.32. If the ROOT on `PATH` is older:

```bash
export PGENFIT3_ROOT=/data01/userdata/ROOT/root-6.40.02_patch
```

## Usage

```bash
# real data (ROOT file with "tree": x, optional ionT)
./run_real_data.sh data.root parms.txt [effparms|none] [ncpu=8] [boundaryMarginSec=0] \
    [out_prefix] [startTime=0] [timeRange=10] [linBinFactor=4] [noplot]

# one simulated toy
./run_single_toy.sh parms.txt simparmsex.txt <seed> [effparms|none] [boundaryMarginSec=900] \
    [ncpu=4] [out_prefix=toy] [startTime=0] [timeRange=10]

# pull study
MEM_PER_WORKER_GB=14 ./run_pull_study.sh parms.txt simparmsex.txt <startSeed> <nToys> \
    [effparms|none] [nWorkers=4] [boundaryMarginSec=900] [ncpu_per_fit=1] [prefix=pull] \
    [startTime=0] [timeRange=10]
./plot_pulls.py pull_study_results/<prefix>_all.txt
```

Every run builds in its own `builds/<tmp>` directory (deleted at the end),
so several runs can run at the same time without corrupting each other's
build. pgenfit2 built in shared directories and did not stop on errors.

### Parameters

The parmsex format is the same as pgenfit/pgenfit2; a negative value means float.

- half-lives of any species: float if negative
- isomer ratio: floats the ground-state feeding `py<gs>` if negative
- P1n/P2n: **held fixed by default** (they only shape daughter feeding here);
  `FLOAT_PN=1` honours their float flags
- `be` (effparms column 1): negative = free, positive with error = Gaussian
  constraint, else fixed; `none` = 1 fixed
- always free: `N0raw`, `bkg`

`SCAN_PARAM=l0` (and `SCAN_POINTS/SCAN_MIN/SCAN_MAX`) prints a
profile-likelihood scan as `SCANDATA` lines.

### Alpha decay

This follows pgenfit's alpha extension (`Pb218parms.txt` format).

- **parmsex columns 24–28**: `AlphaBR errLo errHi lower upper`, with AlphaBR in
  percent of all decays and a negative value meaning float. Rows without alpha
  may give 0 or omit the block. On an isomer row the block comes after the
  isomer columns, and a second block may follow for the isomer itself.
- **paths**: an alpha link is Z−2, A−4 (code 100 in `path.txt`). Any number of
  alpha decays can appear in a chain, mixed with β/β1n/β2n. List every
  species you want followed. Only the parent has to be the first row; the
  other rows can be in any order, since paths are built outward from the
  parent. pgenfit built paths in file order, which breaks below an alpha
  decay. A species that is not listed ends the chain, but its decay is still
  counted.
- **model**: `pa<k>` is the alpha fraction of species k. Beta links carry
  `(1−pa)·(p0n|p1n|p2n)` (P1n/P2n are fractions of the beta decays), and alpha
  links carry `pa`. A decay of species k is detected with weight
  `(1−pa)+pa·ab`; for the parent the weight is `(1−pa)·be+pa·ab`.
  - `ab` = ε<sub>α</sub>/ε<sub>β</sub>, the alpha/beta detection-efficiency ratio,
    from the optional effparms columns 10–11 (`ab err`, negative = float,
    default 1).
  - `be` keeps its pgenfit2 meaning (parent-only beta factor). In pgenfit's T12
    alpha mode, `be` was instead ε<sub>β</sub>/ε<sub>α</sub> for all species.
  - Species without an alpha branch get no extra factor, so beta-only
    networks generate byte-identical models to before (checked on
    As92, the isomer example and pgenfit's example).
- **simulator**: a member alpha-decays with probability AlphaBR. The hit is
  detected with `alphaeff`, placed with its own spread (`deltaxylimitalpha`,
  `dx/dyalphamean`, `dx/dyalphasigma`), and stored with the betas (mode 10,
  no neutrons). The correlation `tree` has a new `alpha` branch (1 = alpha).
  pgenfit's `z` tag used the member id, so the parent's own alphas looked
  untagged. Beta-only chains keep the same random sequence, so their toys are
  identical to before.

Examples:
- `examples/alpha_chain_parms.txt`: the 218Pb chain followed through
  218Po/218At α → 214Pb/214Bi → 214Po α → 210Tl (illustrative values).
- `examples/alpha_test_parms.txt`, `alpha_test_effparms.txt` (`ab=2`) and
  `simparms_alpha_test.txt`: a toy chain with a floated 60 % alpha branch, a
  100 % alpha emitter, and a species reached through both an alpha and a beta
  path. One toy gave T<sub>1/2</sub> = 0.2041 ± 0.0070 s (truth 0.2) and
  pa = 0.647 ± 0.066 (truth 0.60).

### Correlation window (long-lived parents)

`simulation.cc` correlates each decay with implants from `ionbetawindowlow` s
before to `ionbetawindowup` s after it. These are simparms keys, default 10 / 20 s,
which was hard-coded in pgenfit/pgenfit2. The window is written into the output file (`TParameter`s
`ionbetawindowlow/up`), and `fit_decay_curve` refuses a `timeRange` larger than
it, since the forward and backward regions would otherwise be silently truncated.
The accidental background grows in proportion to the window. Keep the
run-boundary margin (`boundaryMarginSec`, default 900 s) longer than the window.

Example: `parms/At224parms.txt` (224At → 224Rn → 224Fr → 224Ra α → 220Rn α →
216Po α → 212Pb; values from `pgenfit/parms/decaydatanudat3.csv`) with
`examples/simparms_At224.txt` (±600 s window, 100000 s beam at 0.05 /s,
backgrounds ÷10): one simulation takes 3.5 GB and 7.5 s. Fit with
`timeRange=600`:

```bash
./run_single_toy.sh parms/At224parms.txt examples/simparms_At224.txt <seed> \
    parms/At224parms.txt_effparms 900 4 at224 0 600
```

### Mixture of implanted species

When several ion species are implanted and not separated (e.g. neighbouring
nuclides in the same PID gate), their decay curves add up.

- **simulation**: `simulation_mix <out.root> <seed> <parmsexA> <simparmsA> <parmsexB> <simparmsB> ...`
  runs each species with its own network and beam/detector settings (implant
  rate = its simparms `beamrate`) over the same beam time, merges all hits,
  and correlates them together, so each decay is also accidentally correlated
  with the other species' implants. The correlation settings come from the
  first simparms, and the backgrounds of all simparms add up, so give them
  once and set `betabkgrateg/betabkgrateu/neubkgrate/r2neubkgrate 0` in the
  others. The file stores the true implant counts (`TParameter nimplant_s<i>`)
  and the correlation tree has a branch `ionspecies` (truth).
- **fit**: pass a comma-separated parmsex list (`A.txt,B.txt`) to `build.sh`
  and `fit_decay_curve`. Each species keeps its own network and parameters,
  prefixed `s<i>_` (`s0_l0`, `s1_l0`, …), while `bkg`, `bkga`, `be` and `ab` are
  shared. The parent normalisations become `s<i>_N0raw = Nimp·frac_s<i>·s<i>_l0`:
  `Nimp` is the total number of implants and **`frac_s<i>` the implant fraction
  of species i**, the new fitted parameter (the last species gets 1 − the
  others). The implant ratio `frac_s0/frac_s1` is printed with its error
  (`ratio01` on the `PULLDATA` line; `frac_s<i>true`/`ratio01true` from the
  file). Alpha gating works the same way. A single parmsex file gives exactly
  the single-species fit.
- **one toy**: `./run_mix_toy.sh <parmsexA,parmsexB> <simparmsA,simparmsB> <seed> [effparms] ...`

Example: `examples/alpha_test_parms.txt` (214Pb chain, 1/s, all
backgrounds) + `examples/mix_B_parms.txt` (100Rb 0.05 s → 100Sr → 100Y, 0.5/s,
backgrounds 0 in `examples/simparms_mix_B.txt`); true fraction 2/3. One toy
(0.76 GB, 1 s simulation):

| | total fit | alpha-gated fit |
|---|---|---|
| frac_s0 (truth 0.668) | 0.671 ± 0.021 | 0.667 ± 0.005 |
| implant ratio (truth 2.01) | 2.04 ± 0.19 | 2.005 ± 0.046 |

Only species 0 emits alphas, so the alpha gate separates the two much
better.

## Memory safety

A single `simulation` process with `simparmsex.txt` (36000 s beam, 130/s
beta background) peaks at **11.2 GB RSS**, since it keeps every hit in memory.
Running ~40 of them at once took the server down on 2026-09-25. Now:

- every script limits itself to **50 % of MemTotal**
  (`PGENFIT3_MEM_BUDGET_GB` may lower this, never raise it);
- each process gets a hard cap with `ulimit -v`. `run_pull_study.sh` gives each
  worker `MEM_PER_WORKER_GB` (default 8; use 14 for the As92 settings) and
  reduces the worker count so that workers × cap ≤ budget;
- a watchdog sums the RSS of every child process every 2 s and kills them all
  if the total exceeds the budget.

Measure one toy before a large study:
`/usr/bin/time -v ./simulation ...` → "Maximum resident set size".

## Simulator fixes relative to pgenfit/pgenfit2

1. **Isomer population** (`simulation.cc`, path-flow test): `&&` → `||`.
   Before, the unpopulated partner of an isomer pair also decayed in every
   event, and a chain with a granddaughter below an isomer aborted with
   "something wrong".
2. **Decay time of multi-path species**: the time of the path that actually
   flowed is used. Before, it was the time of the last path examined, which
   is wrong for any species reachable by more than one path (e.g. Br91 in the
   As92 chain, reached through both Se91 and Se92).
3. `dybetamean` was written into `dxbetamean`.

`decaypath.cc`: an isomer row's P1n/P2n/neutron-efficiency upper errors were
read into the ground-state entry. This was cosmetic, since neither the fit
nor the simulation uses those fields.

Test of fixes 1+2 (parent → ground state/isomer with ratio 0.3 → granddaughter,
ratio floated): fitted ground-state fraction 0.69–0.71 on 6 toys (truth 0.70),
decay-curve χ²/ndof ≈ 1. Before the fixes the same fit gave 0.35.
