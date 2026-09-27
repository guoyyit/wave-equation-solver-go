# 2D Scalar Wave Equation Solver (Go)

A self-contained, dependency-free Go service that numerically solves the
two-dimensional scalar wave equation on a rectangular domain, with fixed
and absorbing boundaries, stability enforcement, crash-safe checkpoint /
resume, and deterministic single-/multi-threaded execution.

    ∂²u/∂t² = C² ( ∂²u/∂x² + ∂²u/∂y² )      on [0,Lx]×[0,Ly], t∈[0,TMax]

Inputs are constructed **in memory** (no network, no credentials). A local
command line entry point reproduces every run. Outputs are the discrete
solution at grid/probe points plus diagnostics (max amplitude, discrete
energy, Courant number, step/time, worker count).

## Build and test

    go build ./...
    go test ./... -count=1
    go run ./cmd/wave -h

## Quick start

    # all-absorbing run, checkpoints every 10 steps
    go run ./cmd/wave -nx 80 -ny 80 -tmax 8 -boundaries absorb \
        -checkpoint-dir ./.wave-checkpoint -checkpoint-every 10

    # continue the same run after a crash / restart (no step is repeated)
    go run ./cmd/wave -nx 80 -ny 80 -tmax 8 -boundaries absorb \
        -checkpoint-dir ./.wave-checkpoint -checkpoint-every 10 -resume

    # mixed boundaries: west,east,south,north ∈ {fixed, absorb}
    go run ./cmd/wave -boundaries 'fixed,absorb,fixed,absorb' -workers 8

`absorb` / `allfixed` select all four edges at once. Results for
`-workers 1` and `-workers N` are bit-identical.

## Discretization

Uniform grid `x_i = i·dx`, `y_j = j·dy`, `dx = Lx/Nx`, `dy = Ly/Ny`
(the grid has `(Nx+1)·(Ny+1)` points). Flat arrays are indexed
`i*ny + j`.

Time integration is the explicit **leapfrog** (second-order centered)
scheme:

    u^(n+1) = 2u^n − u^(n−1)
            + rx² (u^n[i−1,j] − 2u^n[i,j] + u^n[i+1,j])
            + ry² (u^n[i,  j−1] − 2u^n[i,j] + u^n[i,  j+1]),

    rx = C·dt/dx,   ry = C·dt/dy.

The first step uses a second-order Taylor start so the leapfrog only ever
runs with two valid time levels:

    u^1 = u^0 + dt·v0 + ½·dt²·C²·∇²u^0.

### Boundaries

* **Fixed** (`BoundaryFixed`): homogeneous Dirichlet, `u = 0`.
* **Absorbing** (`BoundaryAbsorbing`): first-order **Mur** condition. For
  an edge normal to x, with edge value `e` and its interior neighbor `v`,

      e^(n+1) = v^n + (r−1)/(r+1) · (v^(n+1) − e^n),   r = C·dt/dx

  (analogously with `r = C·dt/dy` for edges normal to y).

The four kinds can be mixed freely. Corner cells are shared by two edges
and are **always pinned to zero**, a single rule that makes every
fixed/absorbing combination self-consistent (no two edges ever write
conflicting corner values).

## Stability control

Leapfrog is stable only under the 2D CFL condition

    C·dt·sqrt(1/dx² + 1/dy²) < 1   (Courant number < 1).

`Config.Validate()` computes the Courant number and **rejects** any
configuration with `courant ≥ 1 − ε` (`ε = 1e-9`), returning an error
wrapping `solver.ErrInvalidConfig` that names the offending `dt`, grid
and speed. There is therefore no parameter combination that can drive the
explicit update into unbounded growth; `Solver.Step` additionally guards
against non-finite/runaway values (`solver.ErrUnstable`) as a defensive
backstop.

Time-step resolution:

* If `Dt > 0` it is used as given and `TMax/Dt` must be an integer number
  of steps (otherwise the run is rejected, so restarts and probe times are
  unambiguous).
* If `Dt == 0`, it is derived from `CFLFactor ∈ (0,1]`: the number of
  steps is `n = ceil(TMax / (CFLFactor·dt_max))` and `dt = TMax/n`, so the
  run lands exactly on `TMax` while remaining strictly inside the CFL
  limit.

`CFLStep()` exposes `dt_max = 1 / (C·sqrt(1/dx² + 1/dy²))`.

## Initial conditions

Both fields must be supplied explicitly:

    Initial: InitialConditions{
        Displacement: FieldFunc,   // u0
        Velocity:     FieldFunc,   // v0; use ZeroField() for zero
    }

A nil displacement or velocity is rejected (`ErrInvalidConfig`), never
silently defaulted. A prescribed nonzero displacement or velocity on a
**fixed** edge is also rejected as contradictory. Helpers:
`ZeroField()`, `GaussianPulse(...)`, `GaussianPulseXY(...)`, and
`RandomField(seed, amp)` (deterministic per `(seed,i,j)`).

## Long runs and recovery

* With `CheckpointDir` set, a checkpoint is atomically committed every
  `CheckpointInterval` steps. Writes go to a temp file, are fsynced, then
  renamed over `checkpoint-latest.json`, so a crash always leaves either
  the previous complete checkpoint or the new one — never a torn file.
* `solver.Resume(dir, cfg)` loads the latest checkpoint, verifies a hash
  of the governing parameters (mismatched physics is rejected), and
  rebuilds a solver positioned exactly at the recorded step. The recorded
  step is not re-advanced; a subsequent `Run()` continues from `step+1`.
* A resumed run produces the same final field as one uninterrupted run
  (verified byte-for-byte in tests).

## Reproducibility and parallelism

* Interior rows are split into contiguous stripes, one per goroutine
  (`Workers`). Each element is written by exactly one goroutine and every
  read comes from `u^n/u^(n−1)`; the floating-point arithmetic is the
  same as the serial path, so 1/2/4/… workers give **bit-identical**
  results. `Workers ≤ 1` selects the serial path.
* `RandomField` derives each value solely from `(seed, i, j)` (no shared
  RNG stream), so it is independent of evaluation order and worker count.
* Wall time improves measurably with more workers on large grids
  (`TestParallelSpeedup` asserts the parallel run is not slower and logs
  the speedup).

## Configuration reference (`solver.Config`)

| Field | Meaning |
|---|---|
| `Lx`, `Ly` | Domain lengths (must be > 0) |
| `Nx`, `Ny` | Grid intervals; `(Nx+1)×(Ny+1)` points, `≥ 2` |
| `C` | Wave speed (> 0) |
| `Dt` | Time step; `0` derives it from `CFLFactor` |
| `TMax` | Final time (> 0); integer multiple of `Dt` |
| `CFLFactor` | Safety factor in `(0,1]` used when `Dt == 0` |
| `Boundary` | `West/East/South/North` ∈ `{Fixed, Absorbing}` |
| `Initial` | Explicit `Displacement` and `Velocity` fields |
| `Workers` | Goroutines for interior updates (`≥ 1`) |
| `RandomSeed` | Seed for in-memory random input |
| `CheckpointDir` | Checkpoint directory; empty disables checkpointing |
| `CheckpointInterval` | Commit every N steps (`≥ 1` when dir set) |

## Test coverage

| Area | Test |
|---|---|
| Boundary combinations / fixed edges / corners | `TestBoundaryCombinations` |
| Absorption vs. fixed reference | `TestAbsorbingBoundaryReducesReflection` |
| CFL / unstable parameter rejection | `TestCFLViolationRejected`, `TestInvalidGeometryRejected` |
| Missing / contradictory initial conditions | `TestMissingInitialConditionsRejected`, `TestContradictoryInitialConditionsRejected` |
| Crash recovery, no repeated step | `TestResumeMatchesUninterrupted`, `TestResumeRejectsMismatchedConfig` |
| Atomic checkpoint | `TestCheckpointAtomicNoTmp` |
| Reproducibility / serial vs. parallel | `TestReproducibleSingleVsParallel`, `TestRandomFieldDeterministic` |
| Parallel scaling | `TestParallelSpeedup` |
| Boundedness | `TestSolutionStaysBounded` |

The absorption test compares residual discrete energy after a pulse has
had time to reach the walls; the absorbing run's energy ratio versus the
all-fixed reference is asserted below 0.35 (observed ≈ 0.12).
