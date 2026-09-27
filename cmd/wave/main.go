// Command wave runs an in-memory 2D scalar wave equation simulation and
// prints the discrete solution at probe points together with diagnostics.
//
// Inputs are constructed entirely in memory (no network, no files apart
// from the local checkpoint directory). Passing -resume continues the run
// from the latest committed checkpoint instead of starting at t=0.
package main

import (
	"flag"
	"fmt"
	"strings"

	solver "github.com/guoyyit/wave-equation-solver-go"
)

func main() {
	if err := run(); err != nil {
		fmt.Println("error:", err)
	}
}

func run() error {
	var (
		nx        = flag.Int("nx", 80, "grid intervals in x")
		ny        = flag.Int("ny", 80, "grid intervals in y")
		lx        = flag.Float64("lx", 10, "domain length in x")
		ly        = flag.Float64("ly", 10, "domain length in y")
		c         = flag.Float64("c", 1, "wave speed")
		tmax      = flag.Float64("tmax", 8, "final physical time")
		cfl       = flag.Float64("cfl", 0.9, "CFL safety factor in (0,1]")
		workers   = flag.Int("workers", 1, "number of worker goroutines")
		seed      = flag.Int64("seed", 1, "random seed for in-memory input")
		random    = flag.Bool("random-ic", false, "use a seeded random initial displacement")
		bounds    = flag.String("boundaries", "absorb", "comma/space list of west,east,south,north kinds (fixed|absorb), or 'allfixed'/'absorb'")
		ckptDir   = flag.String("checkpoint-dir", "./.wave-checkpoint", "checkpoint directory (empty disables)")
		ckptEvery = flag.Int("checkpoint-every", 10, "commit a checkpoint every N steps")
		resume    = flag.Bool("resume", false, "continue from the latest committed checkpoint")
	)
	flag.Parse()

	west, east, south, north, err := parseBoundaries(*bounds)
	if err != nil {
		return err
	}

	dx := *lx / float64(*nx)
	dy := *ly / float64(*ny)
	var disp solver.FieldFunc
	if *random {
		disp = solver.RandomField(*seed, 0.5)
	} else {
		// Offset narrow pulse so reflections / absorption are clearly
		// visible. Narrow enough to be negligible (< 1e-8) on fixed walls
		// for the default grid, keeping any boundary combination valid.
		disp = solver.GaussianPulseXY(dx, dy, *lx*0.25, *ly*0.5, 1.5*dx, 1.0)
	}
	vel := solver.ZeroField()

	cfg := solver.Config{
		Lx: *lx, Ly: *ly, Nx: *nx, Ny: *ny,
		C: *c, TMax: *tmax, CFLFactor: *cfl, Workers: *workers,
		RandomSeed: *seed,
		Boundary:   solver.Boundaries{West: west, East: east, South: south, North: north},
		Initial:    solver.InitialConditions{Displacement: disp, Velocity: vel},
	}
	if *ckptDir != "" {
		cfg.CheckpointDir = *ckptDir
		cfg.CheckpointInterval = *ckptEvery
	}

	var s *solver.Solver
	if *resume {
		s, err = solver.Resume(*ckptDir, cfg)
		if err != nil {
			return err
		}
		fmt.Printf("resumed from step %d t=%.6f, continuing to %d\n", s.StepCount(), s.Time(), s.TargetSteps())
	} else {
		s, err = solver.NewSolver(cfg)
		if err != nil {
			return err
		}
	}

	if err := s.Run(); err != nil {
		return err
	}

	d := s.Diagnostics()
	fmt.Printf("done: steps=%d t=%.6f max|u|=%.6e energy=%.6e courant=%.4f workers=%d\n",
		d.Step, d.Time, d.MaxAbsU, d.Energy, d.Courant, d.Workers)

	// Probe a few fixed grid points (center and quarter points).
	nxPts := *nx + 1
	nyPts := *ny + 1
	probes := [][2]int{
		{nxPts / 2, nyPts / 2},
		{nxPts / 4, nyPts / 2},
		{3 * nxPts / 4, nyPts / 2},
	}
	for _, p := range probes {
		fmt.Printf("u(%d,%d) = %.8e\n", p[0], p[1], s.At(p[0], p[1]))
	}
	return nil
}

func parseBoundaries(spec string) (w, e, so, n solver.BoundaryKind, err error) {
	switch strings.TrimSpace(strings.ToLower(spec)) {
	case "absorb", "absorbing", "allabsorbing":
		return solver.BoundaryAbsorbing, solver.BoundaryAbsorbing,
			solver.BoundaryAbsorbing, solver.BoundaryAbsorbing, nil
	case "allfixed", "fixed":
		return solver.BoundaryFixed, solver.BoundaryFixed,
			solver.BoundaryFixed, solver.BoundaryFixed, nil
	}
	parts := strings.FieldsFunc(spec, func(r rune) bool { return r == ',' || r == ' ' })
	if len(parts) != 4 {
		err = fmt.Errorf("invalid -boundaries %q: use 'absorb', 'allfixed', or 4 comma-separated kinds", spec)
		return
	}
	parse := func(p string) (solver.BoundaryKind, error) {
		switch strings.ToLower(strings.TrimSpace(p)) {
		case "fixed":
			return solver.BoundaryFixed, nil
		case "absorb", "absorbing":
			return solver.BoundaryAbsorbing, nil
		default:
			return 0, fmt.Errorf("unknown boundary kind %q", p)
		}
	}
	if w, err = parse(parts[0]); err != nil {
		return
	}
	if e, err = parse(parts[1]); err != nil {
		return
	}
	if so, err = parse(parts[2]); err != nil {
		return
	}
	n, err = parse(parts[3])
	return
}
