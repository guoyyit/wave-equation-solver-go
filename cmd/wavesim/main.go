// Command wavesim runs a 2D scalar wave equation simulation locally.
// All input is constructed in memory from flags; no network or
// credentials are involved. The same flags always reproduce the same
// result, and a run interrupted before completion resumes from its
// newest committed checkpoint.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"strings"

	"github.com/guoyyit/wave-equation-solver-go/internal/wave"
)

type samplePoints [][2]int

func (s *samplePoints) String() string { return fmt.Sprint(*s) }
func (s *samplePoints) Set(v string) error {
	for _, p := range strings.Split(v, ";") {
		var i, j int
		if _, err := fmt.Sscanf(p, "%d,%d", &i, &j); err != nil {
			return fmt.Errorf("bad sample point %q, want i,j", p)
		}
		*s = append(*s, [2]int{i, j})
	}
	return nil
}

func parseBoundary(s string) (wave.BoundaryType, error) {
	switch s {
	case "fixed":
		return wave.BoundaryFixed, nil
	case "absorbing":
		return wave.BoundaryAbsorbing, nil
	}
	return 0, fmt.Errorf("unknown boundary type %q (fixed|absorbing)", s)
}

// gaussian builds a smooth deterministic pulse, used as initial condition.
func gaussian(nx, ny int, lx, ly, amp, sigma float64) []float64 {
	out := make([]float64, nx*ny)
	cx, cy := lx/2, ly/2
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			x := float64(i) * lx / float64(nx-1)
			y := float64(j) * ly / float64(ny-1)
			r2 := (x-cx)*(x-cx) + (y-cy)*(y-cy)
			out[j*nx+i] = amp * math.Exp(-r2/(2*sigma*sigma))
		}
	}
	return out
}

func main() {
	var (
		nx, ny     = flag.Int("nx", 64, "grid points in x"), flag.Int("ny", 64, "grid points in y")
		lx, ly     = flag.Float64("lx", 1, "extent in x (m)"), flag.Float64("ly", 1, "extent in y (m)")
		c          = flag.Float64("c", 1, "wave speed (m/s)")
		dt         = flag.Float64("dt", 0.005, "time step (s); must satisfy CFL")
		steps      = flag.Int("steps", 100, "number of time steps")
		bl         = flag.String("boundary-left", "absorbing", "left boundary: fixed|absorbing")
		br         = flag.String("boundary-right", "absorbing", "right boundary: fixed|absorbing")
		bt         = flag.String("boundary-top", "absorbing", "top boundary: fixed|absorbing")
		bb         = flag.String("boundary-bottom", "absorbing", "bottom boundary: fixed|absorbing")
		ic         = flag.String("ic", "gaussian", "initial displacement: gaussian|zero")
		vc         = flag.String("vc", "zero", "initial velocity: gaussian|zero")
		workers    = flag.Int("workers", 0, "parallel workers (0 = NumCPU)")
		ckpt       = flag.String("checkpoint", "", "checkpoint file path (enables resume)")
		ckptEvery  = flag.Int("checkpoint-every", 10, "commit a checkpoint every N steps")
		seed       = flag.Int64("seed", 0, "reserved random seed")
	)
	var samples samplePoints
	flag.Var(&samples, "sample", "sample points i,j separated by ';' (default: grid center)")
	flag.Parse()

	left, err := parseBoundary(*bl)
	if err != nil {
		fatal(err)
	}
	right, err := parseBoundary(*br)
	if err != nil {
		fatal(err)
	}
	top, err := parseBoundary(*bt)
	if err != nil {
		fatal(err)
	}
	bottom, err := parseBoundary(*bb)
	if err != nil {
		fatal(err)
	}

	cfg := wave.Config{
		Nx: *nx, Ny: *ny, Lx: *lx, Ly: *ly, C: *c, Dt: *dt, Steps: *steps,
		Boundaries:      wave.Boundaries{Left: left, Right: right, Top: top, Bottom: bottom},
		Seed:            *seed,
		Workers:         *workers,
		CheckpointPath:  *ckpt,
		CheckpointEvery: *ckptEvery,
	}
	sigma := 0.05 * math.Min(*lx, *ly)
	switch *ic {
	case "gaussian":
		cfg.Init.Displacement = gaussian(*nx, *ny, *lx, *ly, 1, sigma)
	case "zero":
		cfg.Init.DisplacementZero = true
	default:
		fatal(fmt.Errorf("unknown -ic %q (gaussian|zero)", *ic))
	}
	switch *vc {
	case "gaussian":
		cfg.Init.Velocity = gaussian(*nx, *ny, *lx, *ly, 1, sigma)
	case "zero":
		cfg.Init.VelocityZero = true
	default:
		fatal(fmt.Errorf("unknown -vc %q (gaussian|zero)", *vc))
	}

	solver, err := wave.NewSolver(cfg)
	if err != nil {
		fatal(err)
	}
	resumed := solver.StepCount() > 0
	if err := solver.Run(); err != nil {
		fatal(err)
	}

	if len(samples) == 0 {
		samples = samplePoints{{*nx / 2, *ny / 2}}
	}
	snap := solver.Snapshot()
	type sample struct {
		I, J  int
		Value float64
	}
	out := struct {
		Resumed     bool
		Diagnostics wave.Diagnostics
		Samples     []sample
	}{Resumed: resumed, Diagnostics: solver.Diagnostics()}
	for _, p := range samples {
		if p[0] < 0 || p[0] >= *nx || p[1] < 0 || p[1] >= *ny {
			fatal(fmt.Errorf("sample point %v outside grid %dx%d", p, *nx, *ny))
		}
		out.Samples = append(out.Samples, sample{p[0], p[1], snap.At(p[0], p[1])})
	}
	enc := json.NewEncoder(os.Stdout)
	enc.SetIndent("", "  ")
	if err := enc.Encode(out); err != nil {
		fatal(err)
	}
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "wavesim:", err)
	os.Exit(1)
}
