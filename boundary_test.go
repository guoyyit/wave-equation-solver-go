package solver

import (
	"math"
	"testing"
)

// runWith is a small helper that builds and completes a simulation.
func runWith(t *testing.T, cfg Config) *Solver {
	t.Helper()
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	return s
}

func TestBoundaryCombinations(t *testing.T) {
	kinds := []BoundaryKind{BoundaryFixed, BoundaryAbsorbing}
	// Exercise a representative set of mixed combinations, including all
	// four edges independently toggled.
	combos := []Boundaries{
		{BoundaryFixed, BoundaryAbsorbing, BoundaryFixed, BoundaryAbsorbing},
		{BoundaryAbsorbing, BoundaryFixed, BoundaryAbsorbing, BoundaryFixed},
		{BoundaryFixed, BoundaryFixed, BoundaryAbsorbing, BoundaryAbsorbing},
		{BoundaryAbsorbing, BoundaryAbsorbing, BoundaryFixed, BoundaryFixed},
	}
	_ = kinds
	for bi, b := range combos {
		cfg := baseConfig()
		cfg.Boundary = b
		s := runWith(t, cfg)

		// Fixed edges must be exactly zero.
		nx, ny := s.nx, s.ny
		u := s.Solution()
		checkZero := func(name string, cond func(i, j int) bool) {
			for i := 0; i < nx; i++ {
				for j := 0; j < ny; j++ {
					if cond(i, j) && u[i*ny+j] != 0 {
						t.Fatalf("combo %d: fixed %s edge nonzero at (%d,%d)=%v", bi, name, i, j, u[i*ny+j])
					}
				}
			}
		}
		if b.West == BoundaryFixed {
			checkZero("west", func(i, j int) bool { return i == 0 && j > 0 && j < ny-1 })
		}
		if b.East == BoundaryFixed {
			checkZero("east", func(i, j int) bool { return i == nx-1 && j > 0 && j < ny-1 })
		}
		if b.South == BoundaryFixed {
			checkZero("south", func(i, j int) bool { return j == 0 && i > 0 && i < nx-1 })
		}
		if b.North == BoundaryFixed {
			checkZero("north", func(i, j int) bool { return j == ny-1 && i > 0 && i < nx-1 })
		}
		// Corners always zero regardless of combination.
		for _, idx := range []int{0, ny - 1, (nx - 1) * ny, (nx-1)*ny + ny - 1} {
			if u[idx] != 0 {
				t.Fatalf("combo %d: corner %d not zero: %v", bi, idx, u[idx])
			}
		}
		// Everything finite.
		for _, v := range u {
			if math.IsNaN(v) || math.IsInf(v, 0) {
				t.Fatalf("combo %d: non-finite value", bi)
			}
		}
	}
}

// TestAbsorbingBoundaryReducesReflection compares a pulse near one edge
// between an all-fixed reference and an absorbing run. The reflected
// energy measured after the wave has had time to travel to the boundary
// and back must be substantially lower in the absorbing case.
func TestAbsorbingBoundaryReducesReflection(t *testing.T) {
	build := func(b BoundaryKind) Config {
		cfg := baseConfig()
		cfg.Nx, cfg.Ny = 120, 80
		cfg.Lx, cfg.Ly = 30, 10
		cfg.TMax = 14
		cfg.CFLFactor = 0.8
		dx := cfg.Lx / float64(cfg.Nx)
		// Pulse close to the west edge so it hits west/east quickly.
		// Narrow enough to be negligible on fixed walls (IC consistency),
		// wide enough to resolve on the grid.
		cfg.Initial.Displacement = GaussianPulseXY(dx, cfg.Ly/float64(cfg.Ny), 3.2, 5, 0.45, 1)
		cfg.Boundary = Boundaries{b, b, b, b}
		return cfg
	}

	fixed := runWith(t, build(BoundaryFixed)).Diagnostics()
	abs := runWith(t, build(BoundaryAbsorbing)).Diagnostics()

	t.Logf("fixed energy=%v max|u|=%v | absorbing energy=%v max|u|=%v",
		fixed.Energy, fixed.MaxAbsU, abs.Energy, abs.MaxAbsU)

	// By t=14 a speed-1 pulse has crossed the 30-long domain (to the far
	// wall and partway back). Reflections from fixed walls must keep far
	// more residual energy than the absorbing run.
	if abs.Energy >= fixed.Energy {
		t.Fatalf("absorbing energy %v not lower than fixed %v", abs.Energy, fixed.Energy)
	}
	ratio := abs.Energy / fixed.Energy
	if ratio > 0.35 {
		t.Fatalf("absorbing energy ratio %.3f too high (reflection not sufficiently suppressed)", ratio)
	}
}
