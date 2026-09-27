package solver

import "testing"

func fieldMaxDiff(a, b []float64) (int, float64) {
	worst := 0
	var d float64
	for i := range a {
		x := a[i] - b[i]
		if x < 0 {
			x = -x
		}
		if x > d {
			d = x
			worst = i
		}
	}
	return worst, d
}

// TestReproducibleSingleVsParallel runs the identical seeded problem with
// 1, 2, 4 and 8 workers and requires bit-identical results, both for a
// smooth pulse and a seeded random initial field.
func TestReproducibleSingleVsParallel(t *testing.T) {
	workers := []int{1, 2, 4, 8}
	var ref []float64
	for ri, random := range []bool{false, true} {
		ref = nil
		for _, w := range workers {
			cfg := baseConfig()
			cfg.Workers = w
			cfg.Nx, cfg.Ny = 70, 60
			cfg.TMax = 5
			cfg.RandomSeed = 42
			if random {
				// A random field touches boundaries with nonzero values,
				// which contradicts fixed walls; use absorbing edges.
				cfg.Initial.Displacement = RandomField(42, 0.5)
				cfg.Boundary = Boundaries{
					BoundaryAbsorbing, BoundaryAbsorbing,
					BoundaryAbsorbing, BoundaryAbsorbing,
				}
			}
			s := runWith(t, cfg)
			got := s.Solution()
			if ref == nil {
				ref = got
				continue
			}
			if idx, d := fieldMaxDiff(ref, got); d != 0 {
				t.Fatalf("case %d random=%v workers=%d: max diff %.3e at %d (want bit-identical)",
					ri, random, w, d, idx)
			}
		}
	}
}

// TestRandomFieldDeterministic checks the seeded field is order-independent.
func TestRandomFieldDeterministic(t *testing.T) {
	f := RandomField(7, 1)
	g := RandomField(7, 1)
	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			if f(i, j) != g(i, j) {
				t.Fatalf("random field not reproducible at %d,%d", i, j)
			}
		}
	}
	if f(0, 0) == f(1, 0) && f(1, 0) == f(0, 1) {
		t.Fatal("random field appears constant")
	}
}
