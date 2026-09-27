package wave

import (
	"testing"
)

// runEnergy runs a centered Gaussian pulse with the given boundaries and
// returns the residual total energy after steps steps.
func runEnergy(t *testing.T, b Boundaries, steps int) float64 {
	t.Helper()
	nx, ny := 64, 64
	cfg := Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.008, Steps: steps,
		Boundaries: b,
		Init:       InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
	}
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	return s.Diagnostics().Energy
}

// TestAbsorbingBoundaryReducesReflection compares the residual field
// energy long after the pulse has reached the boundary: with absorbing
// edges the energy must drain to a small fraction of the all-fixed
// reference, demonstrating suppressed reflection.
func TestAbsorbingBoundaryReducesReflection(t *testing.T) {
	fixed := Boundaries{}
	absorbing := Boundaries{
		Left: BoundaryAbsorbing, Right: BoundaryAbsorbing,
		Top: BoundaryAbsorbing, Bottom: BoundaryAbsorbing,
	}
	const steps = 200 // ~1.6 s: pulse crosses the domain several times

	eFixed := runEnergy(t, fixed, steps)
	eAbs := runEnergy(t, absorbing, steps)
	t.Logf("residual energy: fixed=%.6f absorbing=%.6f ratio=%.4f",
		eFixed, eAbs, eAbs/eFixed)

	if eFixed <= 0 {
		t.Fatal("fixed-boundary reference lost all energy; test setup broken")
	}
	if ratio := eAbs / eFixed; ratio > 0.2 {
		t.Fatalf("absorbing boundary did not significantly reduce reflected energy: ratio=%.4f > 0.2", ratio)
	}
}
