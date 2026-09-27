package wave

import (
	"math"
	"testing"
)

// gaussianIC builds a centered Gaussian displacement pulse.
func gaussianIC(nx, ny int) []float64 {
	d := make([]float64, nx*ny)
	ci, cj := nx/2, ny/2
	sigma := float64(nx) / 10
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			di, dj := float64(i-ci), float64(j-cj)
			d[j*nx+i] = math.Exp(-(di*di + dj*dj) / (2 * sigma * sigma))
		}
	}
	return d
}

// TestMixedBoundariesStaySelfConsistent runs a field combining fixed and
// absorbing edges and checks: the run stays finite, fixed edges hold
// u == 0 exactly, and corners touching a fixed edge are 0.
func TestMixedBoundariesStaySelfConsistent(t *testing.T) {
	nx, ny := 48, 40
	cfg := Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.008, Steps: 120,
		Boundaries: Boundaries{
			Left:   BoundaryFixed,
			Right:  BoundaryAbsorbing,
			Top:    BoundaryFixed,
			Bottom: BoundaryAbsorbing,
		},
		Init: InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
	}
	if err := cfg.Validate(); err != nil {
		t.Fatalf("config invalid: %v", err)
	}
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatalf("run with mixed boundaries failed: %v", err)
	}
	g := s.Snapshot()
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			if v := g.At(i, j); math.IsNaN(v) || math.IsInf(v, 0) {
				t.Fatalf("non-finite value at (%d,%d)", i, j)
			}
		}
	}
	for j := 0; j < ny; j++ {
		if v := g.At(0, j); v != 0 {
			t.Fatalf("fixed left edge not zero at j=%d: %v", j, v)
		}
		if v := g.At(nx-1, j); math.IsNaN(v) {
			t.Fatalf("absorbing right edge NaN at j=%d", j)
		}
	}
	for i := 0; i < nx; i++ {
		if v := g.At(i, ny-1); v != 0 {
			t.Fatalf("fixed top edge not zero at i=%d: %v", i, v)
		}
	}
	// Corners touching at least one fixed edge must be exactly 0; the
	// bottom-right corner touches two absorbing edges and follows Mur.
	fixedCorners := [][2]int{{0, 0}, {0, ny - 1}, {nx - 1, ny - 1}}
	for _, p := range fixedCorners {
		if v := g.At(p[0], p[1]); v != 0 {
			t.Fatalf("corner %v touches a fixed edge but is %v", p, v)
		}
	}
	if v := g.At(nx-1, 0); math.IsNaN(v) || math.IsInf(v, 0) {
		t.Fatalf("absorbing-absorbing corner not finite: %v", v)
	}
}

// TestAllBoundaryCombinationsRunStably enumerates every combination of
// edge types (2^4 = 16) and requires each run to stay finite and bounded.
func TestAllBoundaryCombinationsRunStably(t *testing.T) {
	nx, ny := 32, 32
	for mask := 0; mask < 16; mask++ {
		b := Boundaries{}
		edges := []*BoundaryType{&b.Left, &b.Right, &b.Top, &b.Bottom}
		for k, e := range edges {
			if mask&(1<<k) != 0 {
				*e = BoundaryAbsorbing
			}
		}
		cfg := Config{
			Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.01, Steps: 200,
			Boundaries: b,
			Init:       InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
		}
		s, err := NewSolver(cfg)
		if err != nil {
			t.Fatalf("mask %04b: %v", mask, err)
		}
		if err := s.Run(); err != nil {
			t.Fatalf("mask %04b: run failed: %v", mask, err)
		}
		d := s.Diagnostics()
		if math.IsNaN(d.Energy) || math.IsInf(d.Energy, 0) {
			t.Fatalf("mask %04b: non-finite energy %v", mask, d.Energy)
		}
		if d.MaxAbs > 10 { // initial amplitude is 1; leapfrog is stable under CFL
			t.Fatalf("mask %04b: unbounded growth, MaxAbs=%v", mask, d.MaxAbs)
		}
	}
}
