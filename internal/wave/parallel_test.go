package wave

import (
	"testing"
)

// perfConfig returns a valid config for reproducibility/parallel tests.
func perfConfig(workers int) Config {
	nx, ny := 64, 64
	return Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.008, Steps: 60,
		Boundaries: Boundaries{
			Left: BoundaryAbsorbing, Right: BoundaryAbsorbing,
			Top: BoundaryFixed, Bottom: BoundaryAbsorbing,
		},
		Init:    InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
		Seed:    42,
		Workers: workers,
	}
}

func runSnapshot(t *testing.T, cfg Config) *Grid {
	t.Helper()
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	return s.Snapshot()
}

// TestReproducibleWithFixedSeed runs the same configuration twice and
// requires bit-identical results.
func TestReproducibleWithFixedSeed(t *testing.T) {
	a := runSnapshot(t, perfConfig(1))
	b := runSnapshot(t, perfConfig(1))
	if !snapshotsEqual(a, b) {
		t.Fatal("identical configurations produced different results")
	}
}

// TestParallelMatchesSerial requires multi-worker execution to produce
// exactly the same result as single-worker execution (the contract only
// demands agreement within a tolerance; the row-partitioned update is
// deterministic, so the results are bit-identical).
func TestParallelMatchesSerial(t *testing.T) {
	serial := runSnapshot(t, perfConfig(1))
	for _, workers := range []int{2, 4, 8, 0 /* 0 = NumCPU */} {
		par := runSnapshot(t, perfConfig(workers))
		if !snapshotsEqual(serial, par) {
			t.Fatalf("workers=%d result differs from serial", workers)
		}
	}
}

// benchConfig builds a larger grid so parallel speedup is observable.
func benchConfig(workers int) Config {
	nx, ny := 512, 512
	return Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.001, Steps: 20,
		Init:    InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
		Workers: workers,
	}
}

func BenchmarkRun1Worker(b *testing.B)  { benchmarkRun(b, 1) }
func BenchmarkRun4Workers(b *testing.B) { benchmarkRun(b, 4) }
func BenchmarkRun8Workers(b *testing.B) { benchmarkRun(b, 8) }

func benchmarkRun(b *testing.B, workers int) {
	cfg := benchConfig(workers)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		s, err := NewSolver(cfg)
		if err != nil {
			b.Fatal(err)
		}
		if err := s.Run(); err != nil {
			b.Fatal(err)
		}
	}
}
