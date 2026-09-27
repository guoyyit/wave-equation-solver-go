package solver

import (
	"os"
	"path/filepath"
	"testing"
	"time"
)

// TestCheckpointAtomicNoTmp ensures a committed checkpoint leaves no
// partial temp files behind and the checkpoint file is readable.
func TestCheckpointAtomicNoTmp(t *testing.T) {
	dir := t.TempDir()
	cfg := baseConfig()
	cfg.CheckpointDir = dir
	cfg.CheckpointInterval = 3
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	if tmps := ListTempFiles(dir); len(tmps) != 0 {
		t.Fatalf("partial checkpoint temp files remain: %v", tmps)
	}
	if _, err := os.Stat(filepath.Join(dir, checkpointName)); err != nil {
		t.Fatalf("committed checkpoint missing: %v", err)
	}
}

// TestExplicitDtMustDivideTMax guards the integer-step requirement when
// the caller supplies Dt directly.
func TestExplicitDtMustDivideTMax(t *testing.T) {
	cfg := baseConfig()
	cfg.Dt = 0.05
	// TMax=6 divides evenly by 0.05 -> valid.
	if _, err := NewSolver(cfg); err != nil {
		t.Fatalf("expected valid explicit dt, got %v", err)
	}
	cfg.Dt = 0.07 // 6/0.07 not integer -> rejected
	if _, err := NewSolver(cfg); err == nil {
		t.Fatal("expected rejection for non-dividing explicit dt")
	}
}

// TestParallelSpeedup checks a measurable wall-clock improvement on a
// large enough grid. It is a coarse, machine-sensitive smoke benchmark:
// multi-threaded must not be dramatically slower than single-threaded.
func TestParallelSpeedup(t *testing.T) {
	if testing.Short() {
		t.Skip()
	}
	build := func(w int) Config {
		cfg := baseConfig()
		cfg.Nx, cfg.Ny = 600, 600
		cfg.Lx, cfg.Ly = 60, 60
		cfg.TMax = 1
		cfg.CFLFactor = 0.5
		cfg.Workers = w
		return cfg
	}
	// Warmup.
	runWith(t, build(4))

	s1 := runWithTiming(t, build(1))
	s8 := runWithTiming(t, build(8))
	t.Logf("single=%v parallel(8)=%v speedup=%.2fx",
		s1, s8, float64(s1)/float64(s8))
	if s8 > s1 {
		t.Fatalf("parallel run %v slower than serial %v", s8, s1)
	}
}

func runWithTiming(t *testing.T, cfg Config) time.Duration {
	t.Helper()
	start := time.Now()
	runWith(t, cfg)
	return time.Since(start)
}
