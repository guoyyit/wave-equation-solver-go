package solver

import (
	"path/filepath"
	"testing"
)

// TestResumeMatchesUninterrupted verifies long-run + checkpoint resume:
// run A completes uninterrupted with checkpointing; run B runs only to a
// midpoint, then resumes from the latest committed checkpoint to the end.
// Final fields and steps must be identical, with no step re-advanced.
func TestResumeMatchesUninterrupted(t *testing.T) {
	dirA := t.TempDir()
	dirB := t.TempDir()

	cfgA := baseConfig()
	cfgA.TMax = 6
	cfgA.Nx, cfgA.Ny = 50, 50
	cfgA.Initial.Displacement = GaussianPulse(25, 25, 3, 1)
	cfgA.CheckpointDir = dirA
	cfgA.CheckpointInterval = 7
	sA, err := NewSolver(cfgA)
	if err != nil {
		t.Fatal(err)
	}
	if err := sA.Run(); err != nil {
		t.Fatal(err)
	}
	want := sA.Solution()

	// Fresh run B, stop at the same checkpoint step as A commits.
	cfgB := cfgA
	cfgB.CheckpointDir = dirB
	sB, err := NewSolver(cfgB)
	if err != nil {
		t.Fatal(err)
	}
	mid := 14
	for sB.StepCount() < mid {
		if err := sB.Step(); err != nil {
			t.Fatal(err)
		}
	}
	ckFile := filepath.Join(dirB, checkpointName)
	data, err := loadLatestCheckpoint(dirB)
	if err != nil {
		t.Fatal(err)
	}
	if data.Step != mid {
		t.Fatalf("checkpoint step=%d want %d", data.Step, mid)
	}

	// Simulate a process restart: rebuild via Resume from the committed
	// checkpoint and finish the run. No part-state file should be present.
	sB2, err := Resume(dirB, cfgB)
	if err != nil {
		t.Fatal(err)
	}
	if sB2.StepCount() != mid {
		t.Fatalf("resume must not re-advance: step=%d want %d", sB2.StepCount(), mid)
	}
	if err := sB2.Run(); err != nil {
		t.Fatal(err)
	}
	if sB2.StepCount() != sA.TargetSteps() {
		t.Fatalf("resumed steps=%d want %d", sB2.StepCount(), sA.TargetSteps())
	}

	got := sB2.Solution()
	if len(got) != len(want) {
		t.Fatalf("field length mismatch %d vs %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("resume mismatch at %d: %.17g vs %.17g", i, got[i], want[i])
		}
	}
	_ = ckFile
}

// TestResumeRejectsMismatchedConfig ensures a checkpoint is not silently
// continued under different physics parameters.
func TestResumeRejectsMismatchedConfig(t *testing.T) {
	dir := t.TempDir()
	cfg := baseConfig()
	cfg.CheckpointDir = dir
	cfg.CheckpointInterval = 5
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	for s.StepCount() < 10 {
		if err := s.Step(); err != nil {
			t.Fatal(err)
		}
	}

	bad := cfg
	bad.C = 2.5
	if _, err := Resume(dir, bad); err == nil {
		t.Fatal("expected error resuming with mismatched wave speed")
	}
}
