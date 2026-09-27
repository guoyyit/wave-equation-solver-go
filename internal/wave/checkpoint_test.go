package wave

import (
	"os"
	"path/filepath"
	"testing"
)

// ckptConfig returns a valid checkpointed config rooted at dir.
func ckptConfig(dir string, steps, every int) Config {
	nx, ny := 32, 32
	return Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.01, Steps: steps,
		Init:            InitialCondition{Displacement: gaussianIC(nx, ny), VelocityZero: true},
		CheckpointPath:  filepath.Join(dir, "state.ckpt"),
		CheckpointEvery: every,
	}
}

func snapshotsEqual(a, b *Grid) bool {
	if a.Nx != b.Nx || a.Ny != b.Ny || len(a.Data) != len(b.Data) {
		return false
	}
	for i := range a.Data {
		if a.Data[i] != b.Data[i] {
			return false
		}
	}
	return true
}

// TestCrashRecoveryReproducesUninterruptedRun simulates a crash after a
// committed checkpoint, resumes in a fresh solver, and requires the
// final state to be bit-identical to an uninterrupted run.
func TestCrashRecoveryReproducesUninterruptedRun(t *testing.T) {
	dir := t.TempDir()

	// Reference: no checkpointing at all.
	refCfg := ckptConfig(dir, 120, 10)
	refCfg.CheckpointPath = ""
	ref, err := NewSolver(refCfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := ref.Run(); err != nil {
		t.Fatal(err)
	}

	// Crashed run: advance 37 steps by hand, commit, then "crash".
	cfg := ckptConfig(dir, 120, 10)
	crashed, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 37; i++ {
		crashed.Step()
	}
	if err := crashed.saveCheckpoint(); err != nil {
		t.Fatal(err)
	}
	// Leave a partial temp file behind, as an interrupted write would.
	if err := os.WriteFile(cfg.CheckpointPath+".tmp", []byte("garbage"), 0o644); err != nil {
		t.Fatal(err)
	}

	// Resume in a fresh solver from the committed checkpoint.
	resumed, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if got := resumed.StepCount(); got != 37 {
		t.Fatalf("resumed at step %d, want 37", got)
	}
	if err := resumed.Run(); err != nil {
		t.Fatal(err)
	}
	if got := resumed.StepCount(); got != 120 {
		t.Fatalf("completed %d steps, want 120", got)
	}
	if !snapshotsEqual(ref.Snapshot(), resumed.Snapshot()) {
		t.Fatal("resumed run diverges from uninterrupted reference")
	}
}

// TestResumeDoesNotRepeatCompletedWork checks that Run on a finished
// solver is a no-op: no repeated advancement, no state change.
func TestResumeDoesNotRepeatCompletedWork(t *testing.T) {
	dir := t.TempDir()
	cfg := ckptConfig(dir, 40, 10)

	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	if got := s.StepCount(); got != 40 {
		t.Fatalf("step count %d, want 40", got)
	}
	done := s.Snapshot()

	// A brand-new solver on the same path resumes the final checkpoint
	// and must not advance further.
	again, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if got := again.StepCount(); got != 40 {
		t.Fatalf("reopened solver at step %d, want 40", got)
	}
	if err := again.Run(); err != nil {
		t.Fatal(err)
	}
	if got := again.StepCount(); got != 40 {
		t.Fatalf("Run advanced a finished solver to %d", got)
	}
	if !snapshotsEqual(done, again.Snapshot()) {
		t.Fatal("state changed after re-running a finished solver")
	}
}

// TestCheckpointCommitIsAtomic verifies that after Run every commit left
// a complete checkpoint file and no temporary residue.
func TestCheckpointCommitIsAtomic(t *testing.T) {
	dir := t.TempDir()
	cfg := ckptConfig(dir, 25, 5)
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(cfg.CheckpointPath); err != nil {
		t.Fatalf("committed checkpoint missing: %v", err)
	}
	if _, err := os.Stat(cfg.CheckpointPath + ".tmp"); !os.IsNotExist(err) {
		t.Fatalf("temporary checkpoint file left behind: %v", err)
	}
	// The committed file must be loadable on its own.
	fresh, err := NewSolver(cfg)
	if err != nil {
		t.Fatalf("committed checkpoint not loadable: %v", err)
	}
	if got := fresh.StepCount(); got != 25 {
		t.Fatalf("loaded step %d, want 25", got)
	}
}

// TestCheckpointRejectsMismatchedConfig ensures a checkpoint written for
// one input cannot be resumed with a different one.
func TestCheckpointRejectsMismatchedConfig(t *testing.T) {
	dir := t.TempDir()
	cfg := ckptConfig(dir, 20, 5)
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	other := cfg
	other.Dt = 0.009 // different numerical input, same checkpoint path
	if _, err := NewSolver(other); err == nil {
		t.Fatal("resumed from a checkpoint written with a different configuration")
	}
}

// TestPeriodicCheckpointMatchesFinalStep verifies Run commits exactly at
// the configured interval and at completion.
func TestPeriodicCheckpointMatchesFinalStep(t *testing.T) {
	dir := t.TempDir()
	cfg := ckptConfig(dir, 30, 7) // 30 is not a multiple of 7
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	fresh, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if got := fresh.StepCount(); got != 30 {
		t.Fatalf("final checkpoint at step %d, want 30", got)
	}
}
