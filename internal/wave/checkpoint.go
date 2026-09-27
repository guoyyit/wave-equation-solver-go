package wave

import (
	"bytes"
	"encoding/gob"
	"errors"
	"fmt"
	"hash/fnv"
	"io/fs"
	"os"
)

// checkpointMagic identifies the checkpoint format version.
const checkpointMagic = "wave-equation-solver-go/checkpoint/v1"

// checkpointFile is the on-disk, committed state of a run. Two time
// levels suffice to resume the leapfrog scheme exactly.
type checkpointFile struct {
	Magic       string
	Fingerprint uint64
	Step        int
	Prev        []float64
	Curr        []float64
}

// fingerprint hashes every configuration field that influences the
// numerical result. Fields that only affect scheduling (Workers) or
// checkpoint logistics (CheckpointPath, CheckpointEvery) are excluded,
// so changing them never invalidates a committed checkpoint.
func (c *Config) fingerprint() uint64 {
	numeric := struct {
		Nx, Ny     int
		Lx, Ly     float64
		C          float64
		Dt         float64
		Steps      int
		Boundaries Boundaries
		Init       InitialCondition
		Seed       int64
	}{c.Nx, c.Ny, c.Lx, c.Ly, c.C, c.Dt, c.Steps, c.Boundaries, c.Init, c.Seed}
	var buf bytes.Buffer
	if err := gob.NewEncoder(&buf).Encode(numeric); err != nil {
		panic(fmt.Errorf("wave: fingerprint encode: %w", err))
	}
	h := fnv.New64a()
	h.Write(buf.Bytes())
	return h.Sum64()
}

// saveCheckpoint atomically commits the solver state to
// cfg.CheckpointPath: the payload is written to a temporary file, synced,
// then renamed over the previous checkpoint. A crash at any point leaves
// either the old or the new committed checkpoint, never a partial one.
func (s *Solver) saveCheckpoint() error {
	cp := checkpointFile{
		Magic:       checkpointMagic,
		Fingerprint: s.cfg.fingerprint(),
		Step:        s.step,
		Prev:        s.prev.Data,
		Curr:        s.curr.Data,
	}
	var buf bytes.Buffer
	if err := gob.NewEncoder(&buf).Encode(cp); err != nil {
		return fmt.Errorf("wave: encode checkpoint: %w", err)
	}
	tmp := s.cfg.CheckpointPath + ".tmp"
	f, err := os.OpenFile(tmp, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
	if err != nil {
		return fmt.Errorf("wave: create checkpoint temp file: %w", err)
	}
	if _, err := f.Write(buf.Bytes()); err != nil {
		f.Close()
		return fmt.Errorf("wave: write checkpoint: %w", err)
	}
	if err := f.Sync(); err != nil {
		f.Close()
		return fmt.Errorf("wave: sync checkpoint: %w", err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("wave: close checkpoint: %w", err)
	}
	if err := os.Rename(tmp, s.cfg.CheckpointPath); err != nil {
		return fmt.Errorf("wave: commit checkpoint: %w", err)
	}
	return nil
}

// loadCheckpoint restores the newest committed checkpoint, if any.
// It reports whether a checkpoint was found. A checkpoint whose
// fingerprint does not match the current configuration is an error:
// resuming with different input would silently corrupt results.
func (s *Solver) loadCheckpoint() (bool, error) {
	raw, err := os.ReadFile(s.cfg.CheckpointPath)
	if errors.Is(err, fs.ErrNotExist) {
		return false, nil
	}
	if err != nil {
		return false, fmt.Errorf("wave: read checkpoint: %w", err)
	}
	var cp checkpointFile
	if err := gob.NewDecoder(bytes.NewReader(raw)).Decode(&cp); err != nil {
		return false, fmt.Errorf("wave: decode checkpoint: %w", err)
	}
	if cp.Magic != checkpointMagic {
		return false, fmt.Errorf("wave: checkpoint %q has unknown format %q",
			s.cfg.CheckpointPath, cp.Magic)
	}
	if want := s.cfg.fingerprint(); cp.Fingerprint != want {
		return false, fmt.Errorf("wave: checkpoint %q was created with a different configuration "+
			"(fingerprint %d, want %d); delete it or restore the original input",
			s.cfg.CheckpointPath, cp.Fingerprint, want)
	}
	n := s.cfg.Nx * s.cfg.Ny
	if len(cp.Prev) != n || len(cp.Curr) != n {
		return false, fmt.Errorf("wave: checkpoint %q is corrupt: field lengths %d/%d, want %d",
			s.cfg.CheckpointPath, len(cp.Prev), len(cp.Curr), n)
	}
	if cp.Step < 0 || cp.Step > s.cfg.Steps {
		return false, fmt.Errorf("wave: checkpoint %q is corrupt: step %d outside [0, %d]",
			s.cfg.CheckpointPath, cp.Step, s.cfg.Steps)
	}
	copy(s.prev.Data, cp.Prev)
	copy(s.curr.Data, cp.Curr)
	s.step = cp.Step
	return true, nil
}
