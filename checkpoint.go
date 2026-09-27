package solver

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
)

// checkpointData is the serialized solver state.
type checkpointData struct {
	Version    int       `json:"version"`
	Nx         int       `json:"nx"`
	Ny         int       `json:"ny"`
	Dx         float64   `json:"dx"`
	Dy         float64   `json:"dy"`
	Dt         float64   `json:"dt"`
	Step       int       `json:"step"`
	Time       float64   `json:"time"`
	UPrev      []float64 `json:"u_prev"`
	UCurr      []float64 `json:"u_curr"`
	ConfigHash string    `json:"config_hash"`
}

const checkpointVersion = 1
const checkpointName = "checkpoint-latest.json"

// configHash summarizes the parts of the configuration that must stay
// identical across a checkpoint/resume cycle. Checkpoint placement and
// worker count are intentionally excluded (they may change on restart).
func configHash(c Config) string {
	b := c.Boundary
	key := fmt.Sprintf("v1|%v|%v|%d|%d|%v|%v|%v|%d|%d|%d|%d|%d",
		c.Lx, c.Ly, c.Nx, c.Ny, c.C, c.Dt, c.TMax,
		b.West, b.East, b.South, b.North, c.RandomSeed)
	sum := sha256.Sum256([]byte(key))
	return hex.EncodeToString(sum[:])
}

// commitCheckpoint atomically writes the current solver state to the
// configured directory: data is written and synced to a temp file first,
// then renamed over the checkpoint name. A crash during writing therefore
// leaves either the previous complete checkpoint or the new one, never a
// torn file.
func (s *Solver) commitCheckpoint() error {
	if err := os.MkdirAll(s.cfg.CheckpointDir, 0o755); err != nil {
		return err
	}
	data := checkpointData{
		Version:    checkpointVersion,
		Nx:         s.nx,
		Ny:         s.ny,
		Dx:         s.dx,
		Dy:         s.dy,
		Dt:         s.cfg.Dt,
		Step:       s.step,
		Time:       s.time,
		UPrev:      s.uPrev,
		UCurr:      s.uCurr,
		ConfigHash: configHash(s.cfg),
	}
	buf, err := json.MarshalIndent(data, "", "  ")
	if err != nil {
		return err
	}

	dir := s.cfg.CheckpointDir
	final := filepath.Join(dir, checkpointName)
	tmp, err := os.CreateTemp(dir, ".checkpoint-*.tmp")
	if err != nil {
		return err
	}
	tmpName := tmp.Name()
	cleanup := func() { _ = os.Remove(tmpName) }

	if _, err := tmp.Write(buf); err != nil {
		_ = tmp.Close()
		cleanup()
		return err
	}
	if err := tmp.Sync(); err != nil {
		_ = tmp.Close()
		cleanup()
		return err
	}
	if err := tmp.Close(); err != nil {
		cleanup()
		return err
	}
	if err := os.Rename(tmpName, final); err != nil {
		cleanup()
		return err
	}
	// Best-effort directory sync so the rename itself is durable.
	if d, err := os.Open(dir); err == nil {
		_ = d.Sync()
		_ = d.Close()
	}
	return nil
}

// loadLatestCheckpoint reads the committed checkpoint from dir.
func loadLatestCheckpoint(dir string) (*checkpointData, error) {
	path := filepath.Join(dir, checkpointName)
	buf, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("no checkpoint found in %q", dir)
		}
		return nil, err
	}
	var data checkpointData
	if err := json.Unmarshal(buf, &data); err != nil {
		return nil, fmt.Errorf("checkpoint is corrupt: %w", err)
	}
	if data.Version != checkpointVersion {
		return nil, fmt.Errorf("unsupported checkpoint version %d", data.Version)
	}
	if len(data.UCurr) != data.Nx*data.Ny || len(data.UPrev) != data.Nx*data.Ny {
		return nil, fmt.Errorf("checkpoint field size does not match grid")
	}
	return &data, nil
}

// HasCheckpoint reports whether a committed checkpoint exists in dir.
func HasCheckpoint(dir string) bool {
	_, err := os.Stat(filepath.Join(dir, checkpointName))
	return err == nil
}

// ListTempFiles returns leftover uncommitted checkpoint temp files in dir.
// After a successful commit (and an atomic crash-safe protocol) this must
// be empty.
func ListTempFiles(dir string) []string {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil
	}
	var names []string
	for _, e := range entries {
		name := e.Name()
		if !e.IsDir() && filepath.Ext(name) == ".tmp" {
			names = append(names, filepath.Join(dir, name))
		}
	}
	return names
}
