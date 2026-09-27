package wave

import (
	"strings"
	"testing"
)

// validConfig returns a minimal valid configuration for test mutation.
func validConfig() Config {
	nx, ny := 16, 16
	disp := make([]float64, nx*ny)
	disp[8*nx+8] = 1
	return Config{
		Nx: nx, Ny: ny, Lx: 1, Ly: 1, C: 1, Dt: 0.005, Steps: 10,
		Init: InitialCondition{Displacement: disp, VelocityZero: true},
	}
}

func TestValidateAcceptsValidConfig(t *testing.T) {
	cfg := validConfig()
	if err := cfg.Validate(); err != nil {
		t.Fatalf("valid config rejected: %v", err)
	}
}

func TestValidateRejectsCFLViolation(t *testing.T) {
	cfg := validConfig()
	cfg.Dt = cfg.MaxDt() * 1.01
	err := cfg.Validate()
	if err == nil {
		t.Fatal("CFL-violating Dt accepted")
	}
	if !strings.Contains(err.Error(), "CFL") {
		t.Fatalf("rejection does not explain the stability reason: %v", err)
	}
}

func TestValidateAcceptsCFLBoundary(t *testing.T) {
	cfg := validConfig()
	cfg.Dt = cfg.MaxDt() // r == 1 exactly is still stable
	if err := cfg.Validate(); err != nil {
		t.Fatalf("r == 1 should be stable, got: %v", err)
	}
}

func TestValidateRejectsMissingInitialConditions(t *testing.T) {
	cfg := validConfig()
	cfg.Init = InitialCondition{DisplacementZero: true} // velocity neither given nor declared zero
	err := cfg.Validate()
	if err == nil || !strings.Contains(err.Error(), "velocity") {
		t.Fatalf("missing velocity not reported clearly: %v", err)
	}

	cfg.Init = InitialCondition{VelocityZero: true} // displacement missing
	err = cfg.Validate()
	if err == nil || !strings.Contains(err.Error(), "displacement") {
		t.Fatalf("missing displacement not reported clearly: %v", err)
	}
}

func TestValidateRejectsContradictoryInitialConditions(t *testing.T) {
	cfg := validConfig()
	cfg.Init.Velocity = make([]float64, cfg.Nx*cfg.Ny)
	cfg.Init.VelocityZero = true // field provided AND declared zero
	err := cfg.Validate()
	if err == nil || !strings.Contains(err.Error(), "contradictory") {
		t.Fatalf("contradictory velocity not reported: %v", err)
	}
}

func TestValidateRejectsWrongInitialFieldLength(t *testing.T) {
	cfg := validConfig()
	cfg.Init.Displacement = make([]float64, cfg.Nx*cfg.Ny-1)
	err := cfg.Validate()
	if err == nil || !strings.Contains(err.Error(), "length") {
		t.Fatalf("wrong field length not reported: %v", err)
	}
}

func TestValidateRejectsTrivialInitialCondition(t *testing.T) {
	cfg := validConfig()
	cfg.Init = InitialCondition{DisplacementZero: true, VelocityZero: true}
	if err := cfg.Validate(); err == nil {
		t.Fatal("all-zero initial condition accepted")
	}
}

func TestValidateRejectsBadGeometryAndSteps(t *testing.T) {
	for _, mutate := range []func(*Config){
		func(c *Config) { c.Nx = 1 },
		func(c *Config) { c.Ny = 0 },
		func(c *Config) { c.Lx = 0 },
		func(c *Config) { c.Ly = -1 },
		func(c *Config) { c.C = 0 },
		func(c *Config) { c.Dt = 0 },
		func(c *Config) { c.Steps = 0 },
	} {
		cfg := validConfig()
		mutate(&cfg)
		if err := cfg.Validate(); err == nil {
			t.Fatalf("bad config accepted: %+v", cfg)
		}
	}
}

func TestValidateRejectsUnknownBoundaryType(t *testing.T) {
	cfg := validConfig()
	cfg.Boundaries.Left = BoundaryType(42)
	if err := cfg.Validate(); err == nil {
		t.Fatal("unknown boundary type accepted")
	}
}

func TestValidateRejectsBadCheckpointInterval(t *testing.T) {
	cfg := validConfig()
	cfg.CheckpointPath = t.TempDir() + "/ckpt.bin"
	cfg.CheckpointEvery = 0
	if err := cfg.Validate(); err == nil {
		t.Fatal("CheckpointEvery = 0 with checkpointing enabled accepted")
	}
}

// TestNewSolverRejectsInvalidConfig ensures the constructor refuses
// unstable or inconsistent input instead of solving with it.
func TestNewSolverRejectsInvalidConfig(t *testing.T) {
	cfg := validConfig()
	cfg.Dt = 1 // grossly CFL-violating
	if _, err := NewSolver(cfg); err == nil {
		t.Fatal("NewSolver accepted a CFL-violating configuration")
	}
}
