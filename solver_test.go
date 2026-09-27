package solver

import (
	"errors"
	"math"
	"testing"
)

func baseConfig() Config {
	return Config{
		Lx: 10, Ly: 10, Nx: 60, Ny: 60,
		C: 1, TMax: 6, CFLFactor: 0.8, Workers: 1,
		Boundary: Boundaries{
			West: BoundaryFixed, East: BoundaryFixed,
			South: BoundaryFixed, North: BoundaryFixed,
		},
		Initial: InitialConditions{
			Displacement: GaussianPulse(30, 30, 4, 1),
			Velocity:     ZeroField(),
		},
	}
}

func TestMissingInitialConditionsRejected(t *testing.T) {
	cfg := baseConfig()
	cfg.Initial.Displacement = nil
	if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("missing displacement: want ErrInvalidConfig, got %v", err)
	}

	cfg = baseConfig()
	cfg.Initial.Velocity = nil
	if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("missing velocity: want ErrInvalidConfig, got %v", err)
	}
}

func TestInvalidGeometryRejected(t *testing.T) {
	cases := map[string]func(*Config){
		"zero Lx":        func(c *Config) { c.Lx = 0 },
		"negative Ly":    func(c *Config) { c.Ly = -1 },
		"tiny Nx":        func(c *Config) { c.Nx = 1 },
		"zero speed":     func(c *Config) { c.C = 0 },
		"zero tmax":      func(c *Config) { c.TMax = 0 },
		"bad workers":    func(c *Config) { c.Workers = 0 },
		"bad cfl factor": func(c *Config) { c.CFLFactor = 1.2 },
	}
	for name, mut := range cases {
		t.Run(name, func(t *testing.T) {
			cfg := baseConfig()
			mut(&cfg)
			if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
				t.Fatalf("want ErrInvalidConfig, got %v", err)
			}
		})
	}
}

func TestCFLViolationRejected(t *testing.T) {
	// dt deliberately larger than the CFL limit.
	cfg := baseConfig()
	cfg.Dt = 0.5 // dx ~= 0.1667 => dt must be < ~0.118
	cfg.CFLFactor = 0.8
	_, err := NewSolver(cfg)
	if !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("unstable dt: want ErrInvalidConfig, got %v", err)
	}

	// Even a dt at/above the exact limit with an aggressive factor is rejected.
	cfg = baseConfig()
	cfg.CFLFactor = 1.0
	cfg.Nx, cfg.Ny = 40, 40
	dx := cfg.Lx / float64(cfg.Nx)
	cfg.Dt = dx / cfg.C * 1.01 // courant = 1.01/sqrt(2) ... use big dt
	cfg.Dt = dx / cfg.C * 2.0
	if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("unstable dt: want rejection, got %v", err)
	}
}

func TestContradictoryInitialConditionsRejected(t *testing.T) {
	cfg := baseConfig()
	// Nonzero displacement right on a fixed wall.
	cfg.Initial.Displacement = func(i, j int) float64 {
		if i == 0 {
			return 0.5
		}
		return 0
	}
	if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("nonzero displacement on fixed edge: want rejection, got %v", err)
	}

	cfg = baseConfig()
	cfg.Initial.Velocity = func(i, j int) float64 {
		if j == cfg.Ny {
			return 1
		}
		return 0
	}
	if _, err := NewSolver(cfg); !errors.Is(err, ErrInvalidConfig) {
		t.Fatalf("nonzero velocity on fixed edge: want rejection, got %v", err)
	}
}

func TestSolutionStaysBounded(t *testing.T) {
	cfg := baseConfig()
	s, err := NewSolver(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if err := s.Run(); err != nil {
		t.Fatal(err)
	}
	d := s.Diagnostics()
	if d.MaxAbsU > 2.0 || math.IsNaN(d.MaxAbsU) || math.IsInf(d.MaxAbsU, 0) {
		t.Fatalf("field unbounded: max|u|=%v", d.MaxAbsU)
	}
	if s.StepCount() != s.TargetSteps() {
		t.Fatalf("steps=%d target=%d", s.StepCount(), s.TargetSteps())
	}
}
