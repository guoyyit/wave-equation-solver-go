// Package wave implements a 2D scalar wave equation solver with
// checkpointing, absorbing boundaries and parallel stepping.
package wave

import (
	"errors"
	"fmt"
	"math"
)

// BoundaryType selects the condition applied on one edge of the rectangle.
type BoundaryType int

const (
	// BoundaryFixed is a Dirichlet boundary: u = 0 on the edge.
	BoundaryFixed BoundaryType = iota
	// BoundaryAbsorbing is a first-order Mur absorbing boundary.
	BoundaryAbsorbing
)

// String renders the boundary type for diagnostics and config dumps.
func (b BoundaryType) String() string {
	switch b {
	case BoundaryFixed:
		return "fixed"
	case BoundaryAbsorbing:
		return "absorbing"
	default:
		return "unknown"
	}
}

// validBoundary reports whether b is a known boundary type.
func validBoundary(b BoundaryType) bool {
	return b == BoundaryFixed || b == BoundaryAbsorbing
}

// Boundaries assigns a boundary type to each edge of the rectangle.
type Boundaries struct {
	Left, Right, Top, Bottom BoundaryType
}

// InitialCondition carries the initial displacement and velocity fields.
// A nil slice is accepted only when the matching Zero flag is set, which
// documents the intent "this field is identically zero". Any other
// combination is rejected by Validate instead of silently defaulting.
type InitialCondition struct {
	// Displacement holds u(x, y, 0), row-major, length Nx*Ny.
	Displacement []float64
	// DisplacementZero declares Displacement is intentionally zero.
	DisplacementZero bool
	// Velocity holds du/dt(x, y, 0), row-major, length Nx*Ny.
	Velocity []float64
	// VelocityZero declares Velocity is intentionally zero.
	VelocityZero bool
}

// Config fully describes one simulation run.
type Config struct {
	Nx, Ny     int     // grid points per dimension (>= 2)
	Lx, Ly     float64 // physical extents (m, > 0)
	C          float64 // wave speed (m/s, > 0)
	Dt         float64 // time step (s, > 0); must satisfy the CFL bound
	Steps      int     // number of time steps to advance (>= 1)
	Boundaries Boundaries
	Init       InitialCondition
	Seed       int64 // reserved for stochastic source extensions
	Workers    int   // parallel workers; <= 0 means runtime.NumCPU()
	// CheckpointPath enables checkpointing when non-empty.
	CheckpointPath string
	// CheckpointEvery commits a checkpoint every this many steps (>= 1).
	CheckpointEvery int
}

// Dx returns the grid spacing in x.
func (c *Config) Dx() float64 { return c.Lx / float64(c.Nx-1) }

// Dy returns the grid spacing in y.
func (c *Config) Dy() float64 { return c.Ly / float64(c.Ny-1) }

// Courant returns r = C*Dt*sqrt(1/Dx^2 + 1/Dy^2). The explicit 5-point
// scheme used by this package is stable iff r <= 1.
func (c *Config) Courant() float64 {
	dx, dy := c.Dx(), c.Dy()
	return c.C * c.Dt * math.Sqrt(1/(dx*dx)+1/(dy*dy))
}

// MaxDt returns the largest time step allowed by the CFL condition.
func (c *Config) MaxDt() float64 {
	dx, dy := c.Dx(), c.Dy()
	return 1 / (c.C * math.Sqrt(1/(dx*dx)+1/(dy*dy)))
}

// Validate checks the configuration for consistency and stability.
// It returns the first problem found; nothing is silently defaulted.
func (c *Config) Validate() error {
	if c.Nx < 2 || c.Ny < 2 {
		return fmt.Errorf("wave: grid too small: Nx=%d Ny=%d, both must be >= 2", c.Nx, c.Ny)
	}
	if !(c.Lx > 0) || !(c.Ly > 0) {
		return fmt.Errorf("wave: extents must be positive: Lx=%v Ly=%v", c.Lx, c.Ly)
	}
	if !(c.C > 0) {
		return fmt.Errorf("wave: wave speed must be positive: C=%v", c.C)
	}
	if !(c.Dt > 0) {
		return fmt.Errorf("wave: time step must be positive: Dt=%v", c.Dt)
	}
	if c.Steps < 1 {
		return fmt.Errorf("wave: Steps must be >= 1, got %d", c.Steps)
	}
	if r := c.Courant(); r > 1 {
		return fmt.Errorf("wave: CFL stability violated: r=%.6f > 1 (Dt=%.6g, max stable Dt=%.6g)",
			r, c.Dt, c.MaxDt())
	}
	if b := c.Boundaries; !validBoundary(b.Left) || !validBoundary(b.Right) ||
		!validBoundary(b.Top) || !validBoundary(b.Bottom) {
		return fmt.Errorf("wave: unknown boundary type in %+v", c.Boundaries)
	}
	if err := c.validateInitial(); err != nil {
		return err
	}
	if c.CheckpointPath != "" && c.CheckpointEvery < 1 {
		return fmt.Errorf("wave: CheckpointEvery must be >= 1 when CheckpointPath is set, got %d",
			c.CheckpointEvery)
	}
	return nil
}

// validateInitial enforces the initial-condition contract: each of
// displacement and velocity must be either fully provided with the right
// length, or explicitly declared zero. Missing or contradictory input is
// an error, never a silent default.
func (c *Config) validateInitial() error {
	n := c.Nx * c.Ny
	check := func(name string, field []float64, zero bool) error {
		switch {
		case field == nil && !zero:
			return fmt.Errorf("wave: initial %s missing: provide a field of length %d or set its Zero flag",
				name, n)
		case field != nil && zero:
			return fmt.Errorf("wave: initial %s contradictory: field provided but Zero flag is set", name)
		case field != nil && len(field) != n:
			return fmt.Errorf("wave: initial %s has length %d, want %d (Nx*Ny)", name, len(field), n)
		}
		return nil
	}
	if err := check("displacement", c.Init.Displacement, c.Init.DisplacementZero); err != nil {
		return err
	}
	if err := check("velocity", c.Init.Velocity, c.Init.VelocityZero); err != nil {
		return err
	}
	if c.Init.DisplacementZero && c.Init.VelocityZero {
		return errors.New("wave: trivial initial condition: displacement and velocity both zero")
	}
	return nil
}
