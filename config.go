package solver

import (
	"errors"
	"fmt"
	"math"
)

// BoundaryKind identifies the boundary condition applied on one side
// of the rectangular domain.
type BoundaryKind int

const (
	// BoundaryFixed is the homogeneous Dirichlet boundary: u = 0.
	BoundaryFixed BoundaryKind = iota
	// BoundaryAbsorbing is a first-order Mur absorbing boundary.
	BoundaryAbsorbing
)

// String returns a human readable name of the boundary kind.
func (k BoundaryKind) String() string {
	switch k {
	case BoundaryFixed:
		return "fixed"
	case BoundaryAbsorbing:
		return "absorbing"
	default:
		return fmt.Sprintf("boundary(%d)", int(k))
	}
}

// FieldFunc evaluates a scalar field at grid point (i, j).
type FieldFunc func(i, j int) float64

// Boundaries selects the boundary kind on each of the four sides.
type Boundaries struct {
	West  BoundaryKind
	East  BoundaryKind
	South BoundaryKind
	North BoundaryKind
}

// InitialConditions holds the initial displacement and velocity.
// Both fields must be provided explicitly; use ZeroField() when a zero
// initial velocity (or displacement) is intended. A nil value is never
// silently replaced by a default.
type InitialConditions struct {
	Displacement FieldFunc
	Velocity     FieldFunc
}

// Config describes a 2D scalar wave equation run:
//
//	d2u/dt2 = C^2 (d2u/dx2 + dy2)
//
// on [0, Lx] x [0, Ly] discretized with (Nx+1) x (Ny+1) grid points.
type Config struct {
	Lx       float64 // domain length in x
	Ly       float64 // domain length in y
	Nx       int     // number of grid intervals in x (grid has Nx+1 points)
	Ny       int     // number of grid intervals in y (grid has Ny+1 points)
	C        float64 // wave speed
	Dt       float64 // time step; if zero it is derived with CFLFactor
	TMax     float64 // final physical time
	Boundary Boundaries
	Initial  InitialConditions

	// CFLFactor is the safety multiplier in (0, 1] applied to the CFL
	// limit when Dt is derived automatically.
	CFLFactor float64

	// Workers is the number of goroutines used for interior updates.
	// Values <= 1 select the single-threaded code path.
	Workers int

	// RandomSeed makes memory-constructed random inputs reproducible.
	RandomSeed int64

	// CheckpointDir is where committed checkpoints are written. Empty
	// disables checkpointing.
	CheckpointDir string

	// CheckpointInterval commits a checkpoint every N advanced steps.
	CheckpointInterval int
}

// ErrInvalidConfig is returned by Validate when the configuration is
// incomplete, contradictory or numerically unstable.
var ErrInvalidConfig = errors.New("invalid configuration")

// cflMargin rejects Courant numbers within this distance of 1.0 to keep
// a safety margin against roundoff-driven growth.
const cflMargin = 1e-9

// isFinite reports whether x is neither NaN nor infinite.
func isFinite(x float64) bool {
	return !math.IsNaN(x) && !math.IsInf(x, 0)
}

// Validate checks completeness and numerical stability of the config.
// It returns a wrapped ErrInvalidConfig explaining the rejected reason.
func (c *Config) Validate() error {
	bad := func(format string, args ...any) error {
		return fmt.Errorf("%w: %s", ErrInvalidConfig, fmt.Sprintf(format, args...))
	}

	if !(c.Lx > 0) || !isFinite(c.Lx) {
		return bad("Lx must be a positive finite number, got %v", c.Lx)
	}
	if !(c.Ly > 0) || !isFinite(c.Ly) {
		return bad("Ly must be a positive finite number, got %v", c.Ly)
	}
	if c.Nx < 2 {
		return bad("Nx must be >= 2 grid intervals, got %d", c.Nx)
	}
	if c.Ny < 2 {
		return bad("Ny must be >= 2 grid intervals, got %d", c.Ny)
	}
	if !(c.C > 0) || !isFinite(c.C) {
		return bad("wave speed C must be a positive finite number, got %v", c.C)
	}
	if !(c.TMax > 0) || !isFinite(c.TMax) {
		return bad("TMax must be a positive finite number, got %v", c.TMax)
	}

	validKind := func(k BoundaryKind) bool {
		return k == BoundaryFixed || k == BoundaryAbsorbing
	}
	if !validKind(c.Boundary.West) || !validKind(c.Boundary.East) ||
		!validKind(c.Boundary.South) || !validKind(c.Boundary.North) {
		return bad("all four boundary kinds must be Fixed or Absorbing")
	}

	if c.Initial.Displacement == nil {
		return bad("initial displacement is missing: provide a FieldFunc or ZeroField()")
	}
	if c.Initial.Velocity == nil {
		return bad("initial velocity is missing: provide a FieldFunc or ZeroField()")
	}

	if c.Workers < 1 {
		return bad("Workers must be >= 1, got %d", c.Workers)
	}
	if c.CheckpointDir != "" && c.CheckpointInterval < 1 {
		return bad("CheckpointInterval must be >= 1 when CheckpointDir is set")
	}

	// Resolve the time step.
	if c.Dt < 0 || (c.Dt != 0 && !isFinite(c.Dt)) {
		return bad("Dt must be positive finite or zero (derive from CFL), got %v", c.Dt)
	}
	if c.Dt == 0 {
		if c.CFLFactor <= 0 || c.CFLFactor > 1 {
			return bad("when Dt is zero, CFLFactor must be in (0,1], got %v", c.CFLFactor)
		}
		dt0 := c.CFLFactor * c.CFLStep()
		// Choose an integer number of steps and derive dt = TMax/n so the
		// run lands exactly on TMax (required for deterministic resume and
		// probe comparisons). dt <= dt0, hence the CFL condition still holds.
		n := int(math.Ceil(c.TMax / dt0))
		if n < 1 {
			n = 1
		}
		c.Dt = c.TMax / float64(n)
	}
	if c.Dt <= 0 {
		return bad("derived time step is not positive, got %v", c.Dt)
	}

	// Explicit CFL stability check for the leapfrog scheme.
	courant := c.courant(c.Dt)
	if courant >= 1.0-cflMargin {
		return bad("CFL stability violated: courant number %.12g >= 1 (dt=%v, dx=%v, dy=%v, C=%v); reduce dt or refine the grid",
			courant, c.Dt, c.Lx/float64(c.Nx), c.Ly/float64(c.Ny), c.C)
	}
	if c.Dt > c.TMax {
		return bad("Dt %v exceeds TMax %v: at least one step is required", c.Dt, c.TMax)
	}

	// TMax must line up with an integer number of steps so resuming and
	// probe-time comparisons are unambiguous. With an explicit Dt this
	// means TMax must be a whole multiple of Dt; with Dt left at zero the
	// step size is derived above to guarantee this automatically.
	ratio := c.TMax / c.Dt
	nSteps := int(math.Round(ratio))
	if math.Abs(ratio-float64(nSteps)) > 1e-9*math.Max(1, ratio) {
		return bad("TMax/Dt must be an integer number of steps, got %.12g (set Dt=0 to derive a fitting step)", ratio)
	}

	return nil
}

// StepCount returns the number of time steps covering [0, TMax].
func (c *Config) StepCount() int {
	return int(math.Round(c.TMax / c.Dt))
}

// CFLStep returns the largest stable time step for the spatial grid:
//
//	dt_max = 1 / (C * sqrt(1/dx^2 + 1/dy^2)).
func (c *Config) CFLStep() float64 {
	dx := c.Lx / float64(c.Nx)
	dy := c.Ly / float64(c.Ny)
	return 1.0 / (c.C * math.Sqrt(1.0/(dx*dx)+1.0/(dy*dy)))
}

// courant returns the Courant number C*dt*sqrt(1/dx^2+1/dy^2).
func (c *Config) courant(dt float64) float64 {
	dx := c.Lx / float64(c.Nx)
	dy := c.Ly / float64(c.Ny)
	return c.C * dt * math.Sqrt(1.0/(dx*dx)+1.0/(dy*dy))
}
