package solver

import "fmt"

// Solver advances the 2D scalar wave equation
//
//	d2u/dt2 = C^2 (d2u/dx2 + d2u/dy2)
//
// with a leapfrog finite difference scheme on a uniform rectangular
// grid. Arrays are flat and indexed as i*ny + j with i in [0,nx),
// j in [0,ny).
type Solver struct {
	cfg Config
	nx  int // grid points in x (Nx+1)
	ny  int // grid points in y (Ny+1)
	dx  float64
	dy  float64
	rx  float64 // C*dt/dx
	ry  float64 // C*dt/dy
	rx2 float64 // rx*rx
	ry2 float64 // ry*ry

	// uPrev = u^(n-1), uCurr = u^n, uNext = u^(n+1).
	uPrev []float64
	uCurr []float64
	uNext []float64

	step int
	time float64

	lastCheckpointStep int
}

// prepareConfig validates cfg and returns a normalized copy with the time
// step resolved (e.g. derived from the CFL factor) so downstream code and
// checkpoint hashing always see the concrete Dt.
func prepareConfig(cfg Config) (Config, error) {
	if err := cfg.Validate(); err != nil {
		return Config{}, err
	}
	return cfg, nil
}

// NewSolver constructs a solver and initializes the first two time levels
// from the initial displacement u0 and velocity v0. It returns an error
// if the configuration is invalid or the initial conditions are missing.
func NewSolver(cfg Config) (*Solver, error) {
	pc, err := prepareConfig(cfg)
	if err != nil {
		return nil, err
	}
	cfg = pc

	nx := cfg.Nx + 1
	ny := cfg.Ny + 1
	n := nx * ny
	dx := cfg.Lx / float64(cfg.Nx)
	dy := cfg.Ly / float64(cfg.Ny)
	rx := cfg.C * cfg.Dt / dx
	ry := cfg.C * cfg.Dt / dy

	s := &Solver{
		cfg:   cfg,
		nx:    nx,
		ny:    ny,
		dx:    dx,
		dy:    dy,
		rx:    rx,
		ry:    ry,
		rx2:   rx * rx,
		ry2:   ry * ry,
		uPrev: make([]float64, n),
		uCurr: make([]float64, n),
		uNext: make([]float64, n),
	}

	// u^0 = u0.
	for i := 0; i < nx; i++ {
		for j := 0; j < ny; j++ {
			v := cfg.Initial.Displacement(i, j)
			s.uCurr[i*ny+j] = v
		}
	}
	// Contradictory initial conditions: a fixed (u=0) boundary must not
	// carry a nonzero prescribed displacement or velocity. Reject instead
	// of silently clamping it.
	if err := s.checkInitialConsistency(); err != nil {
		return nil, err
	}
	// Enforce boundary consistency on u^0.
	s.zeroFixedEdges(s.uCurr)

	// First step via Taylor expansion:
	// u^1 = u^0 + dt*v0 + 0.5*dt^2*C^2*laplacian(u^0).
	s.computeFirstStep()
	s.zeroFixedEdges(s.uNext)
	if s.hasAbsorbing() {
		s.applyFirstMurEdges()
	}

	s.uPrev = s.uCurr
	s.uCurr = s.uNext
	s.uNext = make([]float64, n)
	s.step = 1
	s.time = cfg.Dt

	return s, nil
}

// Step advances the solution by exactly one time step.
func (s *Solver) Step() error {
	if s.step >= s.cfg.StepCount() {
		return fmt.Errorf("already at final step %d (t=%v)", s.step, s.time)
	}

	// Clear uNext so corner/edge handling never reads stale values.
	for i := range s.uNext {
		s.uNext[i] = 0
	}
	s.updateInterior()
	s.applyBoundaries()
	if err := s.guardFinite(); err != nil {
		return err
	}

	s.uPrev, s.uCurr, s.uNext = s.uCurr, s.uNext, s.uPrev
	s.step++
	s.time = float64(s.step) * s.cfg.Dt

	if s.cfg.CheckpointDir != "" && s.step%s.cfg.CheckpointInterval == 0 {
		if err := s.commitCheckpoint(); err != nil {
			return fmt.Errorf("checkpoint at step %d: %w", s.step, err)
		}
		s.lastCheckpointStep = s.step
	}
	return nil
}

// Run advances until physical time TMax, committing checkpoints along
// the way when configured.
func (s *Solver) Run() error {
	total := s.cfg.StepCount()
	for s.step < total {
		if err := s.Step(); err != nil {
			return err
		}
	}
	return nil
}

// Resume loads the latest committed checkpoint from dir and rebuilds a
// solver positioned exactly at the recorded step. The recorded step is
// never re-advanced; a subsequent Run/Step continues from step+1.
func Resume(dir string, cfg Config) (*Solver, error) {
	data, err := loadLatestCheckpoint(dir)
	if err != nil {
		return nil, fmt.Errorf("resume: %w", err)
	}
	pc, err := prepareConfig(cfg)
	if err != nil {
		return nil, err
	}
	cfg = pc
	if data.ConfigHash != configHash(cfg) {
		return nil, fmt.Errorf("%w: checkpoint configuration does not match the supplied config", ErrInvalidConfig)
	}
	if data.Nx != cfg.Nx+1 || data.Ny != cfg.Ny+1 {
		return nil, fmt.Errorf("%w: checkpoint grid %dx%d does not match config", ErrInvalidConfig, data.Nx, data.Ny)
	}

	nx := cfg.Nx + 1
	ny := cfg.Ny + 1
	dx := cfg.Lx / float64(cfg.Nx)
	dy := cfg.Ly / float64(cfg.Ny)
	rx := cfg.C * cfg.Dt / dx
	ry := cfg.C * cfg.Dt / dy
	s := &Solver{
		cfg:                cfg,
		nx:                 nx,
		ny:                 ny,
		dx:                 dx,
		dy:                 dy,
		rx:                 rx,
		ry:                 ry,
		rx2:                rx * rx,
		ry2:                ry * ry,
		uPrev:              append([]float64(nil), data.UPrev...),
		uCurr:              append([]float64(nil), data.UCurr...),
		uNext:              make([]float64, nx*ny),
		step:               data.Step,
		time:               data.Time,
		lastCheckpointStep: data.Step,
	}
	return s, nil
}

// Solution returns a copy of the current field u^n.
func (s *Solver) Solution() []float64 {
	out := make([]float64, len(s.uCurr))
	copy(out, s.uCurr)
	return out
}

// At returns u at grid indices (i, j).
func (s *Solver) At(i, j int) float64 {
	return s.uCurr[i*s.ny+j]
}

// StepCount returns the number of steps already advanced.
func (s *Solver) StepCount() int { return s.step }

// TargetSteps returns the total number of steps required to reach TMax.
func (s *Solver) TargetSteps() int { return s.cfg.StepCount() }

// Time returns the current physical time.
func (s *Solver) Time() float64 { return s.time }

// Diagnostics returns run diagnostics at the current state.
func (s *Solver) Diagnostics() Diagnostics {
	return s.ComputeDiagnostics()
}

// icTolerance is the allowed magnitude of initial values on fixed edges
// (generous tolerance for user-supplied analytic fields).
const icTolerance = 1e-8

// checkInitialConsistency rejects nonzero initial displacement or
// velocity prescribed on a fixed boundary, since that contradicts u=0
// there. Corners are always fixed.
func (s *Solver) checkInitialConsistency() error {
	nx, ny := s.nx, s.ny
	b := s.cfg.Boundary
	onFixed := func(i, j int) bool {
		// Only non-corner cells on a fixed edge are required to carry a
		// zero initial value. Corners meet two edges (which may differ in
		// kind) and are reconciled separately by the solver.
		isCorner := (i == 0 || i == nx-1) && (j == 0 || j == ny-1)
		if isCorner {
			return false
		}
		switch {
		case i == 0:
			return b.West == BoundaryFixed
		case i == nx-1:
			return b.East == BoundaryFixed
		case j == 0:
			return b.South == BoundaryFixed
		case j == ny-1:
			return b.North == BoundaryFixed
		default:
			return false
		}
	}
	for i := 0; i < nx; i++ {
		for j := 0; j < ny; j++ {
			if !onFixed(i, j) {
				continue
			}
			if v := s.uCurr[i*ny+j]; v > icTolerance || v < -icTolerance {
				return fmt.Errorf("%w: initial displacement %.6g at (%d,%d) contradicts fixed boundary (must be 0)",
					ErrInvalidConfig, v, i, j)
			}
			if v := s.cfg.Initial.Velocity(i, j); v > icTolerance || v < -icTolerance {
				return fmt.Errorf("%w: initial velocity %.6g at (%d,%d) contradicts fixed boundary (must be 0)",
					ErrInvalidConfig, v, i, j)
			}
		}
	}
	return nil
}

// hasAbsorbing reports whether at least one edge is absorbing.
func (s *Solver) hasAbsorbing() bool {
	b := s.cfg.Boundary
	return b.West == BoundaryAbsorbing || b.East == BoundaryAbsorbing ||
		b.South == BoundaryAbsorbing || b.North == BoundaryAbsorbing
}

// zeroFixedEdges sets fixed boundary cells of u to zero. Corners are
// always zero regardless of the neighboring edge kinds.
func (s *Solver) zeroFixedEdges(u []float64) {
	nx, ny := s.nx, s.ny
	b := s.cfg.Boundary
	u[0] = 0
	u[ny-1] = 0
	u[(nx-1)*ny] = 0
	u[(nx-1)*ny+ny-1] = 0
	zeroRow := func(base int) {
		for j := 1; j < ny-1; j++ {
			u[base+j] = 0
		}
	}
	if b.West == BoundaryFixed {
		zeroRow(0)
	}
	if b.East == BoundaryFixed {
		zeroRow((nx - 1) * ny)
	}
	if b.South == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			u[i*ny] = 0
		}
	}
	if b.North == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			u[i*ny+ny-1] = 0
		}
	}
}

// computeFirstStep fills s.uNext (u^1) from u0 and v0 using the second
// order Taylor start on interior points.
func (s *Solver) computeFirstStep() {
	nx, ny := s.nx, s.ny
	dt := s.cfg.Dt
	halfDt2C2 := 0.5 * dt * dt * s.cfg.C * s.cfg.C
	invDx2 := 1.0 / (s.dx * s.dx)
	invDy2 := 1.0 / (s.dy * s.dy)
	curr := s.uCurr
	for i := 1; i < nx-1; i++ {
		row := i * ny
		for j := 1; j < ny-1; j++ {
			idx := row + j
			lap := (curr[idx-ny]-2*curr[idx]+curr[idx+ny])*invDx2 +
				(curr[idx-1]-2*curr[idx]+curr[idx+1])*invDy2
			v0 := s.cfg.Initial.Velocity(i, j)
			s.uNext[idx] = curr[idx] + dt*v0 + halfDt2C2*lap
		}
	}
}

// applyFirstMurEdges sets u^1 on absorbing edges by applying the Mur
// update with u^0 as the old field. It is only used for the Taylor
// start; subsequent steps go through applyBoundaries.
func (s *Solver) applyFirstMurEdges() {
	// Temporarily arrange state so applyBoundaries sees prev/curr/next
	// exactly like a normal step: it only reads curr and next.
	oldPrev := s.uPrev
	s.uPrev = s.uCurr // not read by applyBoundaries, kept for clarity
	defer func() { s.uPrev = oldPrev }()
	s.applyBoundaries()
}

// ErrUnstable is returned when the update produces a non-finite or
// runaway value despite CFL validation (defensive backstop).
var ErrUnstable = fmt.Errorf("numerical instability detected")

// guardFinite turns numerical blow-ups into a hard error instead of
// letting NaN/Inf propagate. (CFL validation should make this
// unreachable; it is a defensive backstop.)
func (s *Solver) guardFinite() error {
	for _, v := range s.uNext {
		if v != v || v > 1e300 || v < -1e300 {
			return fmt.Errorf("%w: non-finite or runaway field value at step %d", ErrUnstable, s.step+1)
		}
	}
	return nil
}
