package wave

import (
	"fmt"
	"runtime"
	"sync"
)

// Solver advances the 2D scalar wave equation
//
//	u_tt = C^2 (u_xx + u_yy)
//
// with the explicit second-order leapfrog (centered-time, centered-space)
// scheme on a regular grid. Three time levels are kept so the scheme and
// the first-order Mur absorbing boundary only need previous and current
// fields.
type Solver struct {
	cfg  Config
	prev *Grid // u^{n-1}
	curr *Grid // u^{n}
	next *Grid // u^{n+1} (scratch)
	step int   // number of completed steps
}

// NewSolver validates cfg, builds the initial time levels and restores
// the newest committed checkpoint when cfg.CheckpointPath is set.
func NewSolver(cfg Config) (*Solver, error) {
	if err := cfg.Validate(); err != nil {
		return nil, err
	}
	s := &Solver{
		cfg:  cfg,
		prev: NewGrid(cfg.Nx, cfg.Ny),
		curr: NewGrid(cfg.Nx, cfg.Ny),
		next: NewGrid(cfg.Nx, cfg.Ny),
	}
	s.initLevels()
	if cfg.CheckpointPath != "" {
		found, err := s.loadCheckpoint()
		if err != nil {
			return nil, err
		}
		_ = found // absence of a checkpoint simply means "start from t=0"
	}
	return s, nil
}

// initLevels builds u^0 and u^{-1} from the initial condition.
// u^{-1} uses the Taylor start u^{-1} = u^0 - Dt*v^0 + (Dt^2/2) C^2 lap(u^0),
// which keeps the leapfrog scheme second-order accurate from step one.
func (s *Solver) initLevels() {
	c := &s.cfg
	if !c.Init.DisplacementZero {
		copy(s.curr.Data, c.Init.Displacement)
	}
	// Fixed boundaries hold u = 0 for all times, including t = 0.
	s.zeroFixedEdges(s.curr)

	dt, cc := c.Dt, c.C
	dx, dy := c.Dx(), c.Dy()
	vel := c.Init.Velocity
	for j := 1; j < c.Ny-1; j++ {
		for i := 1; i < c.Nx-1; i++ {
			idx := j*c.Nx + i
			u := s.curr.Data
			lap := (u[idx+1]-2*u[idx]+u[idx-1])/(dx*dx) +
				(u[idx+c.Nx]-2*u[idx]+u[idx-c.Nx])/(dy*dy)
			v := 0.0
			if !c.Init.VelocityZero {
				v = vel[idx]
			}
			s.prev.Data[idx] = u[idx] - dt*v + 0.5*dt*dt*cc*cc*lap
		}
	}
	s.zeroFixedEdges(s.prev)
}

// Step advances the solution by one time step.
func (s *Solver) Step() {
	c := &s.cfg
	dt, cc := c.Dt, c.C
	dx, dy := c.Dx(), c.Dy()
	rx := (cc * dt / dx) * (cc * dt / dx)
	ry := (cc * dt / dy) * (cc * dt / dy)

	workers := c.Workers
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	if workers > c.Ny-2 {
		workers = c.Ny - 2
	}
	if workers < 1 {
		workers = 1
	}

	// Interior update, partitioned by rows. Each output cell depends only
	// on curr and prev, so the parallel result is bit-identical to serial.
	var wg sync.WaitGroup
	rows := c.Ny - 2
	chunk := (rows + workers - 1) / workers
	for w := 0; w < workers; w++ {
		j0 := 1 + w*chunk
		j1 := j0 + chunk
		if j1 > c.Ny-1 {
			j1 = c.Ny - 1
		}
		if j0 >= j1 {
			break
		}
		wg.Add(1)
		go func(j0, j1 int) {
			defer wg.Done()
			for j := j0; j < j1; j++ {
				base := j * c.Nx
				for i := 1; i < c.Nx-1; i++ {
					idx := base + i
					u := s.curr.Data
					s.next.Data[idx] = 2*u[idx] - s.prev.Data[idx] +
						rx*(u[idx+1]-2*u[idx]+u[idx-1]) +
						ry*(u[idx+c.Nx]-2*u[idx]+u[idx-c.Nx])
				}
			}
		}(j0, j1)
	}
	wg.Wait()

	s.applyBoundaries()
	s.prev, s.curr, s.next = s.curr, s.next, s.prev
	s.step++
}

// Run advances until cfg.Steps steps have completed in total, committing
// checkpoints as configured. A solver resumed from a checkpoint continues
// from the committed step and never re-advances completed work: calling
// Run on a finished solver is a no-op.
func (s *Solver) Run() error {
	for s.step < s.cfg.Steps {
		s.Step()
		if !s.finite() {
			return fmt.Errorf("wave: numerical divergence detected at step %d "+
				"(non-finite field value); check CFL constraint and input magnitudes", s.step)
		}
		if s.cfg.CheckpointPath != "" && s.step%s.cfg.CheckpointEvery == 0 && s.step < s.cfg.Steps {
			if err := s.saveCheckpoint(); err != nil {
				return fmt.Errorf("wave: checkpoint at step %d: %w", s.step, err)
			}
		}
	}
	if s.cfg.CheckpointPath != "" {
		if err := s.saveCheckpoint(); err != nil {
			return fmt.Errorf("wave: final checkpoint: %w", err)
		}
	}
	return nil
}

// StepCount returns the number of completed steps.
func (s *Solver) StepCount() int { return s.step }

// Snapshot returns a copy of the current displacement field u^n.
func (s *Solver) Snapshot() *Grid { return s.curr.Clone() }
