package wave

import "math"

// Diagnostics aggregates per-run health information.
type Diagnostics struct {
	// Energy is the discrete total energy: kinetic (from the leapfrog
	// velocity estimate) plus elastic (from centered gradients), summed
	// over the cell area. Absorbing boundaries drain it; fixed boundaries
	// conserve it up to scheme accuracy.
	Energy float64
	// MaxAbs is max |u| over the current field.
	MaxAbs float64
	// StepCount is the number of completed steps.
	StepCount int
}

// Diagnostics computes a snapshot of the solver's diagnostic counters.
func (s *Solver) Diagnostics() Diagnostics {
	c := &s.cfg
	dx, dy := c.Dx(), c.Dy()
	dt, cc := c.Dt, c.C
	u, up := s.curr.Data, s.prev.Data
	nx, ny := c.Nx, c.Ny

	var kin, pot, maxAbs float64
	for j := 1; j < ny-1; j++ {
		base := j * nx
		for i := 1; i < nx-1; i++ {
			idx := base + i
			v := (u[idx] - up[idx]) / dt
			kin += 0.5 * v * v
			gx := (u[idx+1] - u[idx-1]) / (2 * dx)
			gy := (u[idx+nx] - u[idx-nx]) / (2 * dy)
			pot += 0.5 * cc * cc * (gx*gx + gy*gy)
			if a := math.Abs(u[idx]); a > maxAbs {
				maxAbs = a
			}
		}
	}
	return Diagnostics{
		Energy:    (kin + pot) * dx * dy,
		MaxAbs:    maxAbs,
		StepCount: s.step,
	}
}

// finite reports whether every cell of the current field is finite.
func (s *Solver) finite() bool {
	for _, v := range s.curr.Data {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return false
		}
	}
	return true
}
