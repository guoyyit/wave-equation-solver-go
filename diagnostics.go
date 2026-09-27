package solver

import "math"

// Diagnostics summarizes the current solution and run state.
type Diagnostics struct {
	Step       int     `json:"step"`
	Time       float64 `json:"time"`
	MaxAbsU    float64 `json:"max_abs_u"`
	Energy     float64 `json:"energy"`
	Courant    float64 `json:"courant"`
	Workers    int     `json:"workers"`
	Checkpoint string  `json:"checkpoint"`
}

// ComputeDiagnostics evaluates max-amplitude and discrete-energy metrics
// for the current field.
//
// Energy is a discrete proxy for interior kinetic + potential energy,
// averaging over cells and using one-sided differences where needed:
//
//	E = (1/N) * sum [ 0.5*C^2*( (du/dx)^2 + (du/dy)^2 ) ]
//
// It is used comparatively (absorbing vs. fixed reference), not as a
// strictly conserved Hamiltonian for the bounded domain.
func (s *Solver) ComputeDiagnostics() Diagnostics {
	nx, ny := s.nx, s.ny
	u := s.uCurr
	c2 := s.cfg.C * s.cfg.C
	invDx := 1.0 / s.dx
	invDy := 1.0 / s.dy

	var maxAbs, energy float64
	for i := 0; i < nx; i++ {
		row := i * ny
		for j := 0; j < ny; j++ {
			idx := row + j
			v := u[idx]
			if a := math.Abs(v); a > maxAbs {
				maxAbs = a
			}

			var dudx, dudy float64
			switch {
			case i == 0:
				dudx = (u[idx+ny] - v) * invDx
			case i == nx-1:
				dudx = (v - u[idx-ny]) * invDx
			default:
				dudx = (u[idx+ny] - u[idx-ny]) * 0.5 * invDx
			}
			switch {
			case j == 0:
				dudy = (u[idx+1] - v) * invDy
			case j == ny-1:
				dudy = (v - u[idx-1]) * invDy
			default:
				dudy = (u[idx+1] - u[idx-1]) * 0.5 * invDy
			}
			energy += 0.5 * c2 * (dudx*dudx + dudy*dudy)
		}
	}
	energy /= float64(nx * ny)

	return Diagnostics{
		Step:       s.step,
		Time:       s.time,
		MaxAbsU:    maxAbs,
		Energy:     energy,
		Courant:    s.cfg.courant(s.cfg.Dt),
		Workers:    s.cfg.Workers,
		Checkpoint: s.cfg.CheckpointDir,
	}
}
