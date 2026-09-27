package solver

// applyBoundaries fills the boundary values of uNext after the interior
// update, according to the configured boundary kinds.
//
//   - Fixed edges: u = 0 (homogeneous Dirichlet).
//
//   - Absorbing edges: first-order Mur condition. For an edge normal to
//     x with adjacent interior value v and edge value e,
//
//     e^{n+1} = v^n + (r-1)/(r+1) * (v^{n+1} - e^n),  r = C*dt/dx
//
//     (analogously with r = C*dt/dy for edges normal to y). The update
//     only references u^n and already-computed u^(n+1) interior values,
//     so no extra history storage is required.
//
// The four corner cells are shared by two edges. They are pinned to zero
// before edge updates so that any combination of fixed/absorbing edges
// yields a single, self-consistent corner value (the two edges never
// race to write different numbers).
func (s *Solver) applyBoundaries() {
	nx, ny := s.nx, s.ny
	next, curr := s.uNext, s.uCurr
	ax := (s.rx - 1) / (s.rx + 1)
	ay := (s.ry - 1) / (s.ry + 1)

	// Corners meet two edges and are always pinned to zero. This single
	// rule makes every fixed/absorbing combination self-consistent: the
	// edge loops below skip corner columns/rows.
	next[0] = 0
	next[ny-1] = 0
	next[(nx-1)*ny] = 0
	next[(nx-1)*ny+ny-1] = 0

	// West edge (i = 0); interior neighbor i = 1.
	if s.cfg.Boundary.West == BoundaryFixed {
		for j := 1; j < ny-1; j++ {
			next[j] = 0
		}
	} else {
		for j := 1; j < ny-1; j++ {
			next[j] = curr[ny+j] + ax*(next[ny+j]-curr[j])
		}
	}

	// East edge (i = nx-1); interior neighbor i = nx-2.
	ib := (nx - 1) * ny
	ii := (nx - 2) * ny
	if s.cfg.Boundary.East == BoundaryFixed {
		for j := 1; j < ny-1; j++ {
			next[ib+j] = 0
		}
	} else {
		for j := 1; j < ny-1; j++ {
			next[ib+j] = curr[ii+j] + ax*(next[ii+j]-curr[ib+j])
		}
	}

	// South edge (j = 0); interior neighbor j = 1.
	if s.cfg.Boundary.South == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			next[i*ny] = 0
		}
	} else {
		for i := 1; i < nx-1; i++ {
			row := i * ny
			next[row] = curr[row+1] + ay*(next[row+1]-curr[row])
		}
	}

	// North edge (j = ny-1); interior neighbor j = ny-2.
	if s.cfg.Boundary.North == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			next[i*ny+ny-1] = 0
		}
	} else {
		for i := 1; i < nx-1; i++ {
			row := i * ny
			next[row+ny-1] = curr[row+ny-2] + ay*(next[row+ny-2]-curr[row+ny-1])
		}
	}
}
