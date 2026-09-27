package wave

// zeroFixedEdges sets u = 0 on every edge configured BoundaryFixed.
// Absorbing edges are left untouched: they are enforced by a time-update
// rule (Mur) during stepping, not by constraining the field value.
func (s *Solver) zeroFixedEdges(g *Grid) {
	c := &s.cfg
	nx, ny := c.Nx, c.Ny
	if c.Boundaries.Left == BoundaryFixed {
		for j := 0; j < ny; j++ {
			g.Data[j*nx] = 0
		}
	}
	if c.Boundaries.Right == BoundaryFixed {
		for j := 0; j < ny; j++ {
			g.Data[j*nx+nx-1] = 0
		}
	}
	if c.Boundaries.Bottom == BoundaryFixed {
		for i := 0; i < nx; i++ {
			g.Data[i] = 0
		}
	}
	if c.Boundaries.Top == BoundaryFixed {
		for i := 0; i < nx; i++ {
			g.Data[(ny-1)*nx+i] = 0
		}
	}
}

// applyBoundaries enforces the configured boundary conditions on the next
// time level after the interior update. It runs before the level rotation
// in Step, so s.next is u^{n+1} and s.curr is u^n.
//
// Absorbing edges use the first-order Mur condition, e.g. on the left edge
//
//	u^{n+1}_{0,j} = u^n_{1,j} + (r-1)/(r+1) * (u^{n+1}_{1,j} - u^n_{0,j})
//
// with r = C*Dt/Dx. The coefficient magnitude is < 1 whenever the CFL
// condition holds, so the boundary update cannot introduce instability.
//
// Corners are handled explicitly so mixed edge types stay consistent:
// a corner touching any fixed edge is held at 0; a corner where both
// edges absorb takes the average of the two one-dimensional Mur updates.
func (s *Solver) applyBoundaries() {
	c := &s.cfg
	nx, ny := c.Nx, c.Ny
	next, curr := s.next.Data, s.curr.Data
	kx := (c.C*c.Dt/c.Dx() - 1) / (c.C*c.Dt/c.Dx() + 1)
	ky := (c.C*c.Dt/c.Dy() - 1) / (c.C*c.Dt/c.Dy() + 1)

	// Edges, corners excluded.
	if c.Boundaries.Left == BoundaryFixed {
		for j := 1; j < ny-1; j++ {
			next[j*nx] = 0
		}
	} else {
		for j := 1; j < ny-1; j++ {
			next[j*nx] = curr[j*nx+1] + kx*(next[j*nx+1]-curr[j*nx])
		}
	}
	if c.Boundaries.Right == BoundaryFixed {
		for j := 1; j < ny-1; j++ {
			next[j*nx+nx-1] = 0
		}
	} else {
		for j := 1; j < ny-1; j++ {
			next[j*nx+nx-1] = curr[j*nx+nx-2] + kx*(next[j*nx+nx-2]-curr[j*nx+nx-1])
		}
	}
	if c.Boundaries.Bottom == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			next[i] = 0
		}
	} else {
		for i := 1; i < nx-1; i++ {
			next[i] = curr[nx+i] + ky*(next[nx+i]-curr[i])
		}
	}
	if c.Boundaries.Top == BoundaryFixed {
		for i := 1; i < nx-1; i++ {
			next[(ny-1)*nx+i] = 0
		}
	} else {
		for i := 1; i < nx-1; i++ {
			next[(ny-1)*nx+i] = curr[(ny-2)*nx+i] + ky*(next[(ny-2)*nx+i]-curr[(ny-1)*nx+i])
		}
	}

	// Corners.
	s.corner(0, 0, c.Boundaries.Left, c.Boundaries.Bottom,
		curr[1]+kx*(next[1]-curr[0]),       // from the left edge
		curr[nx]+ky*(next[nx]-curr[0]))     // from the bottom edge
	s.corner(nx-1, 0, c.Boundaries.Right, c.Boundaries.Bottom,
		curr[nx-2]+kx*(next[nx-2]-curr[nx-1]),           // from the right edge
		curr[2*nx-1]+ky*(next[2*nx-1]-curr[nx-1]))       // from the bottom edge
	s.corner(0, ny-1, c.Boundaries.Left, c.Boundaries.Top,
		curr[(ny-1)*nx+1]+kx*(next[(ny-1)*nx+1]-curr[(ny-1)*nx]), // from the left edge
		curr[(ny-2)*nx]+ky*(next[(ny-2)*nx]-curr[(ny-1)*nx]))       // from the top edge
	s.corner(nx-1, ny-1, c.Boundaries.Right, c.Boundaries.Top,
		curr[(ny-1)*nx+nx-2]+kx*(next[(ny-1)*nx+nx-2]-curr[(ny-1)*nx+nx-1]), // right
		curr[(ny-2)*nx+nx-1]+ky*(next[(ny-2)*nx+nx-1]-curr[(ny-1)*nx+nx-1])) // top
}

// corner writes the boundary value at corner (i, j) given the types of the
// two adjacent edges and their respective one-dimensional updates a and b.
func (s *Solver) corner(i, j int, e1, e2 BoundaryType, a, b float64) {
	switch {
	case e1 == BoundaryFixed || e2 == BoundaryFixed:
		s.next.Data[j*s.cfg.Nx+i] = 0
	default:
		s.next.Data[j*s.cfg.Nx+i] = 0.5 * (a + b)
	}
}
