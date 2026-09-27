package solver

import "sync"

// updateInterior computes
//
//	u^(n+1) = 2u^n - u^(n-1)
//	  + rx^2 (u^n[i-1,j] - 2u^n[i,j] + u^n[i+1,j])
//	  + ry^2 (u^n[i,j-1] - 2u^n[i,j] + u^n[i,j+1])
//
// for all interior points, where rx = C*dt/dx, ry = C*dt/dy.
//
// Interior rows are split into contiguous stripes, one per worker. Each
// element is written by exactly one goroutine and all reads are from
// uCurr/uPrev, so floating point results are bit-identical for any
// Worker count (including 1); only scheduling differs.
func (s *Solver) updateInterior() {
	if s.cfg.Workers <= 1 {
		s.updateInteriorRange(1, s.nx-1)
		return
	}

	workers := s.cfg.Workers
	interiorRows := s.nx - 2 // rows i = 1 .. nx-2
	if workers > interiorRows {
		workers = interiorRows
	}
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		i0 := 1 + w*interiorRows/workers
		i1 := 1 + (w+1)*interiorRows/workers
		if i1 <= i0 {
			continue
		}
		wg.Add(1)
		go func() {
			defer wg.Done()
			s.updateInteriorRange(i0, i1)
		}()
	}
	wg.Wait()
}

// updateInteriorRange updates interior rows i in [i0, i1).
func (s *Solver) updateInteriorRange(i0, i1 int) {
	ny := s.ny
	rx2 := s.rx2
	ry2 := s.ry2
	prev := s.uPrev
	curr := s.uCurr
	next := s.uNext
	for i := i0; i < i1; i++ {
		row := i * ny
		for j := 1; j < ny-1; j++ {
			idx := row + j
			lap := rx2*(curr[idx-ny]-2*curr[idx]+curr[idx+ny]) +
				ry2*(curr[idx-1]-2*curr[idx]+curr[idx+1])
			next[idx] = 2*curr[idx] - prev[idx] + lap
		}
	}
}
