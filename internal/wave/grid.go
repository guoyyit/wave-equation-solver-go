package wave

// Grid stores a scalar field on a regular Nx by Ny lattice, row-major.
type Grid struct {
	Nx, Ny int
	Data   []float64
}

// NewGrid allocates a zeroed Nx by Ny grid.
func NewGrid(nx, ny int) *Grid { return &Grid{Nx: nx, Ny: ny, Data: make([]float64, nx*ny)} }

// At returns the value at (i, j) without bounds checks.
func (g *Grid) At(i, j int) float64 { return g.Data[j*g.Nx+i] }

// Set writes the value at (i, j) without bounds checks.
func (g *Grid) Set(i, j int, v float64) { g.Data[j*g.Nx+i] = v }

// Clone deep-copies the grid.
func (g *Grid) Clone() *Grid {
	out := &Grid{Nx: g.Nx, Ny: g.Ny, Data: make([]float64, len(g.Data))}
	copy(out.Data, g.Data)
	return out
}
