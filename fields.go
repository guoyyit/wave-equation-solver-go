package solver

import "math"

// GaussianPulseXY returns a FieldFunc sampled at grid indices, where
// (cx, cy) is the pulse center in physical coordinates, dx/dy are the
// grid spacings and s is the Gaussian half-width.
func GaussianPulseXY(dx, dy, cx, cy, s, amp float64) FieldFunc {
	inv2s2 := 1.0 / (2.0 * s * s)
	return func(i, j int) float64 {
		x := float64(i) * dx
		y := float64(j) * dy
		dx2 := x - cx
		dy2 := y - cy
		return amp * math.Exp(-(dx2*dx2+dy2*dy2)*inv2s2)
	}
}

// GaussianPulse returns a pulse centered at grid index (ci, cj) with
// half-width s measured in grid cells.
func GaussianPulse(ci, cj int, s, amp float64) FieldFunc {
	inv2s2 := 1.0 / (2.0 * s * s)
	return func(i, j int) float64 {
		di := float64(i - ci)
		dj := float64(j - cj)
		return amp * math.Exp(-(di*di+dj*dj)*inv2s2)
	}
}

// ZeroField is the explicitly chosen zero displacement/velocity field.
func ZeroField() FieldFunc {
	return func(i, j int) float64 { return 0 }
}

// RandomField builds a reproducible pseudo-random field from seed, with
// values in [-amp, amp]. The value at each grid index depends only on
// (seed, i, j), so it is independent of evaluation order and of the
// number of worker goroutines.
func RandomField(seed int64, amp float64) FieldFunc {
	return func(i, j int) float64 {
		// SplitMix64-style mixing of (seed, i, j).
		z := uint64(seed) ^ mix64index(uint64(uint32(i)), uint64(uint32(j)))
		z += 0x9e3779b97f4a7c15
		z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
		z = (z ^ (z >> 27)) * 0x94d049bb133111eb
		z ^= z >> 31
		// Map [0,1) to [-amp, amp).
		f := float64(z>>11) / float64(uint64(1)<<53)
		return amp * (2*f - 1)
	}
}

func mix64index(i, j uint64) uint64 {
	x := i*0x9e3779b97f4a7c15 ^ j*0xc2b2ae3d27d4eb4f
	x ^= x >> 33
	x *= 0xff51afd7ed558ccd
	x ^= x >> 33
	return x
}
