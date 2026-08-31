package main

import (
	"math"
	"math/rand"
)

// ============================================================================
// 4元運動量。単位は自然単位系 (c = 1), GeV。
// この節の計算はすべて厳密（近似なし）。
// ============================================================================

type P4 struct{ Px, Py, Pz, E float64 }

func (a P4) Add(b P4) P4 { return P4{a.Px + b.Px, a.Py + b.Py, a.Pz + b.Pz, a.E + b.E} }

func (a P4) Scale(s float64) P4 { return P4{a.Px * s, a.Py * s, a.Pz * s, a.E * s} }

func (a P4) P() float64 { return math.Sqrt(a.Px*a.Px + a.Py*a.Py + a.Pz*a.Pz) }

func (a P4) Pt() float64 { return math.Hypot(a.Px, a.Py) }

func (a P4) Phi() float64 { return math.Atan2(a.Py, a.Px) }

// M は不変質量 sqrt(E^2 - |p|^2)。数値誤差で負になった場合は 0 を返す。
func (a P4) M() float64 {
	m2 := a.E*a.E - a.Px*a.Px - a.Py*a.Py - a.Pz*a.Pz
	if m2 <= 0 {
		return 0
	}
	return math.Sqrt(m2)
}

// Eta は擬ラピディティ -ln(tan(theta/2))。
func (a P4) Eta() float64 {
	p := a.P()
	if p-math.Abs(a.Pz) < 1e-12 {
		if a.Pz >= 0 {
			return 10
		}
		return -10
	}
	return 0.5 * math.Log((p+a.Pz)/(p-a.Pz))
}

// Rapidity は真のラピディティ 0.5*ln((E+pz)/(E-pz))。
func (a P4) Rapidity() float64 {
	if a.E-math.Abs(a.Pz) < 1e-12 {
		return math.Copysign(10, a.Pz)
	}
	return 0.5 * math.Log((a.E+a.Pz)/(a.E-a.Pz))
}

// BoostVector は beta = p/E。
func (a P4) BoostVector() (float64, float64, float64) {
	return a.Px / a.E, a.Py / a.E, a.Pz / a.E
}

// Boost は速度 beta = (bx,by,bz) のローレンツブースト（厳密）。
func (a P4) Boost(bx, by, bz float64) P4 {
	b2 := bx*bx + by*by + bz*bz
	if b2 <= 0 || b2 >= 1 {
		return a
	}
	gamma := 1 / math.Sqrt(1-b2)
	bp := bx*a.Px + by*a.Py + bz*a.Pz
	g2 := (gamma - 1) / b2
	return P4{
		Px: a.Px + g2*bp*bx + gamma*bx*a.E,
		Py: a.Py + g2*bp*by + gamma*by*a.E,
		Pz: a.Pz + g2*bp*bz + gamma*bz*a.E,
		E:  gamma * (a.E + bp),
	}
}

// DeltaPhi は [-pi, pi] に折り返した方位角差。
func DeltaPhi(p1, p2 float64) float64 {
	d := math.Mod(p1-p2+3*math.Pi, 2*math.Pi) - math.Pi
	return d
}

// DeltaR は eta-phi 平面上の距離。
func DeltaR(a, b P4) float64 {
	de := a.Eta() - b.Eta()
	dp := DeltaPhi(a.Phi(), b.Phi())
	return math.Hypot(de, dp)
}

// ============================================================================
// 崩壊
// ============================================================================

// TwoBodyDecay は親の静止系で等方的に 2 体崩壊させ、実験室系へブーストする。
// 4元運動量は厳密に保存する（浮動小数点精度の範囲で）。
func TwoBodyDecay(parent P4, m1, m2 float64) (P4, P4) {
	M := parent.M()
	if M <= m1+m2 {
		// 運動学的に禁止。親を質量比で分けて返す（保存のみ守る）。
		f := m1 / (m1 + m2)
		return parent.Scale(f), parent.Scale(1 - f)
	}

	// 静止系での娘粒子の運動量（Källén 関数）
	pstar := math.Sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2))) / (2 * M)

	cosT := 2*rand.Float64() - 1
	sinT := math.Sqrt(1 - cosT*cosT)
	phi := 2 * math.Pi * rand.Float64()

	dx := pstar * sinT * math.Cos(phi)
	dy := pstar * sinT * math.Sin(phi)
	dz := pstar * cosT

	p1 := P4{dx, dy, dz, math.Hypot(pstar, m1)}
	p2 := P4{-dx, -dy, -dz, math.Hypot(pstar, m2)}

	bx, by, bz := parent.BoostVector()
	return p1.Boost(bx, by, bz), p2.Boost(bx, by, bz)
}

// ============================================================================
// 共鳴の質量分布
// ============================================================================

// BreitWigner は相対論的 Breit-Wigner から質量をサンプルする。
// s = M^2 + M*Gamma*tan(theta) の変数変換で逆関数法を使う（区間で切断）。
func BreitWigner(M, Gamma, lo, hi float64) float64 {
	tLo := math.Atan((lo*lo - M*M) / (M * Gamma))
	tHi := math.Atan((hi*hi - M*M) / (M * Gamma))
	t := tLo + rand.Float64()*(tHi-tLo)
	s := M*M + M*Gamma*math.Tan(t)
	if s <= 0 {
		return M
	}
	return math.Sqrt(s)
}

// ExpSample は平均 mean の指数分布。
func ExpSample(mean float64) float64 {
	return -mean * math.Log(1-rand.Float64())
}
