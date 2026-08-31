package main

import (
	"math"
	"math/rand"
)

// ============================================================================
// 終状態粒子
// ============================================================================

type Kind int

const (
	KindEM       Kind = iota // e, gamma  -> 電磁カロリメータで止まる
	KindHadron               // pi, K, p, n -> ハドロンカロリメータ
	KindMuon                 // ミューオン -> 全層を貫通
	KindNeutrino             // 検出されない -> MET になる
)

type Particle struct {
	Name   string
	Charge float64
	Kind   Kind
	Mass   float64
	P      P4
}

// 粒子質量 [GeV]
const (
	mElectron = 0.000511
	mMuon     = 0.105658
	mPion     = 0.139570
	mKaon     = 0.493677
	mProton   = 0.938272

	mZ, gZ = 91.1876, 2.4952
	mW, gW = 80.377, 2.085
	mH, gH = 125.25, 0.0041
)

// ============================================================================
// 生成される事象
// ============================================================================

type Process int

const (
	ProcZmumu Process = iota
	ProcZee
	ProcZjets
	ProcHgg
	ProcWmunu
	ProcMinBias
)

var processNames = []string{
	"Z -> mu+ mu-",
	"Z -> e+ e-",
	"Z -> q qbar (2 jets)",
	"H -> gamma gamma",
	"W+ -> mu+ nu",
	"minimum bias (soft QCD)",
}

type Event struct {
	Proc      Process
	Particles []Particle
	TruthMass float64 // 生成した共鳴の真の質量
}

// bosonP4 は pp 衝突で生成されるボゾンの 4元運動量。
//
// 【近似】pT スペクトルとラピディティ分布は PDF から計算した本物ではなく、
// 実測に大まかに合わせた経験的な形（pT ~ Gamma(2, 6GeV)、y ~ Gauss(0,1.6)）。
// ここだけは PDF + 行列要素の計算をしていないので、厳密ではない。
func bosonP4(m float64) P4 {
	pt := -6.0 * (math.Log(1-rand.Float64()) + math.Log(1-rand.Float64()))
	phi := 2 * math.Pi * rand.Float64()
	y := rand.NormFloat64() * 1.6

	mt := math.Hypot(m, pt) // 横質量
	return P4{
		Px: pt * math.Cos(phi),
		Py: pt * math.Sin(phi),
		Pz: mt * math.Sinh(y),
		E:  mt * math.Cosh(y),
	}
}

// softRecoil はボゾンの横運動量を打ち消すソフトなハドロン群を作る。
// 横方向の運動量保存は厳密に成立させる（sum pT = -pT(boson)）。
// 縦方向はビーム残骸に逃げるので保存させない（実際の pp 衝突と同じ扱い）。
func softRecoil(boson P4, n int) []Particle {
	type v2 struct{ x, y float64 }
	raw := make([]v2, n)
	var sx, sy float64
	for i := range raw {
		pt := ExpSample(0.9)
		ph := 2 * math.Pi * rand.Float64()
		raw[i] = v2{pt * math.Cos(ph), pt * math.Sin(ph)}
		sx += raw[i].x
		sy += raw[i].y
	}
	// 残差を等分して厳密に釣り合わせる
	cx := (-boson.Px - sx) / float64(n)
	cy := (-boson.Py - sy) / float64(n)

	out := make([]Particle, 0, n)
	for i := range raw {
		px := raw[i].x + cx
		py := raw[i].y + cy
		eta := rand.NormFloat64() * 2.0
		pt := math.Hypot(px, py)
		pz := pt * math.Sinh(eta)
		m := mPion
		e := math.Sqrt(px*px + py*py + pz*pz + m*m)

		q, name := 0.0, "n"
		switch rand.Intn(3) {
		case 0:
			q, name = +1, "pi+"
		case 1:
			q, name = -1, "pi-"
		}
		out = append(out, Particle{name, q, KindHadron, m, P4{px, py, pz, e}})
	}
	return out
}

// fragment はパートンをハドロン群に変換する。
//
// 【近似】これは本物のフラグメンテーションではない。Lund 弦模型を真似た
// おもちゃモデル（縦方向の分配 z ~ (1-z)^2、横方向のキック sigma = 0.35 GeV）で、
// 最後に全体をスケールしてパートンの運動量に合わせている。
// ジェットの本数・質量分布は定量的には信用できない。
func fragment(parton P4) []Particle {
	dir := 1.0 / parton.P()
	ux, uy, uz := parton.Px*dir, parton.Py*dir, parton.Pz*dir

	// パートン方向に直交する 2 軸
	ax, ay, az := 0.0, 0.0, 1.0
	if math.Abs(uz) > 0.9 {
		ax, ay, az = 1.0, 0.0, 0.0
	}
	e1x, e1y, e1z := uy*az-uz*ay, uz*ax-ux*az, ux*ay-uy*ax
	n1 := math.Sqrt(e1x*e1x + e1y*e1y + e1z*e1z)
	e1x, e1y, e1z = e1x/n1, e1y/n1, e1z/n1
	e2x, e2y, e2z := uy*e1z-uz*e1y, uz*e1x-ux*e1z, ux*e1y-uy*e1x

	remaining := parton.P()
	var hadrons []Particle
	for remaining > 1.0 && len(hadrons) < 30 {
		u := rand.Float64()
		z := 1 - math.Cbrt(1-u) // (1-z)^2 分布
		pl := z * remaining
		if pl < 0.3 {
			pl = 0.3
		}
		remaining -= pl

		kt := math.Abs(rand.NormFloat64()) * 0.35
		kphi := 2 * math.Pi * rand.Float64()
		k1, k2 := kt*math.Cos(kphi), kt*math.Sin(kphi)

		px := pl*ux + k1*e1x + k2*e2x
		py := pl*uy + k1*e1y + k2*e2y
		pz := pl*uz + k1*e1z + k2*e2z

		q, name, m := 0.0, "pi0", mPion
		switch rand.Intn(4) {
		case 0:
			q, name = +1, "pi+"
		case 1:
			q, name = -1, "pi-"
		case 2:
			q, name, m = +1, "K+", mKaon
		}
		e := math.Sqrt(px*px + py*py + pz*pz + m*m)
		kind := KindHadron
		if name == "pi0" {
			kind = KindEM // pi0 -> gamma gamma として電磁シャワー扱い
		}
		hadrons = append(hadrons, Particle{name, q, kind, m, P4{px, py, pz, e}})
	}

	// 生成したハドロンの運動量をパートンに合わせて再スケール
	var sum P4
	for _, h := range hadrons {
		sum = sum.Add(h.P)
	}
	if sum.P() > 0 {
		s := parton.P() / sum.P()
		for i := range hadrons {
			p := hadrons[i].P.Scale(s)
			m := hadrons[i].Mass
			p.E = math.Sqrt(p.Px*p.Px + p.Py*p.Py + p.Pz*p.Pz + m*m)
			hadrons[i].P = p
		}
	}
	return hadrons
}

// Generate は 1 事象を生成する。
func Generate(proc Process) Event {
	ev := Event{Proc: proc}

	switch proc {
	case ProcZmumu, ProcZee:
		m := BreitWigner(mZ, gZ, 60, 120)
		boson := bosonP4(m)
		ev.TruthMass = m
		dm, name, kind := mMuon, "mu", KindMuon
		if proc == ProcZee {
			dm, name, kind = mElectron, "e", KindEM
		}
		a, b := TwoBodyDecay(boson, dm, dm)
		ev.Particles = append(ev.Particles,
			Particle{name + "+", +1, kind, dm, a},
			Particle{name + "-", -1, kind, dm, b},
		)
		ev.Particles = append(ev.Particles, softRecoil(boson, 8+rand.Intn(8))...)

	case ProcZjets:
		m := BreitWigner(mZ, gZ, 60, 120)
		boson := bosonP4(m)
		ev.TruthMass = m
		// Z -> q qbar。質量はコンスティチュエント質量で近似。
		q1, q2 := TwoBodyDecay(boson, 0.33, 0.33)
		ev.Particles = append(ev.Particles, fragment(q1)...)
		ev.Particles = append(ev.Particles, fragment(q2)...)
		ev.Particles = append(ev.Particles, softRecoil(boson, 10+rand.Intn(10))...)

	case ProcHgg:
		m := BreitWigner(mH, gH, 120, 130)
		boson := bosonP4(m)
		ev.TruthMass = m
		a, b := TwoBodyDecay(boson, 0, 0)
		ev.Particles = append(ev.Particles,
			Particle{"gamma", 0, KindEM, 0, a},
			Particle{"gamma", 0, KindEM, 0, b},
		)
		ev.Particles = append(ev.Particles, softRecoil(boson, 10+rand.Intn(10))...)

	case ProcWmunu:
		m := BreitWigner(mW, gW, 50, 110)
		boson := bosonP4(m)
		ev.TruthMass = m
		mu, nu := TwoBodyDecay(boson, mMuon, 0)
		ev.Particles = append(ev.Particles,
			Particle{"mu+", +1, KindMuon, mMuon, mu},
			Particle{"nu", 0, KindNeutrino, 0, nu},
		)
		ev.Particles = append(ev.Particles, softRecoil(boson, 8+rand.Intn(8))...)

	case ProcMinBias:
		ev.Particles = softRecoil(P4{}, 20+rand.Intn(20))
	}

	return ev
}
