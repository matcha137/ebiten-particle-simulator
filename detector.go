package main

import (
	"math"
	"math/rand"
	"sort"
)

// ============================================================================
// 検出器ジオメトリ（ATLAS を大まかに模した r-phi 断面、単位 m）
// ============================================================================

const (
	RBeamPipe = 0.034
	RTrkIn    = 0.050
	RTrkOut   = 1.100
	RSolenoid = 1.250
	REcalIn   = 1.500
	REcalOut  = 2.050
	RHcalIn   = 2.100
	RHcalOut  = 3.900
	RMuIn     = 4.300
	RMuOut    = 7.000

	NPhiTowers = 72
)

// ============================================================================
// 飛跡（ソレノイド磁場中のヘリックスの横断面射影 = 円弧）
// 曲率 1/R = 0.3 q B / pT は厳密な関係式。
// ============================================================================

type Vec2 struct{ X, Y float64 }

// bendRadius は横運動量 pt [GeV] の荷電粒子がソレノイド磁場 B [T] 中で描く
// 円弧の半径 [m]。1/R = 0.3 q B / pT は厳密な関係式。
func bendRadius(pt, q, B float64) float64 {
	return pt / (0.299792458 * math.Abs(q) * B)
}

// TrackPath は原点から出た荷電粒子の円弧を返す（半径 rStop まで、または反転点まで）。
func TrackPath(p P4, q, B, rStop float64) []Vec2 {
	pt := p.Pt()
	if pt < 1e-6 {
		return nil
	}
	vx, vy := p.Px/pt, p.Py/pt

	if q == 0 || B == 0 {
		return []Vec2{{0, 0}, {vx * rStop, vy * rStop}}
	}

	R := bendRadius(pt, q, B)
	sgn := math.Copysign(1, q*B)

	// 中心は向心力の向き: q*(v x B)
	cx, cy := R*sgn*vy, -R*sgn*vx
	// 中心から粒子へのベクトル
	ux, uy := -cx, -cy

	pts := []Vec2{{0, 0}}
	const dl = 0.03
	maxArc := 2 * math.Pi * R
	for s := dl; s <= maxArc; s += dl {
		a := -sgn * s / R
		sn, cs := math.Sincos(a)
		x := cx + ux*cs - uy*sn
		y := cy + ux*sn + uy*cs
		pts = append(pts, Vec2{x, y})
		if math.Hypot(x, y) >= rStop {
			break
		}
		if len(pts) > 4000 {
			break
		}
	}
	return pts
}

// ImpactPhi はカロリメータ面 (r = REcalIn) に到達する方位角を返す。
// 磁場で曲げられるため、生成時の phi とは異なる（低 pT ほど大きくずれる）。
func ImpactPhi(p P4, q, B float64) (float64, bool) {
	pt := p.Pt()
	if pt < 1e-6 {
		return 0, false
	}
	if q == 0 || B == 0 {
		return p.Phi(), true
	}
	R := bendRadius(pt, q, B)
	if 2*R < REcalIn {
		return 0, false // 内部トラッカー内でループして到達しない
	}
	// 弦の長さ REcalIn に対する回転角
	alpha := 2 * math.Asin(REcalIn/(2*R))
	sgn := math.Copysign(1, q*B)
	return p.Phi() - sgn*alpha/2, true
}

// ReachesTracker は荷電粒子が内部トラッカーにヒットを残すかを返す。
// 円弧の最大到達半径は 2R なので、2R > RTrkIn が幾何学的な最低条件。
// カロリメータに届かずにループで落ちる粒子でも、トラッカー内では飛跡になる
// （この判定を ImpactPhi と共有してはいけない）。
//
// 【近似】実験が課す最小ヒット数の要求は入れていないので、最内層をかすめる
// だけの飛跡も 1 本と数える。飛跡多重度は数割の過大評価になりうる。
func ReachesTracker(p P4, q, B float64) bool {
	if q == 0 {
		return false
	}
	pt := p.Pt()
	if pt < 1e-6 {
		return false
	}
	if B == 0 {
		return true
	}
	return 2*bendRadius(pt, q, B) > RTrkIn
}

// ============================================================================
// エネルギー分配と分解能
// ============================================================================

// 分解能パラメータ sigma/E = a/sqrt(E) (+) b   ((+) は二乗和)
const (
	ecalA, ecalB = 0.10, 0.007
	hcalA, hcalB = 0.50, 0.030
	// 内部トラッカー: sigma(pT)/pT = c1 (+) c2*pT
	trkC1, trkC2 = 0.010, 0.00040
	// ミューオンスペクトロメータ
	muC1, muC2 = 0.021, 0.00015
)

func smear(v, rel float64) float64 {
	x := v * (1 + rand.NormFloat64()*rel)
	if x < 0 {
		return 0
	}
	return x
}

func ecalRes(E float64) float64 { return math.Hypot(ecalA/math.Sqrt(math.Max(E, 0.1)), ecalB) }
func hcalRes(E float64) float64 { return math.Hypot(hcalA/math.Sqrt(math.Max(E, 0.1)), hcalB) }

// depositFractions は電磁 / ハドロンカロリメータへの分配。
// 【近似】シャワーの縦方向発達は解いていない。典型値を使った定数分配。
func depositFractions(k Kind, E float64) (em, had float64) {
	switch k {
	case KindEM:
		return 0.98 * E, 0.02 * E
	case KindHadron:
		return 0.30 * E, 0.65 * E // 残りは不可視エネルギー（非補償性）
	case KindMuon:
		return math.Min(0.3, E), math.Min(2.0, E) // MIP
	}
	return 0, 0
}

// ============================================================================
// 再構成
// ============================================================================

type RecoObj struct {
	Label  string // "mu", "e/gamma", "jet"
	P      P4
	Charge float64
}

type Towers struct {
	Ecal [NPhiTowers]float64 // 各 phi ビンの ET [GeV]
	Hcal [NPhiTowers]float64
}

type Reco struct {
	Objects  []RecoObj
	Jets     []P4
	Towers   Towers
	METx     float64
	METy     float64
	MET      float64
	InvMass  float64 // 主要 2 オブジェクトの不変質量
	MT       float64 // 横質量（W 事象用）
	MassOK   bool
	MTOK     bool
	NTracks  int
	SumEt    float64
	TruthMET float64
}

// phiBin は方位角をタワー番号に落とす。
// ImpactPhi は磁場による曲がりを phi に足し引きするため [-pi, pi] を
// はみ出した値を返しうる。クランプすると bin 0 / 71 に偽のホットタワーが
// 立つので、2pi で折り返す。
func phiBin(phi float64) int {
	b := int(math.Floor((phi + math.Pi) / (2 * math.Pi) * NPhiTowers))
	return ((b % NPhiTowers) + NPhiTowers) % NPhiTowers
}

// Reconstruct は真の粒子リストから検出器応答と再構成量を計算する。
func Reconstruct(ev Event, B float64) *Reco {
	r := &Reco{}

	var truthNux, truthNuy float64

	for _, p := range ev.Particles {
		if p.Kind == KindNeutrino {
			truthNux += p.P.Px
			truthNuy += p.P.Py
			continue
		}

		// --- 飛跡カウント: 内部トラッカーで完結するのでカロリメータ到達とは独立 ---
		if p.Charge != 0 && p.P.Pt() > 0.4 && ReachesTracker(p.P, p.Charge, B) {
			r.NTracks++
		}

		// --- カロリメータへの到達判定（磁場でループする粒子は落ちる） ---
		iphi, reach := ImpactPhi(p.P, p.Charge, B)
		if !reach {
			continue
		}
		et := p.P.Pt()
		em, had := depositFractions(p.Kind, p.P.E)
		emEt := em * et / math.Max(p.P.E, 1e-9)
		hadEt := had * et / math.Max(p.P.E, 1e-9)

		b := phiBin(iphi)
		r.Towers.Ecal[b] += smear(emEt, ecalRes(em))
		r.Towers.Hcal[b] += smear(hadEt, hcalRes(had))

		// --- オブジェクト再構成 ---
		switch p.Kind {
		case KindMuon:
			// ミューオンスペクトロメータ + 内部飛跡の複合測定
			rel := math.Hypot(muC1, muC2*p.P.Pt())
			s := smear(p.P.Pt(), rel) / math.Max(p.P.Pt(), 1e-9)
			r.Objects = append(r.Objects, RecoObj{"mu", p.P.Scale(s), p.Charge})
		case KindEM:
			if p.P.Pt() > 2.0 {
				s := smear(p.P.E, ecalRes(p.P.E)) / math.Max(p.P.E, 1e-9)
				r.Objects = append(r.Objects, RecoObj{"e/gamma", p.P.Scale(s), p.Charge})
			} else {
				r.Objects = append(r.Objects, RecoObj{"soft", p.P, p.Charge})
			}
		default:
			s := smear(p.P.E, hcalRes(p.P.E)) / math.Max(p.P.E, 1e-9)
			r.Objects = append(r.Objects, RecoObj{"had", p.P.Scale(s), p.Charge})
		}
	}

	// --- MET: 再構成した可視オブジェクトの横運動量ベクトル和の逆符号 ---
	for _, o := range r.Objects {
		r.METx -= o.P.Px
		r.METy -= o.P.Py
		r.SumEt += o.P.Pt()
	}
	r.MET = math.Hypot(r.METx, r.METy)
	r.TruthMET = math.Hypot(truthNux, truthNuy)

	// --- ジェット再構成（シード付きコーン、R = 0.4） ---
	r.Jets = coneJets(r.Objects, 0.4, 15.0)

	// --- 不変質量 / 横質量 ---
	switch ev.Proc {
	case ProcZmumu:
		r.InvMass, r.MassOK = leadingPairMass(r.Objects, "mu")
	case ProcZee, ProcHgg:
		r.InvMass, r.MassOK = leadingPairMass(r.Objects, "e/gamma")
	case ProcZjets:
		if len(r.Jets) >= 2 {
			r.InvMass = r.Jets[0].Add(r.Jets[1]).M()
			r.MassOK = true
		}
	case ProcWmunu:
		var lead *RecoObj
		for i := range r.Objects {
			if r.Objects[i].Label != "mu" {
				continue
			}
			if lead == nil || r.Objects[i].P.Pt() > lead.P.Pt() {
				lead = &r.Objects[i]
			}
		}
		if lead != nil {
			dphi := DeltaPhi(lead.P.Phi(), math.Atan2(r.METy, r.METx))
			r.MT = math.Sqrt(2 * lead.P.Pt() * r.MET * (1 - math.Cos(dphi)))
			r.MTOK = true
		}
	}

	return r
}

func leadingPairMass(objs []RecoObj, label string) (float64, bool) {
	var sel []P4
	for _, o := range objs {
		if o.Label == label {
			sel = append(sel, o.P)
		}
	}
	if len(sel) < 2 {
		return 0, false
	}
	sort.Slice(sel, func(i, j int) bool { return sel[i].Pt() > sel[j].Pt() })
	return sel[0].Add(sel[1]).M(), true
}

// coneJets は単純なシード付きコーン法。
// 【近似】実験で使う anti-kt ではない。赤外・共線安全性も保証されない。
func coneJets(objs []RecoObj, R, minPt float64) []P4 {
	idx := make([]int, 0, len(objs))
	for i, o := range objs {
		if o.P.Pt() > 0.5 {
			idx = append(idx, i)
		}
	}
	sort.Slice(idx, func(a, b int) bool { return objs[idx[a]].P.Pt() > objs[idx[b]].P.Pt() })

	used := make(map[int]bool)
	var jets []P4
	for _, i := range idx {
		if used[i] {
			continue
		}
		axis := objs[i].P
		var jet P4
		var members []int
		for iter := 0; iter < 3; iter++ {
			jet = P4{}
			members = members[:0]
			for _, j := range idx {
				if used[j] {
					continue
				}
				if DeltaR(objs[j].P, axis) < R {
					jet = jet.Add(objs[j].P)
					members = append(members, j)
				}
			}
			if jet.Pt() < 1e-6 {
				break
			}
			axis = jet
		}
		if jet.Pt() < minPt {
			continue
		}
		for _, j := range members {
			used[j] = true
		}
		jets = append(jets, jet)
	}
	sort.Slice(jets, func(a, b int) bool { return jets[a].Pt() > jets[b].Pt() })
	return jets
}
