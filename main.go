package main

import (
	"fmt"
	"image/color"
	"math"
	"strings"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
	"github.com/hajimehoshi/ebiten/v2/inpututil"
	"github.com/hajimehoshi/ebiten/v2/vector"
)

const (
	screenW = 1280
	screenH = 780

	viewCX = 385.0
	viewCY = 390.0
	pxPerM = 51.0 // 7.0 m -> 357 px
)

// ============================================================================
// ヒストグラム
// ============================================================================

type Hist struct {
	Lo, Hi float64
	Bins   []int
	N      int
	Sum    float64
	Sum2   float64
	Title  string
}

func NewHist(lo, hi float64, n int, title string) *Hist {
	return &Hist{Lo: lo, Hi: hi, Bins: make([]int, n), Title: title}
}

func (h *Hist) Fill(v float64) {
	h.N++
	h.Sum += v
	h.Sum2 += v * v
	if v < h.Lo || v >= h.Hi {
		return
	}
	i := int((v - h.Lo) / (h.Hi - h.Lo) * float64(len(h.Bins)))
	h.Bins[i]++
}

func (h *Hist) Mean() float64 {
	if h.N == 0 {
		return 0
	}
	return h.Sum / float64(h.N)
}

func (h *Hist) RMS() float64 {
	if h.N < 2 {
		return 0
	}
	m := h.Mean()
	v := h.Sum2/float64(h.N) - m*m
	if v < 0 {
		return 0
	}
	return math.Sqrt(v)
}

func (h *Hist) Max() int {
	m := 1
	for _, b := range h.Bins {
		if b > m {
			m = b
		}
	}
	return m
}

func histRangeFor(p Process) (float64, float64, string) {
	switch p {
	case ProcHgg:
		return 100, 150, "m(gamma gamma) [GeV]"
	case ProcWmunu:
		return 0, 120, "mT(mu, MET) [GeV]"
	case ProcZjets:
		return 20, 160, "m(jet1, jet2) [GeV]"
	case ProcMinBias:
		return 0, 60, "MET [GeV]"
	default:
		return 60, 120, "m(l+ l-) [GeV]"
	}
}

// ============================================================================
// Game
// ============================================================================

type Game struct {
	proc  Process
	B     float64
	ev    Event
	reco  *Reco
	hist  *Hist
	paths [][]Vec2
	cols  []color.RGBA
	kinds []Kind
}

func NewGame() *Game {
	g := &Game{proc: ProcZmumu, B: 2.0}
	g.newHist()
	g.newEvent()
	return g
}

func (g *Game) newHist() {
	lo, hi, t := histRangeFor(g.proc)
	g.hist = NewHist(lo, hi, 50, t)
}

func (g *Game) fillValue() (float64, bool) {
	if g.proc == ProcWmunu {
		return g.reco.MT, g.reco.MTOK
	}
	if g.proc == ProcMinBias {
		return g.reco.MET, true
	}
	return g.reco.InvMass, g.reco.MassOK
}

func (g *Game) newEvent() {
	g.ev = Generate(g.proc)
	g.reco = Reconstruct(g.ev, g.B)
	if v, ok := g.fillValue(); ok {
		g.hist.Fill(v)
	}
	g.buildPaths()
}

// buildPaths は表示用の飛跡を作る。
func (g *Game) buildPaths() {
	g.paths = g.paths[:0]
	g.cols = g.cols[:0]
	g.kinds = g.kinds[:0]
	for _, p := range g.ev.Particles {
		var rStop float64
		var col color.RGBA
		switch p.Kind {
		case KindMuon:
			rStop, col = RMuOut, color.RGBA{255, 60, 60, 255}
		case KindEM:
			rStop, col = REcalIn+0.25, color.RGBA{90, 210, 255, 255}
		case KindHadron:
			rStop, col = RHcalIn+0.6, color.RGBA{110, 230, 150, 255}
		case KindNeutrino:
			continue // 検出されない
		}
		if p.Charge == 0 && p.Kind != KindMuon {
			col = color.RGBA{col.R / 2, col.G / 2, col.B / 2, 200}
		}
		pts := TrackPath(p.P, p.Charge, g.B, rStop)
		if len(pts) < 2 {
			continue
		}
		g.paths = append(g.paths, pts)
		g.cols = append(g.cols, col)
		g.kinds = append(g.kinds, p.Kind)
	}
}

// fastFill は表示なしで大量に事象を生成しヒストグラムに詰める。
func (g *Game) fastFill(n int) {
	for i := 0; i < n; i++ {
		ev := Generate(g.proc)
		r := Reconstruct(ev, g.B)
		var v float64
		var ok bool
		switch g.proc {
		case ProcWmunu:
			v, ok = r.MT, r.MTOK
		case ProcMinBias:
			v, ok = r.MET, true
		default:
			v, ok = r.InvMass, r.MassOK
		}
		if ok {
			g.hist.Fill(v)
		}
	}
}

func (g *Game) Update() error {
	if inpututil.IsKeyJustPressed(ebiten.KeySpace) || inpututil.IsKeyJustPressed(ebiten.KeyR) {
		g.newEvent()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyF) {
		g.fastFill(500)
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyC) {
		g.newHist()
	}
	for i, k := range []ebiten.Key{
		ebiten.KeyDigit1, ebiten.KeyDigit2, ebiten.KeyDigit3,
		ebiten.KeyDigit4, ebiten.KeyDigit5, ebiten.KeyDigit6,
	} {
		if inpututil.IsKeyJustPressed(k) {
			g.proc = Process(i)
			g.newHist()
			g.newEvent()
		}
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyUp) {
		g.B += 0.5
		g.newEvent()
	}
	if inpututil.IsKeyJustPressed(ebiten.KeyDown) {
		if g.B > 0.5 {
			g.B -= 0.5
		}
		g.newEvent()
	}
	return nil
}

func toScreen(x, y float64) (float32, float32) {
	return float32(viewCX + x*pxPerM), float32(viewCY - y*pxPerM)
}

func (g *Game) Draw(dst *ebiten.Image) {
	dst.Fill(color.RGBA{8, 10, 16, 255})
	g.drawDetector(dst)
	g.drawTowers(dst)
	g.drawTracks(dst)
	g.drawMET(dst)
	g.drawPanel(dst)
}

func (g *Game) drawDetector(dst *ebiten.Image) {
	ring := func(r float64, c color.RGBA, w float32) {
		vector.StrokeCircle(dst, float32(viewCX), float32(viewCY), float32(r*pxPerM), w, c, true)
	}
	ring(RBeamPipe, color.RGBA{120, 120, 140, 255}, 2)
	ring(RTrkIn, color.RGBA{60, 70, 95, 255}, 1)
	ring(RTrkOut, color.RGBA{60, 70, 95, 255}, 1)
	ring(RSolenoid, color.RGBA{110, 90, 160, 255}, 2)
	ring(REcalIn, color.RGBA{80, 90, 40, 255}, 1)
	ring(REcalOut, color.RGBA{80, 90, 40, 255}, 1)
	ring(RHcalIn, color.RGBA{95, 55, 40, 255}, 1)
	ring(RHcalOut, color.RGBA{95, 55, 40, 255}, 1)
	ring(RMuIn, color.RGBA{70, 45, 55, 255}, 1)
	ring(RMuOut, color.RGBA{70, 45, 55, 255}, 1)
}

func (g *Game) drawTowers(dst *ebiten.Image) {
	maxEt := 1.0
	for i := 0; i < NPhiTowers; i++ {
		if e := g.reco.Towers.Ecal[i] + g.reco.Towers.Hcal[i]; e > maxEt {
			maxEt = e
		}
	}
	wPx := float32(2*math.Pi*REcalIn*pxPerM/NPhiTowers) * 0.85

	for i := 0; i < NPhiTowers; i++ {
		phi := -math.Pi + (float64(i)+0.5)*2*math.Pi/NPhiTowers
		cx, sy := math.Cos(phi), math.Sin(phi)

		if e := g.reco.Towers.Ecal[i]; e > 0.05 {
			l := (REcalOut - REcalIn) * math.Min(1, e/maxEt)
			x0, y0 := toScreen(REcalIn*cx, REcalIn*sy)
			x1, y1 := toScreen((REcalIn+l)*cx, (REcalIn+l)*sy)
			vector.StrokeLine(dst, x0, y0, x1, y1, wPx, color.RGBA{230, 230, 90, 255}, true)
		}
		if e := g.reco.Towers.Hcal[i]; e > 0.05 {
			l := (RHcalOut - RHcalIn) * math.Min(1, e/maxEt)
			x0, y0 := toScreen(RHcalIn*cx, RHcalIn*sy)
			x1, y1 := toScreen((RHcalIn+l)*cx, (RHcalIn+l)*sy)
			vector.StrokeLine(dst, x0, y0, x1, y1, wPx*1.3, color.RGBA{220, 110, 70, 255}, true)
		}
	}
}

func (g *Game) drawTracks(dst *ebiten.Image) {
	for i, pts := range g.paths {
		c := g.cols[i]
		w := float32(1.3)
		if g.kinds[i] == KindMuon {
			w = 2.2
		}
		for j := 1; j < len(pts); j++ {
			x0, y0 := toScreen(pts[j-1].X, pts[j-1].Y)
			x1, y1 := toScreen(pts[j].X, pts[j].Y)
			vector.StrokeLine(dst, x0, y0, x1, y1, w, c, true)
		}
	}
	vector.StrokeCircle(dst, float32(viewCX), float32(viewCY), 3, 2, color.RGBA{255, 255, 255, 255}, true)
}

func (g *Game) drawMET(dst *ebiten.Image) {
	if g.reco.MET < 5 {
		return
	}
	phi := math.Atan2(g.reco.METy, g.reco.METx)
	l := math.Min(RMuOut, 1.0+g.reco.MET*0.06)
	x0, y0 := toScreen(0, 0)
	x1, y1 := toScreen(l*math.Cos(phi), l*math.Sin(phi))
	vector.StrokeLine(dst, x0, y0, x1, y1, 3, color.RGBA{255, 255, 255, 220}, true)
	// 矢じり
	for _, s := range []float64{2.6, -2.6} {
		a := phi + s
		xa, ya := toScreen(l*math.Cos(phi)+0.35*math.Cos(a), l*math.Sin(phi)+0.35*math.Sin(a))
		vector.StrokeLine(dst, x1, y1, xa, ya, 3, color.RGBA{255, 255, 255, 220}, true)
	}
}

func (g *Game) drawPanel(dst *ebiten.Image) {
	x0 := 790
	var b strings.Builder

	fmt.Fprintf(&b, "PROCESS  %s\n", processNames[g.proc])
	fmt.Fprintf(&b, "B field  %.1f T   (solenoid, +z)\n\n", g.B)
	b.WriteString("1-6 process / SPACE new event / F +500 events\n")
	b.WriteString("C clear hist / UP-DOWN B field\n\n")

	b.WriteString("--- TRUTH ---\n")
	if g.ev.TruthMass > 0 {
		fmt.Fprintf(&b, " generated mass   %8.3f GeV\n", g.ev.TruthMass)
	}
	fmt.Fprintf(&b, " final particles  %8d\n", len(g.ev.Particles))
	fmt.Fprintf(&b, " true MET (nu)    %8.2f GeV\n\n", g.reco.TruthMET)

	b.WriteString("--- RECONSTRUCTED ---\n")
	fmt.Fprintf(&b, " tracks (pT>0.4)  %8d\n", g.reco.NTracks)
	fmt.Fprintf(&b, " sum ET           %8.1f GeV\n", g.reco.SumEt)
	fmt.Fprintf(&b, " MET              %8.2f GeV\n", g.reco.MET)
	if g.reco.MassOK {
		fmt.Fprintf(&b, " invariant mass   %8.3f GeV\n", g.reco.InvMass)
	}
	if g.reco.MTOK {
		fmt.Fprintf(&b, " transverse mass  %8.3f GeV\n", g.reco.MT)
	}
	b.WriteString("\n leading objects:\n")
	shown := 0
	for _, o := range g.reco.Objects {
		if o.Label == "soft" || o.Label == "had" || shown >= 4 {
			continue
		}
		fmt.Fprintf(&b, "  %-8s pT=%7.2f  eta=%+5.2f  phi=%+5.2f\n",
			o.Label, o.P.Pt(), o.P.Eta(), o.P.Phi())
		shown++
	}
	for i, j := range g.reco.Jets {
		if i >= 3 {
			break
		}
		fmt.Fprintf(&b, "  jet%d     pT=%7.2f  eta=%+5.2f  phi=%+5.2f\n",
			i+1, j.Pt(), j.Eta(), j.Phi())
	}

	ebitenutil.DebugPrintAt(dst, b.String(), x0, 14)

	// ---- ヒストグラム ----
	hx, hy, hw, hh := float32(x0), float32(430), float32(455), float32(250)
	vector.StrokeRect(dst, hx, hy, hw, hh, 1, color.RGBA{70, 80, 100, 255}, false)
	h := g.hist
	maxN := float64(h.Max())
	bw := hw / float32(len(h.Bins))
	for i, n := range h.Bins {
		if n == 0 {
			continue
		}
		bh := float32(float64(n) / maxN * float64(hh-4))
		vector.DrawFilledRect(dst, hx+float32(i)*bw+1, hy+hh-bh, bw-1, bh,
			color.RGBA{100, 180, 255, 220}, false)
	}
	info := fmt.Sprintf("%s\nentries %d   mean %.2f   RMS %.2f\nrange [%.0f, %.0f]",
		h.Title, h.N, h.Mean(), h.RMS(), h.Lo, h.Hi)
	ebitenutil.DebugPrintAt(dst, info, x0, 690)
}

func (g *Game) Layout(w, h int) (int, int) { return screenW, screenH }

func main() {
	ebiten.SetWindowSize(screenW, screenH)
	ebiten.SetWindowTitle("Collider Event Display")
	if err := ebiten.RunGame(NewGame()); err != nil {
		panic(err)
	}
}
