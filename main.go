// Copyright 2016 The Particle Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"encoding/binary"
	"flag"
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"sort"

	"github.com/pointlander/compress"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

const (
	// Num number of particles
	Num = 100
	// Steps number of simulation steps
	Steps = 20000
	// Size of unit i.e. particle (charge radius of proton in cm)
	Size = 8.775e-14
	// Mass (mass of proton in g for simplicity)
	Mass = 1.6726219e-24
	// Temp temperature of gas (kelvin)
	Temp = 273.0
	// BoxX box dimensions (cm/r_proton)
	BoxX = 50.0 * Size
	// BoxY box dimensions (cm/r_proton)
	BoxY = 50.0 * Size
	// KB boltzmann constant (ergs/K)
	KB = 1.380658e-16
)

// Particle is a particle with mass
type Particle struct {
	Mass, Size, X, Y, VX, VY float64
}

func (p *Particle) morton() (code [16]byte) {
	state := [...]uint32{
		fixed32(p.X), fixed32(p.Y),
		fixed32(p.VX), fixed32(p.VY),
	}
	length := 8 * len(code)
	for i := 0; i < length; i++ {
		code[15-i/8] >>= 1
		code[15-i/8] |= byte(state[i%4]&1) << 7
		state[i%4] >>= 1
	}
	return
}

// Particles is a gas of particles
type Particles []Particle

func (p Particles) Len() int {
	return len(p)
}

func (p Particles) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

func (p Particles) Less(i, j int) bool {
	a, b := big.Int{}, big.Int{}
	ab, bb := p[i].morton(), p[j].morton()
	a.SetBytes(ab[:])
	b.SetBytes(bb[:])
	return a.Cmp(&b) < 0
}

func fixed32(a float64) uint32 {
	return uint32(a*65536 + 0.5)
}

func dot2D(x1, y1, x2, y2 float64) float64 {
	return x1*x2 + y1*y2
}

func didCollide(a, b *Particle) float64 {
	dx, dy := a.X-b.X, a.Y-b.Y
	distance := dx*dx + dy*dy
	if distance < 4.0*a.Size*a.Size {
		return distance
	}
	return -1
}

func rollbackTime(distance float64, a, b *Particle) float64 {
	vx, vy := a.VX-b.VX, a.VY-b.VY
	va := math.Sqrt(vx*vx + vy*vy)
	return (2.0*a.Size - distance) / va
}

func selectDt(a *Particle) float64 {
	v := math.Sqrt(a.VX*a.VX + a.VY*a.VY)
	return 0.33 * a.Size / v
}

func enforceWallsPeriodic(a *Particle, dt float64) {
	if a.X > BoxX && a.Y > BoxY {
		//Outside top right corner
		a.X = a.X - BoxX
		a.Y = a.Y - BoxY
	} else if a.X < 0.0 && a.Y > BoxY {
		//Outside top left corner
		a.X = a.X + BoxX
		a.Y = a.Y - BoxY
	} else if a.X < 0.0 && a.Y < 0.0 {
		//Outside bottom left corner
		a.X = a.X + BoxX
		a.Y = a.Y + BoxY
	} else if a.X > BoxX && a.Y < 0.0 {
		//Outside bottom right corner
		a.X = a.X - BoxX
		a.Y = a.Y + BoxY
	} else if a.X > BoxX {
		//Outside right wall
		a.X = a.X - BoxX
	} else if a.X < 0.0 {
		//Outside left wall
		a.X = a.X + BoxX
	} else if a.Y > BoxY {
		//Outside top wall
		a.Y = a.Y - BoxY
	} else if a.Y < 0.0 {
		//Outside bottom wall
		a.Y = a.Y + BoxY
	}
}

func collision(a, b *Particle) {
	coeff1 := dot2D(a.VX-b.VX, a.VY-b.VY, a.X-b.X, a.Y-b.Y)
	coeff1 /= dot2D(a.X-b.X, a.Y-b.Y, a.X-b.X, a.Y-b.Y)
	coeff2 := dot2D(b.VX-a.VX, b.VY-a.VY, b.X-a.X, b.Y-a.Y)
	coeff2 /= dot2D(b.X-a.X, b.Y-a.Y, b.X-a.X, b.Y-a.Y)

	aMass, bMass := math.Abs(a.Mass), math.Abs(b.Mass)

	a.VX -= coeff1 * (a.X - b.X) * (2.0 * b.Mass) / (aMass + bMass)
	a.VY -= coeff1 * (a.Y - b.Y) * (2.0 * b.Mass) / (aMass + bMass)
	b.VX -= coeff2 * (b.X - a.X) * (2.0 * a.Mass) / (aMass + bMass)
	b.VY -= coeff2 * (b.Y - a.Y) * (2.0 * a.Mass) / (aMass + bMass)
}

func move(a *Particle, dt float64) {
	a.X += dt * a.VX
	a.Y += dt * a.VY
}

func computeV(temp, mass float64) float64 {
	return math.Sqrt(3.0 * KB * temp / mass)
}

func graphVelocityDistribution(particles []Particle, filename string) {
	velocities, min, max := make([]float64, Num), math.MaxFloat64, 0.0
	for i, particle := range particles {
		velocity := math.Sqrt(particle.VX*particle.VX + particle.VY*particle.VY)
		if velocity < min {
			min = velocity
		}
		if velocity > max {
			max = velocity
		}
		velocities[i] = velocity
	}
	width := (max - min) / 10
	values := make(plotter.Values, 10)
	for _, velocity := range velocities {
		bucket := min + width
		for i := range values {
			if velocity < bucket {
				values[i]++
				break
			}
			bucket += width
		}
	}

	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = "Velocity Distribution"
	p.Y.Label.Text = "Particles"

	w := vg.Points(20)

	bars, err := plotter.NewBarChart(values, w)
	if err != nil {
		panic(err)
	}
	bars.LineStyle.Width = vg.Length(0)
	bars.Color = plotutil.Color(0)

	p.Add(bars)

	if err := p.Save(5*vg.Inch, 3*vg.Inch, filename); err != nil {
		panic(err)
	}
}

func graphComplexity(complexity []int) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	points := make(plotter.XYs, len(complexity))
	for x, y := range complexity {
		points[x].X = float64(x)
		points[x].Y = float64(y)
	}

	s, err := plotter.NewScatter(points)
	if err != nil {
		panic(err)
	}
	s.Radius = vg.Length(1)

	p.Add(s)
	if err := p.Save(16*vg.Inch, 16*vg.Inch, "complexity.png"); err != nil {
		panic(err)
	}

	/*p, err = plot.New()
	if err != nil {
		panic(err)
	}

	l, err := plotter.NewLine(points[len(points)-100:])
	if err != nil {
		panic(err)
	}
	l.LineStyle.Width = vg.Points(1)

	p.Add(l)
	if err := p.Save(16*vg.Inch, 16*vg.Inch, "complexity_zoom.png"); err != nil {
		panic(err)
	}*/
}

var (
	negative    = flag.Bool("negative", false, "simulate negative mass")
	kcomplexity = flag.Bool("kcomplexity", false, "compute k complexity")
)

func main() {
	flag.Parse()

	particles := make(Particles, Num)
	complexity := make([]int, Steps)
	velocity := computeV(Temp, Mass)
	for i := range particles {
		particles[i].Mass = Mass
		if *negative && i > Num/2 {
			particles[i].Mass = -Mass
		}
		particles[i].Size = Size
		particles[i].X = rand.Float64() * 0.95 * BoxX
		particles[i].Y = rand.Float64() * 0.95 * BoxY

		theta := rand.Float64() * 2.0 * math.Pi
		particles[i].VX = velocity * math.Cos(theta)
		particles[i].VY = velocity * math.Sin(theta)
	}

	delta := selectDt(&particles[0])
	in, out := &bytes.Buffer{}, &bytes.Buffer{}
	write := func(n float64) {
		bits := math.Float64bits(n)
		bytes := make([]byte, 8)
		binary.LittleEndian.PutUint64(bytes, bits)
		in.Write(bytes)
	}
	_ = write
	for s := 0; s < Steps; s++ {
		for i := range particles {
			a := &particles[i]
			for j := range particles {
				b := &particles[j]
				if distance := didCollide(a, b); i != j && distance > 0.0 {
					dt := rollbackTime(math.Sqrt(distance), a, b)
					move(a, -dt)
					move(b, -dt)

					collision(a, b)

					move(a, dt)
					move(b, dt)
				}
			}
			enforceWallsPeriodic(a, delta)
			move(a, delta)
		}
		if s == 2000 {
			graphVelocityDistribution(particles, "velocities_start.png")
		}

		in.Reset()
		sort.Sort(particles)
		a := particles[0].morton()
		last := big.Int{}
		last.SetBytes(a[:])
		for _, particle := range particles[1:] {
			/*write(particle.Mass)
			write(particle.Size)
			write(particle.X)
			write(particle.Y)*/
			//write(math.Sqrt(particle.VX*particle.VX + particle.VY*particle.VY))
			//write(math.Atan2(particle.VY, particle.VX))
			/*write(particle.VX)
			write(particle.VY)*/
			b := particle.morton()
			current, diff := big.Int{}, big.Int{}
			current.SetBytes(b[:])
			diff.Sub(&current, &last)
			if diff.Sign() < 0 {
				panic("diff is negative")
			}
			in.Write(diff.Bytes())
			last = current
		}
		out.Reset()
		if *kcomplexity {
			input := make(chan []byte, 1)
			input <- in.Bytes()
			close(input)
			compress.BijectiveBurrowsWheelerCoder(input).MoveToFrontRunLengthCoder().AdaptiveCoder().Code(out)
			complexity[s] = out.Len()
		}
	}
	graphVelocityDistribution(particles, "velocities_end.png")
	if *kcomplexity {
		graphComplexity(complexity)
	}

	meanV, meanVPositive, meanVNegative, countPositive, countNegative :=
		0.0, 0.0, 0.0, 0, 0
	for i := range particles {
		vx, vy := particles[i].VX, particles[i].VY
		v := math.Sqrt(vx*vx + vy*vy)
		meanV += v
		if particles[i].Mass > 0 {
			meanVPositive += v
			countPositive++
		} else {
			meanVNegative += v
			countNegative++
		}
	}
	fmt.Printf("mean velocity=%v\n", meanV/float64(len(particles)))
	fmt.Printf("mean positive velocity=%v\n", meanVPositive/float64(countPositive))
	if *negative {
		fmt.Printf("mean negative velocity=%v\n", meanVNegative/float64(countNegative))
	}
}
