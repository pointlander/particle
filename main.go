package main

import (
	"bytes"
	"encoding/binary"
	"math"
	"math/rand"

	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/plotutil"
	"github.com/gonum/plot/vg"
	"github.com/pointlander/compress"
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

type Particle struct {
	Mass, Size, X, Y, VX, VY float64
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

	a.VX -= coeff1 * (a.X - b.X) * (2.0 * b.Mass) / (a.Mass + b.Mass)
	a.VY -= coeff1 * (a.Y - b.Y) * (2.0 * b.Mass) / (a.Mass + b.Mass)
	b.VX -= coeff2 * (b.X - a.X) * (2.0 * a.Mass) / (a.Mass + b.Mass)
	b.VY -= coeff2 * (b.Y - a.Y) * (2.0 * a.Mass) / (a.Mass + b.Mass)
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

func main() {
	particles := make([]Particle, Num)
	complexity := make([]int, Steps)
	velocity := computeV(Temp, Mass)
	for i := range particles {
		particles[i].Mass = Mass
		particles[i].Size = Size
		particles[i].X = rand.Float64() * 0.95 * BoxX
		particles[i].Y = rand.Float64() * 0.95 * BoxY

		theta := rand.Float64() * 2.0 * math.Pi
		particles[i].VX = velocity * math.Cos(theta)
		particles[i].VY = velocity * math.Sin(theta)
	}

	dt := selectDt(&particles[0])
	in, out := &bytes.Buffer{}, &bytes.Buffer{}
	write := func(n float64) {
		bits := math.Float64bits(n)
		bytes := make([]byte, 8)
		binary.LittleEndian.PutUint64(bytes, bits)
		in.Write(bytes)
	}
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
			enforceWallsPeriodic(a, dt)
			move(a, dt)
		}
		if s == 2000 {
			graphVelocityDistribution(particles, "velocities_start.png")
		}

		in.Reset()
		for _, particle := range particles {
			write(particle.Mass)
			write(particle.Size)
			write(particle.X)
			write(particle.Y)
			//write(math.Sqrt(particle.VX*particle.VX + particle.VY*particle.VY))
			//write(math.Atan2(particle.VY, particle.VX))
			write(particle.VX)
			write(particle.VY)
		}
		out.Reset()
		input := make(chan []byte, 1)
		input <- in.Bytes()
		close(input)
		compress.BijectiveBurrowsWheelerCoder(input).MoveToFrontRunLengthCoder().AdaptiveCoder().Code(out)
		complexity[s] = out.Len()
	}
	graphVelocityDistribution(particles, "velocities_end.png")
	graphComplexity(complexity)
}
