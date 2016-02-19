package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

func dist(x1, y1, z1, x2, y2, z2, ex1, ey1, ez1, ex2, ey2, ez2, leng float64, ch chan [2]float64) {
	var xmu0, xla0 float64
	// distance vector between centers:
	x12 := x2 - x1
	y12 := y2 - y1
	z12 := z2 - z1
	// 	First, find the normal vector en at minimal distance between the
	//  carrier lines g1 and g2:
	//  Since en is normal both to e1 and e2 we have for any normal vector
	//  r12*e1=lam-mu*(e1*e2)
	//  r12*e2=lam*(e1*e2)-mu
	e1r := ex1*x12 + ey1*y12 + ez1*z12
	e2r := ex2*x12 + ey2*y12 + ez2*z12
	e12 := ex1*ex2 + ey1*ey2 + ez1*ez2
	// 	if e1 is parallel to e2 (e1*e2=1), there are infinitely many solutions, just pick one. \
	//  if e1 is not parallel to e2 (e1*e2!=1), there is only one solution
	//  lam0=[(e1*r12)-(e1*e2)(e2*r12)]/[1-(e1*e2)^2]
	//  mu0=-[(e2*r12)-(e1*e2)(e1*r12)]/[1-(e1*e2)^2]
	if math.Abs(math.Abs(e12)-1) < 1e-8 {
		xmu0 = 0
		xla0 = e1r
	} else {
		recipr := 1.0 / (1.0 - e12*e12)
		xla0 = recipr * (e1r - e12*e2r)
		xmu0 = -recipr * (e2r - e12*e1r)
	}
	// Shortest perpendicular distance between carrier lines:
	xx12 := x12 + xmu0*ex2 - xla0*ex1
	yy12 := y12 + xmu0*ey2 - xla0*ey1
	zz12 := z12 + xmu0*ez2 - xla0*ez1
	rnsq := xx12*xx12 + yy12*yy12 + zz12*zz12
	// 	Now for the squared in-plane distance risq between two line segments:
	//  Rectangle half lengths h1=L1/2, h2=L2/2
	h1 := 0.5 * leng
	h2 := 0.5 * leng
	// 	If the origin is contained in the rectangle, life is easy: the origin is the minimum, and
	//  the in-plane distance is zero!
	var risq float64
	if math.Abs(xla0) <= h1 && math.Abs(xmu0) <= h2 {
		risq = 0.0
	} else {
		// Find minimum of f=gamma^2+delta^2-2*gamma*delta*(e1*e2)
		// where gamma, delta are the line parameters reckoned from the intersection (=lam0,mu0)
		// First, find the lines gamm and delm that are nearest to the origin:
		gamm := math.Min(-xla0-h1, -xla0+h1)
		delm := math.Min(-xmu0-h2, -xmu0+h2)
		//  Now choose the line gamma=gamm and optimize delta:
		delms := gamm * e12
		// look if delms is within [-xmu0+/-L/2]:
		var del float64
		switch a := xmu0 + delms; {
		case a >= (-h1) && a <= h1:
			del = delms
		case a < (-h1):
			del = -xmu0 - h1
		case a > h1:
			del = -xmu0 + h1
		}
		// Distance at these gam, del:
		f1 := gamm*gamm + del*del - 2*gamm*del*e12
		// Now choose the line delta=deltam and optimize gamma:
		gamms := delm * e12
		//  look if gamms is within [-xla0+/-L/2]:
		var gam float64
		switch a := xla0 + gamms; {
		case a >= (-h1) && a <= h1:
			gam = gamms
		case a < (-h1):
			gam = -xla0 - h1
		case a > h1:
			gam = -xla0 + h1
		}
		//  Distance at these gam, del:
		f2 := gam*gam + delm*delm - 2*gam*delm*e12
		//  Compare f1 and f2 to find risq:
		risq = math.Min(f1, f2)
	}
	rclsq := rnsq + risq
	rcl := math.Sqrt(rclsq)
	res := [2]float64{rcl, math.Sqrt(1 - e12*e12)}
	ch <- res
}

func generator(r, leng, boxleng, vol float64) (x, y, z, ex, ey, ez []float64) {
	rodvol := math.Pi*(r*r)*leng + 4/3*math.Pi*(r*r*r)
	total := int(boxleng * boxleng * boxleng * vol / rodvol)
	x = make([]float64, total) // x,y,z store the center of each rod
	y = make([]float64, total)
	z = make([]float64, total)
	ex = make([]float64, total) // ex,ey,ez store the orientation of each rod
	ey = make([]float64, total)
	ez = make([]float64, total)
	var S, U, V, V1, V2 float64
	for i := 0; i < total; i++ {
		x[i] = rand.Float64() * boxleng
		y[i] = rand.Float64() * boxleng
		z[i] = rand.Float64() * boxleng
		S = 2
		for S >= 1 {
			U = rand.Float64()
			V = rand.Float64()
			V1 = (U - 0.5) * 2
			V2 = (V - 0.5) * 2
			S = V1*V1 + V2*V2
		}
		ex[i] = 2 * V1 * math.Sqrt(1-S)
		ey[i] = 2 * V2 * math.Sqrt(1-S)
		ez[i] = 1 - 2*S
	}
	return
}

func main() {
	rand.Seed(time.Now().UTC().UnixNano())
	r := 1.6
	leng := 800.0
	boxleng := 5000.0
	vol := 0.001
	x, y, z, ex, ey, ez := generator(r, leng, boxleng, vol)
	ch := make(chan [2]float64, 100)
	count := 0
	for i := 0; i < len(x); i++ {
		for j := i + 1; j < len(x); j++ {
			if (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]) <= (leng+2*8.6)*(leng+2*8.6) {
				go dist(x[i], y[i], z[i], x[j], y[j], z[j], ex[i], ey[i], ez[i], ex[j], ey[j], ez[j], leng, ch)
				count++
			}
		}
	}
	dists := make([]float64, count)
	sinbeta := make([]float64, count)
	var res [2]float64
	for i := 0; i < count; i++ {
		res = <-ch
		dists[i] = res[0]
		sinbeta[i] = res[1]
	}
	fmt.Println(len(dists))
}
