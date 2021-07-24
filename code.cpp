/*

Fire missile assuming quadratic air resistance through isothermal atmosphere on a uniform spherical spinning Earth.
In reality, Earth is more of an oblate spheroid, and the layers of atmosphere are defined by how their temperatures change.
Air density and g are calculated as a function of altitude.

Does not simulate any thrust or change in projectile mass.
To make code more realistic, maybe set initial altitude (alt0) to be above Earth's atmosphere
to correspond to after the multi-stage thrust.

The goal is to minimize a cost function to hit a target via the optimal initial velocity vector.
To do this, I assume that launch speed is adjustable.
To hit the target, I use gradient-descent optimization of cost(),
a function that feels a bit hacked together. Please tinker with and fine tune cost()!
Then, I vary initial velocity by 0.1% in random directions to get a sense of where the rocket would actually hit given
random winds, random thrust, random wobbling, etc.

Runge-Kutta method (RK4) is used with adaptive step size proportional to velocity/acceleration.
This ensures that the dv/v for each time step is constant since
  dt = (dv/v) * v/a
I want to keep dv/v the same because you can use the simple projectile range formula from freshman physics to show that
  dRange / Range is proportional to dv / v
but maybe you'll want to do some other type of adaptive step size.
An adaptive time step is important when k is large and v changes rapidly in the beginning.

All units in this file are MKS (and sometimes angular degrees).

I give most physical variables global scope!



Inputs positions are in spherical coordinates.
Calculations are in rectangular coordinates with center of Earth as origin.
North pole is on the +z axis.
To convert...
  x = r sin(th) cos(phi)
  y = r sin(th) sin(phi)
  z = r cos(th)
Launch position has y=0 (so phi=0).
Acceleration of gravity is
  - G M / r^2 rHat
where
  r = sqrt(x^2 + y^2 + z^2)
  rHat = sin(th) cos(phi) xHat + sin(th) sin(phi) yHat + cos(th) zHat
  phi = atan2(y, x)
  th = atan2(sqrt(x^2 + y^2), z)
To get air-resistance acceleration, use
  - k exp( - (r-re) / 10400.0 ) v^2 vHat
where
  re = 6.27e6
All coordinates are attached to the rotating Earth,
so Newton's 2nd law gains some non-inertial (Coriolis and centrifugal) terms...
  a = F/m + 2 v*Omega + (Omega*r)*Omega
where "*" is a vector cross product, and Omega = 7.292e-5 zHat
for a sidereal day of 86164 s.
  a = F/m + 2 (vy xHat - vx yHat) Omega + (x xHat + y yHat) Omega^2


(c) 2020 Bradley Knockel
*/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;


#define piOver180 1.74532925199432957e-2
#define sixth 0.16666666666666667
#define GM 3.986e14        // G * M_earth
#define re 6.27e6          // average radius of Earth
#define Omega 7.292e-5     // angular velocity of Earth
#define oneOverH 9.615e-5  // 1/H; H is height scale of air's exponential density fall


// set some constants
const double  k = 1e-4;        // strength of air resistance at sea level (depends on projectile)
const double alt0 = 0.0;       // start altitude (0 is sea level)
const double lat0 = 0.0;       // start latitude in degrees (-90 to 90)
const double altEnd = 0.0;     // target altitude (0 is sea level)
const double latEnd = 90.0;    // target latitude in degrees (-90 to 90 from south pole to north pole)
const double longEnd = 180.0;    // target longitude in degrees (-180 to 180; start longitude is 0)
const double dvOverV = 0.01;   // desired fractional change of v each time step
const double tLimit = 1e6;     // time limit in case of orbits
const int iStop = 0.10000001 / dvOverV;      // how many steps to wait to output to file
const string dataFile = "data.dat";    // name of output file for trajectory
const string dataFile2 = "data2.dat";  // name of output file when varying v0



//   alt0, lat0   -->   r0, th0  (phi0=0)   -->   x0, z0  (y0=0)
const double r0 = re + alt0;
const double th0 = (90 - lat0) * piOver180 ;
const double x0 = r0 * sin( th0 );
const double z0 = r0 * cos( th0 );

//   altEnd, latEnd   -->   rEnd, thEnd, phiEnd
const double rEnd = re + altEnd;
const double thEnd = (90 - latEnd) * piOver180 ;
const double phiEnd = longEnd * piOver180;

// target rectangular coordinates
const double xEnd = rEnd * sin(thEnd) * cos(phiEnd);
const double yEnd = rEnd * sin(thEnd) * sin(phiEnd);
const double zEnd = rEnd * cos(thEnd);

// Reference distance, speed, and time
const double dRef = sqrt( r0*r0 + rEnd*rEnd - 2*r0*rEnd*           // straight-line distance to target
    (sin(th0)*sin(thEnd)*cos(phiEnd) + cos(th0)*cos(thEnd)) );
const double vRef = sqrt( GM/(re*re) * dRef); // since range is proportional to v0^2 / g
const double tRef = vRef * re*re / GM;        // since time of flight is proportional to v0 / g



// declare
double g,kz,dt;
double r,th,phi,vr,alt,dAlt;
double x,y,z,vx,vy,vz,v,t;
double v0,v0x,v0y,v0z;
ofstream myfile(dataFile);
ofstream myfile2(dataFile2);




// --------------- cost() -----------------------
// 
// To hit the target with the lowest launch speed in the shortest time.
// Lowest launch speed could save fuel or allow for more maneuverability with unused fuel.
// Shortest time is good because target could move or missile could be shot down.
//
// One should minimize this cost to find the best trajectory.
// To find a gradient of cost(), I first made sure that RK4() interpolates back to hit altEnd exactly.
//
// This should be a smooth function so that finding gradients is a useful approach to minimizing it.
//
// I think it's important to hit the target THEN worry about the rest,
//   so I multiply the targeting term by 100.
// Large k broke it because large changes in v0 are needed to slightly change where missile lands,
//   so I put a sqrt() around the v0 term.
// Very large k still breaks it since the v0 term dominates.



double cost(){
  return 100*sqrt(pow(xEnd-x,2)+pow(yEnd-y,2)+pow(zEnd-z,2))/dRef + sqrt(v0/vRef) + t / tRef;
}





// --------------- Runge-Kutta method: RK4()  --------------------

// Each of the following arrays has the form...
//     {x, y, z, vx, vy, vz}
//      0  1  2  3   4   5
double vars1[6];
double vars2[6];  // to temporarily store things
double k1[6],k2[6],k3[6],k4[6];


// Returns time derivative of v1 as v2
void derivative(double v1[], double v2[]) {

  phi = atan2(v1[1], v1[0]);
  th = atan2(sqrt( pow(v1[0],2) + pow(v1[1],2) ), v1[2]);
  r = sqrt( pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2) );
  v = sqrt( pow(v1[3],2) + pow(v1[4],2) + pow(v1[5],2) );
  g = GM / (r*r);
  kz = k * exp(-(r-re) * oneOverH);  // isothermal Earth atmosphere

  // derivatives of positions are velocities
  v2[0] = v1[3];
  v2[1] = v1[4];
  v2[2] = v1[5];

  // accelerations due to quadratic air resistance and gravity (including non-inertial)
  v2[3] = - kz * v * v1[3] - g * sin(th) * cos(phi) + Omega * Omega * v1[0] + 2 * Omega * v1[4];
  v2[4] = - kz * v * v1[4] - g * sin(th) * sin(phi) + Omega * Omega * v1[1] - 2 * Omega * v1[3];
  v2[5] = - kz * v * v1[5] - g * cos(th);

}


// linear combination of arrays; saves to vars2
void vAdd(double a[], double b[], double c) {
  for (int k=0; k<6; k++)
    vars2[k] = a[k] + b[k]*c;
}


void RK4(bool print, double v0x, double v0y, double v0z) {

  // initialize
  vars1[0] = x0;
  vars1[1] = 0.0;  // y0
  vars1[2] = z0;
  vars1[3] = v0x;
  vars1[4] = v0y;
  vars1[5] = v0z;
  alt = alt0;
  vr = sin(th0)*v0x + cos(th0)*v0y;

  // compute
  t = 0.0;
  int i = 0;
  while( t<tLimit && alt<1.0e8 && (vr > 0.0 || alt>altEnd) ){

    // calculate k1, k2, k3, and k4
    derivative(vars1, k1);    // create k1 from vars1   
    dt = dvOverV * v / sqrt( pow(k1[3],2) + pow(k1[4],2) + pow(k1[5],2) );  // adaptive time step
    vAdd(vars1, k1, dt*0.5);  // create vars2
    derivative(vars2, k2);    // create k2 from vars2
    vAdd(vars1, k2, dt*0.5);  // create vars2
    derivative(vars2, k3);    // create k3 from vars2
    vAdd(vars1, k3, dt);      // create vars2
    derivative(vars2, k4);    // create k4 from vars2

    // update vars1 and t
    for (int j=0; j<6; j++) {
      vars2[j] = dt * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]) * sixth;
      vars1[j] += vars2[j];
    }
    t += dt;

    // update some spherical coordinates that are needed
    alt = sqrt( pow(vars1[0],2) + pow(vars1[1],2) + pow(vars1[2],2) ) - re;
    phi = atan2(vars1[1], vars1[0]);
    th = atan2(sqrt( pow(vars1[0],2) + pow(vars1[1],2) ), vars1[2]);
    vr = sin(th)*cos(phi)*vars1[3] + sin(th)*sin(phi)*vars1[4] + cos(th)*vars1[5] ;  // dot product of rHat and v
    double dAlt = sin(th)*cos(phi)*vars2[0] + sin(th)*sin(phi)*vars2[1] + cos(th)*vars2[2] ;

    // linearly interpolate back to get the missile to hit at exactly altEnd
    if ( alt < altEnd && vr < 0.0 && (altEnd-alt) <= -dAlt ) {
      double fraction = (altEnd - alt) / dAlt;  // a negative number
      for (int j=0; j<6; j++)
        vars1[j] += fraction * vars2[j];
      t += fraction * dt;
      break;
    }

    // print to file?
    i++;
    if ( print && i == iStop ) {
      i = 0;
      myfile << t << ' ' << vars1[0] << ' ' << vars1[1] << ' ' << vars1[2] 
                  << ' ' << vars1[3] << ' ' << vars1[4] << ' ' << vars1[5]
                  << ' ' << alt << ' ' << re*acos(sin(lat0)*cos(th) + cos(lat0)*sin(th)*cos(phi)) << endl;
    }

  }

  // unpack the final values; these, and t, are outputs of RK4()
  x  = vars1[0];
  y  = vars1[1];
  z  = vars1[2];
  vx = vars1[3];
  vy = vars1[4];
  vz = vars1[5];

  // create more outputs of RK4()
  v0 = sqrt( v0x*v0x + v0y*v0y + v0z*v0z );
  v = sqrt( vx*vx + vy*vy + vz*vz );
  r = sqrt( x*x + y*y + z*z );
  alt = r-re;
  phi = atan2(y, x);
  th = atan2(sqrt( x*x + y*y ), z );

}





// ---------- gradient descent optimization: optimize() -------------


// initialize multiplier when taking steps along gradient
double m = 100000000.0;

const double mStop = 0.00001;  // smallest m before stopping
const double f = 0.000000001;      // fraction of v0 to change when finding gradient


////// Find optimal velocity using basic "Gradient descent" method
// That is, choose best v0x, v0y, and v0z to minimize cost()
void optimize(){

  // set initial guess (fire straight up)
  v0x = vRef*sin(th0);
  v0y = 0.0;
  v0z = vRef*cos(th0);

  // get initial cost
  RK4(false,v0x,v0y,v0z);
  double b = cost();


  while(true) {

    double grad[3];      // gradient of cost function
    double deltaV = f*v0;

    // calculate gradient of cost function
    RK4(false, v0x+deltaV, v0y, v0z);
    grad[0] = (cost() - b) / deltaV;
    RK4(false, v0x, v0y+deltaV, v0z);
    grad[1] = (cost() - b) / deltaV;
    RK4(false, v0x, v0y, v0z+deltaV);
    grad[2] = (cost() - b) / deltaV;

    // double m
    m *= 2.0;

    // Take a step in direction of negative gradient.
    // Reduce m until this step actually minimizes b
    while (true) {

      v0x -= grad[0] * m;
      v0y -= grad[1] * m;
      v0z -= grad[2] * m;
      RK4(false,v0x,v0y,v0z);
      double bTemp = cost();

      // don't allow negative v0r
      // don't allow cost to increase
      // don't allow missile to miss in r direction
      double v0r = sin(th0)*v0x + cos(th0)*v0y;
      if (v0r < 0 || bTemp >= b || abs(altEnd-alt)/dRef > 0.00001) {
        v0x += grad[0] * m;
        v0y += grad[1] * m;
        v0z += grad[2] * m;
        m *= 0.5;  // half m
      } else {
        b = bTemp;
        break;
      }

      if (m < mStop)
        break;

    }

    cout << "    v0x = " << v0x
         << "  v0y = " << v0y
         << "  v0z = " << v0z 
         << "    b = " << b 
         << "  m = " << m  << endl;

    if (m < mStop)
      break;

  }
}




// --------------- varyV0() --------------------------------------


double oldV0x,oldV0y,oldV0z,oldX,oldY,oldZ,thOld,phiOld;
double newV0x,newV0y,newV0z;

// deviation from v0 landing location along old final phiHat, -thHat, rHat respectively
double missX,missY,missZ;

void varyV0(int iStop) {

  // 0.1% of v0
  double radius = 0.001 * v0;

  // get the "correct" values before randomly changing them
  oldV0x = v0x;
  oldV0y = v0y;
  oldV0z = v0z;
  oldX = x;
  oldY = y;
  oldZ = z;
  thOld = th;
  phiOld = phi;

  for(int i=0; i<=iStop; i++) {

    // vary the v0 in a way that favors no direction to get newV0x,newV0y,newV0z
    while(true) {

      // generate random numbers between -1 and 1 (uniform distribution)
      newV0x = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2.0) ) - 1.0;
      newV0y = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2.0) ) - 1.0;
      newV0z = static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/2.0) ) - 1.0;

      double temp = newV0x*newV0x + newV0y*newV0y + newV0z*newV0z;
      if ( temp < 1 ) {              // does coordinate fall in sphere of radius = 1 ?
        temp = radius / sqrt(temp);  // scale factor to get to desired radius
        newV0x = oldV0x + temp * newV0x;
        newV0y = oldV0y + temp * newV0y;
        newV0z = oldV0z + temp * newV0z;
        break;
      }

    }

    // launch
    RK4(false,newV0x,newV0y,newV0z);

    // calculate missX,missY,missZ by projecting delta(x) vector onto old phiHat,-thHat,rHat
    missX = - sin(phiOld)*(x-oldX) + cos(phiOld)*(y-oldY);
    missY = - cos(thOld)*cos(phiOld)*(x-oldX) - cos(thOld)*sin(phiOld)*(y-oldY) + sin(thOld)*(z-oldZ) ;
    missZ =   sin(thOld)*cos(phiOld)*(x-oldX) + sin(thOld)*sin(phiOld)*(y-oldY) + cos(thOld)*(z-oldZ);

    // save to file
    myfile2 << newV0x << ' ' << newV0y << ' ' << newV0z << ' '
            << missX  << ' ' << missY  << ' ' << missZ  << endl;

  }

}





// --------------- main() --------------------------------------


int main() {


  // run RK4() many times to optimize v0x, v0y, and v0z to minimize cost()
  cout << "  optimizing trajectory..." << endl;
  optimize();



  /*

  // Use with alt0>0, lat0=0, and k=0 to get circular orbits
  // Increase v0y to get elliptical (won't look elliptical unless Omega = 0)
  v0x = 0.0;
  v0y = sqrt(GM/r0) - Omega*r0;
  v0z = 0.0;

  // use with lat0=0 and k=0 to get escape velocity
  v0x = 0.0;
  v0y = sqrt(2*GM/r0) - Omega*r0;
  v0z = 0.0;

  // initial velocity for some nice Coriolis forces (alt0=0; lat0=0; k=1e-4)
  v0x = 30000.0;
  v0y = 0.0;
  v0z = 16000.0;

  */


  // run RK4() one more time but now save to file
  myfile.precision(10);  // Set precision for numeric output to myfile to 10 digits
  myfile << "# t, x, y, z, vx, vy, vz, alt, distance from start along great circle of Earth's radius  (units are MKS)" << endl;
  RK4(true,v0x,v0y,v0z);



  // project v0 onto rHat, thHat, and phiHat
  double v0r = sin(th0)*v0x + cos(th0)*v0z ;
  double v0th = cos(th0)*v0x - sin(th0)*v0z ;
  double v0phi = v0y;


  cout << "  Optimized:"
       << "    v0x = " << v0x
       << "    v0y = " << v0y
       << "    v0z = " << v0z  << endl
       << "    v0 = " << v0 << endl
       << "    elevation from ground (in degrees) = " << atan2(v0r, sqrt( v0th*v0th + v0phi*v0phi )) / piOver180 << endl
       << "    clockwise angle from North (in degrees) = " << atan2(v0phi, -v0th) / piOver180 << endl;

  cout << "  Target:"
       << "    alt = " << altEnd
       << "    lat = " << latEnd
       << "    long = " << longEnd  << endl;

  cout << "  Actual:"
       << "    alt = " << alt
       << "    lat = " << 90 - th / piOver180
       << "    long = " << phi / piOver180
       << "    t = " << t  << endl;

  cout << "  Straight-line miss distance = "
       << sqrt( pow(xEnd-x,2) + pow(yEnd-y,2) + pow(zEnd-z,2) )
       << endl;



  // Produce another file when changing v0 by 1% in many random ways (all use same altEnd)
  //cout << "  randomly varying v0..." << endl;
  myfile2.precision(8);
  myfile2 << "# v0x, v0y, v0z, miss x, miss y, miss z  (the miss vector projected onto phiHat,-thHat,rHat of target location so that y is North)" << endl;
  srand (static_cast <unsigned> (time(0)));  // "seed" random number generation
  varyV0(1000);  // do 1000 launches


}

 
