#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/*
  kiso.c

  C code for isotropic edge correction 

  Note: coordinates assumed to be normalised to unit radius
        angles are in radians
*/

double kisocapweight();
int IntersectCircles();

void kisocap(n, x1, x2, x3, Dmat, centre, height, wmat) 
     int *n;               /* number of points */
     double *x1, *x2, *x3; /* Cartesian coordinates of points */
     double *Dmat;         /* matrix of distances between each pair of points */
     double *centre;       /* Cartesian coords of centre of cap */
     double *height;       /* cap height */
     double *wmat;         /* output - (inverse) weight matrix */
{
  int N;
  double c1, c2, c3, h;
  int i,j,ijpos,jipos,ncut;
  double dij, cosdij;
  double cutA[3], cutB[3]; 

  N = *n;

  /* cap parameters */
  c1 = centre[0];
  c2 = centre[1];
  c3 = centre[2];
  h = *height;

  for(i = 1; i < N; i++) {
    for(j = 0; j < N; j++) {
      ijpos = N * j + i;
      if(i == j) {
	wmat[ijpos] = 1.0; 
      } else {
	dij = Dmat[ijpos];
	cosdij = cos(dij);
	ncut = IntersectCircles(c1, c2, c3, h, 
				x1[i], x2[i], x3[i], cosdij,
				1e-8, cutA, cutB);
	if(ncut == 2) {
	  wmat[ijpos] = kisocapweight(x1[i], x2[i], x3[i], 
				      x1[j], x2[j], x3[j],
				      cosdij,
				      c1, c2, c3, h,
				      cutA, cutB);
	} else wmat[ijpos] = 1.0;
      } 
    }
  }
}

double kisocapweight(xi1, xi2, xi3, 
		     xj1, xj2, xj3,
		     cosdij,
		     c1, c2, c3, h,
		     y, z) 
     double xi1, xi2, xi3;
     double xj1, xj2, xj3;
     double cosdij;
     double c1, c2, c3; /* centre of cap */
     double h; /* distance from sphere centre to plane of cap */
     double *y, *z;
{
  double v1, v2, v3, s1, s2, s3, vlen;
  double ydotv, zdotv, ydots, zdots, thetay, thetaz, sindij;
  double alpha, sinalpha, cosalpha, m1, m2, m3, mdotc;
  double theta[3];
  double anglefrac;
  /* orthonormalise */
  v1 = xj1 - cosdij * xi1;
  v2 = xj2 - cosdij * xi2;
  v3 = xj3 - cosdij * xi3;
  vlen = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
  v1 = v1/vlen;
  v2 = v2/vlen;
  v3 = v3/vlen;
  s1 = xi2 * v3 - xi3 * v2;
  s2 = -xi1 * v3 + xi3 * v1; 
  s3 = xi1 * v2 - xi2 * v1;
  /* represent intersections using orthonormal basis */
  ydotv = y[0] * v1 + y[1] * v2 + y[2] * v3;
  zdotv = z[0] * v1 + z[1] * v2 + z[2] * v3;
  ydots = y[0] * s1 + y[1] * s2 + y[2] * s3;
  zdots = z[0] * s1 + z[1] * s2 + z[2] * s3;
  /* compute angular positions of intersections */
  thetay = fmod(atan2(ydots, ydotv), M_2PI);
  thetaz = fmod(atan2(zdots, zdotv), M_2PI);
  anglefrac = abs(thetay - thetaz)/M_2PI;
  /* determine whether angle zero is inside cap */
  sindij = sqrt(1 - cosdij * cosdij);
  m1 = sindij * v1 + cosdij * xi1;
  m2 = sindij * v2 + cosdij * xi2;
  m3 = sindij * v3 + cosdij * xi3;
  mdotc = m1 * c1 + m2 * c2 + m3 * c3;
  if(mdotc > h) {
    /* angle zero is inside the cap */
    return(1 - anglefrac);
  } 
  /* angle zero is outside the cap */
  return(anglefrac);
}

int IntersectCircles(a1, a2, a3, ah, b1, b2, b3, bh, tol,
		     y, z)
     double a1, a2, a3, ah, b1, b2, b3, bh, tol; /* input */
     double *y, *z; /* output - vectors of length 3 */
{
  /* Find intersection of two small circles on the sphere */
  double x1, x2, x3, xlen2, adotb, ca, cb, s1, s2, s3, t2, t;
  /* dot product of normals */
  adotb = a1 * b1 + a2 * b2 + a3 * b3;
  /* cross product of normals */
  x1 = a2 * b3 - a3 * b2;
  x2 = -a1 * b3 + a3 * b1;
  x3 = a1 * b2 - a2 * b1;
  /* squared length of cross product */
  xlen2 = x1 * x1 + x2 * x2 + x3 * x3;
  if(xlen2 == 0) {   
    /* Normals are either equal or opposite.
       Intersection is either empty or infinite. */
    if((adotb > 0 && abs(ah - bh) < tol) ||
       (adotb < 0 && abs(ah + bh) < tol))
      return(3); /* intersection is infinite */
    return(0); /* intersection is empty */
  } else {
    /* Either 1 or 2 intersection points */
    ca = (ah - bh * adotb)/xlen2;
    cb = (bh - ah * adotb)/xlen2;
    s1 = ca * a1 + cb * b1;
    s2 = ca * a2 + cb * b2;
    s3 = ca * a3 + cb * b3;
    t2 = (1 - ca * ca - cb * cb - 2 * ca * cb * adotb)/xlen2;
    if(t2 >= tol) {
      /* general position - two intersection points */
      t = sqrt(t2);
      y[0] = s1 + t * x1;
      y[1] = s2 + t * x2;
      y[2] = s3 + t * x3;
      z[0] = s1 - t * x1;
      z[1] = s2 - t * x2;
      z[2] = s3 - t * x3;
      return(2);
    } else if(t2 >= 0.0) {
      /* osculating circles - one intersection point */
      y[0] = s1;
      y[1] = s2;
      y[2] = s3;
      return(1);
    } 
  }
  return(0);
}

