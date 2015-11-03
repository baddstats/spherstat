#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/*
  kiso.c

  C code for isotropic edge correction 

  Note: coordinates assumed to be normalised to unit radius
        angles are in radians
*/

#define JA (0 == 0)
#define NEIN (!(JA))

/* declarations of helper functions */
double isoWcap(), isoWband();
int IntersectCircles();

/* ----------------- R-callable interface functions --------------------- */

void kisocapweights(n, x1, x2, x3, Dmat, centre, height, wmat) 
     int *n;               /* number of points */
     double *x1, *x2, *x3; /* Cartesian coordinates of points */
     double *Dmat;         /* matrix of distances between each pair of points */
     double *centre;       /* Cartesian coords of centre of cap */
     double *height;       /* height of plane defining cap */
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
	  wmat[ijpos] = isoWcap(x1[i], x2[i], x3[i], 
				x1[j], x2[j], x3[j],
				cosdij,
				c1, c2, c3, h,
				cutA, cutB);
	} else wmat[ijpos] = 1.0;
      } 
    }
  }
}

void kisobandweights(n, x1, x2, x3, Dmat, centre, height1, height2, iscomp, 
		     wmat) 
     int *n;               /* number of points */
     double *x1, *x2, *x3; /* Cartesian coordinates of points */
     double *Dmat;         /* matrix of distances between each pair of points */
     double *centre;       /* Cartesian coords of centre of band */
     double *height1,
            *height2;      /* heights of planes defining band */
     int *iscomp;         /* 0 if window is a band, 1 if it's the complement */
     double *wmat;         /* output - (inverse) weight matrix */
{
  int N;
  double c1, c2, c3, hup, hlo;
  int i,j,ijpos,jipos, ncutup, ncutlo, ncuts, freepos, iscomplement;
  double dij, cosdij;
  double cutAup[3], cutBup[3], cutAlo[3], cutBlo[3]; 
  double cuts[12];  /* 3 x 4 matrix */

  N = *n;
  iscomplement = *iscomp;

  /* cap parameters */
  c1 = centre[0];
  c2 = centre[1];
  c3 = centre[2];
  hup = *height1;
  hlo = *height2;

  for(i = 1; i < N; i++) {
    for(j = 0; j < N; j++) {
      ijpos = N * j + i;
      wmat[ijpos] = 1.0; 
      if(i != j) {
	dij = Dmat[ijpos];
	cosdij = cos(dij);
	/* find all crossing points */
	ncutup = IntersectCircles(c1, c2, c3, hup, 
				  x1[i], x2[i], x3[i], cosdij,
				  1e-8, cutAup, cutBup);
	ncutlo = IntersectCircles(c1, c2, c3, hlo, 
				  x1[i], x2[i], x3[i], cosdij,
				  1e-8, cutAlo, cutBlo);
	/* Treat infinite intersection as empty */
	if(ncutup > 2) ncutup == 0;
	if(ncutlo > 2) ncutlo == 0;
	ncuts = ncutup + ncutlo;
	if(ncuts > 0) {
	  /* Pack coordinates of intersection points into array */
	  freepos = 0; 
	  if(ncutup >= 1) {
	    cuts[0] = cutAup[0];
	    cuts[1] = cutAup[1];
	    cuts[2] = cutAup[2];
	    if(ncutup == 2) {
	      cuts[3] = cutBup[0];
	      cuts[4] = cutBup[1];
	      cuts[5] = cutBup[2];
	    }
	    freepos = 3 * ncutup;
	  }
	  if(ncutlo >= 1) {
	    cuts[freepos]     = cutAlo[0];
	    cuts[freepos + 1] = cutAlo[1];
	    cuts[freepos + 2] = cutAlo[2];
	    if(ncutlo == 2) {
	      cuts[freepos + 3] = cutBlo[0];
	      cuts[freepos + 4] = cutBlo[1];
	      cuts[freepos + 5] = cutBlo[2];
	    }
	  }
	  wmat[ijpos] = isoWband(x1[i], x2[i], x3[i], 
				 x1[j], x2[j], x3[j],
				 cosdij,
				 c1, c2, c3, hup, hlo, iscomplement, 
				 cuts, ncuts);
	}
      }
    }
  }
}

/* ------------- internal functions --------------------------- */

double isoWcap(xi1, xi2, xi3, 
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
  double m1, m2, m3, mdotc;
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

double isoWband(xi1, xi2, xi3, 
		xj1, xj2, xj3,
		cosdij,
		c1, c2, c3, hup, hlo, iscomp,
		cuts, ncuts) 
     double xi1, xi2, xi3;
     double xj1, xj2, xj3;
     double cosdij;
     double c1, c2, c3; /* centre of band */
     double hup, hlo; /* distances from sphere centre to planes of caps */
     int iscomp; /* 0 if window is a band, 1 if it's the complement */
     double *cuts;
     int ncuts;
{
  double v1, v2, v3, s1, s2, s3, vlen;
  double sindij, cuti1, cuti2, cuti3, cutdotV, cutdotS;
  double alpha, sinalpha, cosalpha, m1, m2, m3, mdotc;
  double theta[5];
  double totlen, tmp;
  int i, ipos, sorted;
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
  /* represent intersections using orthonormal basis; 
     compute angular positions */
  for(i = 0; i < ncuts; i++) {
    ipos = 3 * i;
    cuti1 = cuts[ipos];
    cuti2 = cuts[ipos+1];
    cuti3 = cuts[ipos+2];
    cutdotV = cuti1 * v1 + cuti2 * v2 + cuti3 * v3;
    cutdotS = cuti1 * s1 + cuti2 * s2 + cuti3 * s3;
    theta[i] = fmod(atan2(cutdotS, cutdotV), M_2PI);
  }
  /* sort angles */
  do {
    sorted = JA;
    for(i = 0; i < ncuts-1; i++) {
      if(theta[i] > theta[i+1]) {
	sorted = NEIN;
	tmp = theta[i];
	theta[i] = theta[i+1];
	theta[i+1] = tmp;
      }
    }
  } while(!sorted);
  /* add additional angle */
  theta[ncuts] = theta[0] + M_2PI;
  /* start measuring intervals */
  totlen = 0.0;
  sindij = sqrt(1 - cosdij * cosdij);
  /* examine intervals */
  for(i = 0; i < ncuts; i++) {
    /* midpoint of interval */
    alpha = (theta[i] + theta[i+1])/2.0;
    cosalpha = cos(alpha);
    sinalpha = sin(alpha);
    m1 = sindij * (v1 * cosalpha + s1 * sinalpha) + cosdij * xi1;
    m2 = sindij * (v2 * cosalpha + s2 * sinalpha) + cosdij * xi2;
    m3 = sindij * (v3 * cosalpha + s3 * sinalpha) + cosdij * xi3;
    /* determine whether midpoint is inside window */
    mdotc = m1 * c1 + m2 * c2 + m3 * c3;
    if((iscomp == 0 && mdotc >= hlo && mdotc <= hup) ||
       (iscomp != 0 && (mdotc <= hlo || mdotc >= hup))) {
	/* midpoint is inside window; add length of corresponding interval */
      totlen += theta[i+1] - theta[i];
    }
  } 
  return(totlen/M_2PI);
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

