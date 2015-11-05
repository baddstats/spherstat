#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/* 
   Given some spherical polygon S and some point X known to be located
   inside S, determine whether an arbitrary point P lies inside S, 
   outside S, or on its boundary.

   REFERENCE: Bevis, M. and Chatelain, J.-L. (1989) 
              Mathematical Geology, vol 21.

   Translated to C by Adrian Baddeley.

   UNITS AND SIGN CONVENTION:
   Latitudes and longitudes are specified in RADIANS.
   Latitudes are positive to the north and negative to the south.
   Longitudes are positive to the east and negative to the west.

   X must not lie on any great circle that includes two vertices of S.

   VERTEX ENUMERATION:
   The vertices of S should be numbered sequentially around the border 
   of the spherical polygon without repetition. Vertex 1 lies between
   vertex nv and vertex 2. Neighbouring vertices must be separated by
   less than pi radians (180 degrees). (In order to generate a polygon 
   side whose arc length equals or exceeds 180 degrees simply introduce 
   an additional (pseudo)vertex). Having chosen vertex 1, the user may
   number the remaining vertices in either direction. However if the user
   wishes to use the suhroutine SPA to determine the area of the polygon S
   (Bevis & Cambareri, 1987, Math. Geol., v.19, p. 335-346) then he or she 
   must follow the convention whereby in moving around the polygon border 
   in the direction of increasing vertex number clockwise bends occur at
   salient vertices. A vertex is salient if the interior angle is less
   than 180 degrees. (In the case of a convex polygon this convention 
   implies that vertices are numbered in clockwise sequence).

   RESULT:
   plocation = 0 implies P is outside of S
   plocation=1 implies P is inside of S
   plocation=2 implies P on boundary of S
   plocation=3 implies user error (P is antipodal to X)

*/

/* R interface function */
void RcallPtInSphPoly(plat, plon, np, vlat, vlon, nv, xlat, xlon, plocation) 
  double *plat, *plon; /* coordinates of query points */
  int *np;             /* number of query points */
  double *vlat, *vlon; /* coordinates of polygon vertices */
  int *nv;             /* number of vertices */
  double *xlat, *xlon; /* coordinates of point X known to be inside S */
  int *plocation;      /* output vector */
{
  double *scratch;
  void PtInSphPoly();
  scratch = (double *) R_alloc(*nv, sizeof(double));
  PtInSphPoly(plat, plon, *np, vlat, vlon, *nv, *xlat, *xlon, 
	      scratch, plocation);
}

/* main function */
void PtInSphPoly(plat, plon, np, vlat, vlon, nv, xlat, xlon, 
		 scratch, plocation) 
  double *plat, *plon; /* coordinates of query points */
  int np;              /* number of query points */
  double *vlat, *vlon; /* coordinates of polygon vertices */
  int nv;              /* number of vertices */
  double xlat, xlon;   /* coordinates of point X known to be inside S */
  double *scratch;     /* intermediate calculation storage */
  int *plocation;      /* output vector */
{
  double platj, plonj;
  double dellon, vAlat,vAlon,vBlat,vBlon,tlonA,tlonB,tlonP; 
  double tlon_X,tlon_P,tlon_B;
  double *tlonv;
  int i, iprev, inext, j, icross, istrike, location;
  int ibrngAB,ibrngAP,ibrngPB, ibrng_BX, ibrng_BP;

  double TrnsfmLon();
  int EastOrWest();

  tlonv = scratch;

  for(i = 0; i < nv; i++) {
    tlonv[i] = TrnsfmLon(xlat, xlon, vlat[i], vlon[i]);
    iprev = ((i > 0) ? i : nv) - 1;
    if(vlat[i] == vlat[iprev] && vlon[i] == vlon[iprev])
      error("Vertices %d and %d are identical", iprev, i);
    if(tlonv[i] == tlonv[iprev])
      error("Vertices %d and %d lie on the same great circle as X", iprev, i);
    if(vlat[i] == -vlat[iprev]) {
      dellon = vlon[i]-vlon[iprev];
      if(dellon > M_PI) dellon -= M_2PI;
      if(dellon < -M_PI) dellon += M_2PI;
      if(abs(dellon) == M_PI)
	error("Vertices %d and %d are antipodal", iprev, i);
    }
  }

  for(j = 0; j < np; j++) {
    platj = plat[j];
    plonj = plon[j];
    if(platj == -xlat) {
      dellon = plonj - xlon;
      if(dellon > M_PI) dellon -= M_2PI;
      if(dellon < -M_PI) dellon += M_2PI;
      if(abs(dellon) == M_PI) {
	warning("A data point is antipodal to X");
	plocation[j] = 3;
	continue;
      }
    }

    if(platj == xlat && plonj == xlon) {
      /* X and P identical, hence P is inside S */
      plocation[j] = 1; 
      continue;
    }

    location = 0; /* default (P is outside S) */
    icross   = 0; /* initialize counter */
    
    tlonP = TrnsfmLon(xlat, xlon, platj, plonj);


    /* ......... loop over sides of S ................ */
    for(i = 0; i < nv; i++) {
      vAlat = vlat[i];
      vAlon = vlon[i];
      tlonA = tlonv[i];

      inext = i + 1;
      if(inext == nv) inext = 0;
      vBlat = vlat[inext];
      vBlon = vlon[inext];
      tlonB = tlonv[inext];

      istrike = 0;

      if(tlonP == tlonA) { 
	istrike = 1;
      } else {
	ibrngAB = EastOrWest(tlonA,tlonB);
        ibrngAP = EastOrWest(tlonA,tlonP);
        ibrngPB = EastOrWest(tlonP,tlonB);
	if(ibrngAP == ibrngAB && ibrngPB == ibrngAB) { istrike=1 ; }
      }

      if(istrike == 1) {

	if(platj == vAlat && plonj == vAlon) {
	  plocation[j] = 2; /* P lies on a vertex of S */
	  break;
	}

	tlon_X = TrnsfmLon(vAlat,vAlon,xlat,xlon);
        tlon_B = TrnsfmLon(vAlat,vAlon,vBlat,vBlon);
        tlon_P = TrnsfmLon(vAlat,vAlon,platj,plonj);

	if(tlon_P == tlon_B) {
	  plocation[j] = 2; /* P lies on side of S */
	  break;
	} else {
	  ibrng_BX = EastOrWest(tlon_B,tlon_X);
          ibrng_BP = EastOrWest(tlon_B,tlon_P);
	  if(ibrng_BX == -ibrng_BP) ++icross;
	}
	
      }

    }
    /* end of loop over edges of S */

    /* 
       If the arc XP crosses the boundary S an even number of times 
       then P is in S
    */

    if(icross % 2 == 0) {
      plocation[j]=1;
    }
  }
}

/*
  This subroutine is required by subroutines DefSPolyBndry & LctPtRelBndry.
  It finds the 'longitude' of point Q in a geographic coordinate system
  for which point P acts as a 'north pole'. 
 
*/				
double TrnsfmLon(plat,plon,qlat,qlon)
				double plat, plon, qlat, qlon;
      {
	double a, t, b, tranlon;
	a = qlon - plon;
	t = sin(a) * cos(qlat);
	b = sin(qlat) * cos(plat) - cos(qlat) * sin(plat) * cos(a);
	tranlon = atan2(t,b);
	return(tranlon);
      }

/*
  This subroutine is required by subroutine LctPtRelBndry.
  This routine determines if in travelling the shortest path
  from point C (at longitude clon) to point Dlon (at longitude dlon)
  one is heading east, west or neither.
*/

int EastOrWest(clon, dlon)
     double clon, dlon;
{
  double del;
  int ibrng;
  del = dlon - clon;
  if(del > M_PI) del -= M_2PI;
  if(del < -M_PI) del += M_2PI;
  if(abs(del) == M_PI) return(0);
  return((del > 0.0) ? -1 : 1);
}

