#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/* 
   Compute the area of a spherical polygon.

   Original (Fortran) code:
   Bevis and Cambareri (1987)
   Mathematical Geology, vol.19

   Translated/adapted to C by Adrian Baddeley

   Computes the area of a spherical polygon with nv vertices and sides.
   All arguments are double precision variables or arrays.

   ARGUMENTS:

   vlat, vlon      vectors containing the latitude and longitude
                   of each vertex. The ith. vertex is located at
                   [vlat(i),vlon(i)].  See notes below.

   nv              The number of vertices and sides in the spherical polygon

   rad             The radius of the sphere

   UNITS:
   Latitudes and longitudes are specified in RADIANS. The user selects
   the units of length in which to specify the radius, and the polygon
   area will be returned in the square of these units.

   SIGN CONVENTION:
   Latitudes are positive to the north and negative to the south.
   Longitudes are positive to the east and negative to the west.

   VERTEX ENUMERATION:
   The vertices are numbered sequentially around the border of the
   spherical polygon. Vertex 1 lies between vertex nv and vertex 2.
   The user must follow the convention whereby in moving around the
   polygon border in the direction of increasing vertex number clockwise
   bends occur at salient vertices. A vertex is salient if the interior
   angle is less than 180 degrees. (In the case of a convex polygon
   this convention implies that the vertices are numbered in clockwise
   sequence).
   Two adjacent vertices may never be exactly 180 degrees apart
   because they could be connected by infinitely many different
   great circle arcs, and thus the border of the spherical
   polygon would not be uniquely defined.
*/

void Rpolyarea(vlat,vlon,nv,rad,area)
     double *vlat, *vlon;
     int *nv;
     double *rad, *area;
{
  double polyarea();
  *area = polyarea(vlat, vlon, *nv, *rad);
}

double polyarea(vlat,vlon,nv,rad)
     double *vlat, *vlon;
     int nv;
     double rad;
{
  double area, sum, flat, flon, blat, blon, fang, bang, fvb;
  int iv, ivnext, ivprev;
  double trnsfrmlon();

  sum=0.0;

  for(iv = 0; iv < nv; iv++) {
    ivnext = iv+1;
    if(ivnext >= nv) ivnext = 0;
    ivprev = ((iv > 0) ? iv : nv) - 1;
    flat = vlat[ivnext];
    flon = vlon[ivnext];
    blat = vlat[ivprev];
    blon = vlon[ivprev];
    fang = trnsfrmlon(vlat[iv],vlon[iv],flat,flon);
    bang = trnsfrmlon(vlat[iv],vlon[iv],blat,blon);
    fvb=bang-fang;
    if(fvb < 0.0) fvb += M_2PI;
    sum += fvb;
  }
  area = (sum-M_PI*((double)(nv-2)))*rad*rad;
  return(area);
}
/*
  Finds the "longitude" of point Q in a geographic coordinate system
  for which point P acts as a "north pole'. 
*/

double trnsfrmlon(plat,plon,qlat,qlon)
     double plat, plon, qlat, qlon;
{
  double tranlon, t, b, dellon;
  dellon = qlon - plon;
  t = sin(dellon) * cos(qlat);
  b = sin(qlat) * cos(plat) - cos(qlat) * sin(plat)* cos(dellon);
  tranlon = atan2(t,b);
  return(tranlon);
}
