#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void dkisoband(int *, double *, double *, double *, double *, double *, double *, int *, int *, double *, double *);
extern void dkisocap(int *, double *, double *, double *, double *, double *, int *, double *, double *);
extern void dwkisoband(int *, double *, double *, double *, double *, double *, double *, double *, int *, int *, double *, double *);
extern void dwkisocap(int *, double *, double *, double *, double *, double *, double *, int *, double *, double *);
extern void kisobandweights(int *, double *, double *, double *, double *, double *, double *, double *, int *, double *);
extern void kisocapweights(int *, double *, double *, double *, double *, double *, double *, double *);
extern void RcallPtInSphPoly(double *, double *, int *, double *, double *, int *, double *, double *, int *);
extern void Rpolyarea(double *, double *, int *, double *, double*);

static const R_CMethodDef CEntries[] = {
    {"dkisoband",        (DL_FUNC) &dkisoband,        11},
    {"dkisocap",         (DL_FUNC) &dkisocap,          9},
    {"dwkisoband",       (DL_FUNC) &dwkisoband,       12},
    {"dwkisocap",        (DL_FUNC) &dwkisocap,        10},
    {"kisobandweights",  (DL_FUNC) &kisobandweights,  10},
    {"kisocapweights",   (DL_FUNC) &kisocapweights,    8},
    {"RcallPtInSphPoly", (DL_FUNC) &RcallPtInSphPoly,  9},
    {"Rpolyarea",        (DL_FUNC) &Rpolyarea,         5},
    {NULL, NULL, 0}
};

void R_init_spherstat(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
