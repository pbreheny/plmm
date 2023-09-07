#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>

extern SEXP mfdr_gaussian(SEXP);

// List accessor function
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
    return elmt;
}


// Cross product of y with jth column of X
double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

static const R_CallMethodDef CallEntries[] = {
  {"mfdr_gaussian",  (DL_FUNC) &mfdr_gaussian,   1},
  {NULL, NULL, 0}
};

void R_init_penalizedLMM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
