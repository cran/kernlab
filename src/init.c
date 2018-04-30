#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP fullsubstringk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smo_optim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stringtv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP subsequencek(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP substringk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tron_optim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"fullsubstringk", (DL_FUNC) &fullsubstringk,  6},
    {"smo_optim",      (DL_FUNC) &smo_optim,      23},
    {"stringtv",       (DL_FUNC) &stringtv,        7},
    {"subsequencek",   (DL_FUNC) &subsequencek,    6},
    {"substringk",     (DL_FUNC) &substringk,      6},
    {"tron_optim",     (DL_FUNC) &tron_optim,     27},
    {NULL, NULL, 0}
};

void R_init_kernlab(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
