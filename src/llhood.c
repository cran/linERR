#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP llhood(SEXP PAR, SEXP DATA, SEXP RSETS, SEXP AGES, SEXP DOSE, SEXP COVARS1, 
            SEXP N, SEXP M, SEXP P, SEXP FAILTIMES, SEXP NCOVS1, SEXP MAXEXP, SEXP LAG){
  double *res, *failtime, *lag, auxDose, auxDoseC, auxCov, auxCCov,
          sum, sumC, b, c, tmp;
  int i, j, k, l, nrisk, ncases, *maxexp = INTEGER(MAXEXP), *n = INTEGER(N), *m = INTEGER(M),
      *ncovs1 = INTEGER(NCOVS1);
  double dose[*n], doseC[*n], par[*ncovs1+1], sumacovs[*n], sumaCcovs[*n];
  SEXP RES;
  
  PROTECT(RES = allocVector(REALSXP, *m));
  res = REAL(RES);
  lag = REAL(LAG);
  for (l=0; l < *m; l++){
    failtime = REAL(VECTOR_ELT(FAILTIMES, l));
    nrisk  = 0;
    ncases = 0;
    for (i=0; i < *n; i++){
      k = 2;
      auxDose  = 0.0;
      auxDoseC = 0.0;
      auxCov   = 0.0;
      auxCCov  = 0.0;
      
      if (REAL(VECTOR_ELT(RSETS, l))[i] == 1.0){
        nrisk = nrisk + 1;
        if (*ncovs1!=0){
          for (j=0; j < *ncovs1; j++){
            par[j+1] = REAL(VECTOR_ELT(PAR, j+1))[0];
            auxCov += par[j+1]*REAL(VECTOR_ELT(COVARS1, j))[i];
          }
        }
        tmp = *failtime - *lag;
        while((REAL(VECTOR_ELT(AGES, k-1))[i] <= tmp) && (k < *maxexp + 1)){
          auxDose += REAL(VECTOR_ELT(DOSE, k-1))[i];
          if (REAL(VECTOR_ELT(DATA, 5))[i] == 1.0 && REAL(VECTOR_ELT(DATA, 4))[i] == *failtime){
            auxDoseC += REAL(VECTOR_ELT(DOSE, k-1))[i];
          }
          k++;
        }

        if (REAL(VECTOR_ELT(DATA, 5))[i] == 1.0 && REAL(VECTOR_ELT(DATA, 4))[i] == *failtime){
          ncases = ncases + 1;
          if (*ncovs1!=0){
            for (j=0; j < *ncovs1; j++){
              par[j+1] = REAL(VECTOR_ELT(PAR, j+1))[0];
               auxCCov += par[j+1]*REAL(VECTOR_ELT(COVARS1, j))[i];
            }
          }
        }
      }
      
      dose[i]      = auxDose;
      doseC[i]     = auxDoseC;
      sumacovs[i]  = auxCov;
      sumaCcovs[i] = auxCCov;
    }
    sum = 0.0; sumC = 0.0;
    // Log-linear term covariates
    par[0] = REAL(VECTOR_ELT(PAR, 0))[0];
    for (j = 0; j < *n; j++){
      b = 0.0; c = 0.0;
      if (*ncovs1 != 0){
        if(REAL(VECTOR_ELT(RSETS, l))[j] == 1.0){
          b = exp(sumacovs[j]);
          if (REAL(VECTOR_ELT(DATA, 5))[j] == 1.0 && REAL(VECTOR_ELT(DATA, 4))[j] == *failtime){
            c = exp(sumaCcovs[j]);
          }
        }
      }
      sum  += par[0]*exp(sumacovs[j])*dose[j]+b;
      sumC += par[0]*exp(sumaCcovs[j])*doseC[j]+c;
    }

    if (*ncovs1 == 0){
      sum = nrisk + sum;
      sumC = ncases + sumC;
    }
    res[l] = log(sumC/sum);
  }
  UNPROTECT(1);
  return RES;
}
