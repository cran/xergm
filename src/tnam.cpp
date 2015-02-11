#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector netLagCppLoop(IntegerMatrix mat, IntegerMatrix pdistmat, 
    IntegerVector pathdist, IntegerVector decay, NumericVector y,
    std::string normalization, bool reciprocal) {
  int normal = 0;
  int recip = 0;
  NumericVector result(y.size());
  if (normalization == "no") {
    normal = 1;
  }
  for (int i = 0; i < mat.nrow(); i++) {
    if (normalization == "row") {
      for (int j = 0; j < mat.ncol(); j++) {
        if (!R_IsNA(mat(i, j))) {
          normal = normal + mat(i, j);
        }
      }
    }
    for (int j = 0; j < mat.ncol(); j++) {
      if (normalization == "column") {
        for (int k = 0; k < mat.nrow(); k++) {
          if (!R_IsNA(mat(k, j))) {
            normal = normal + mat(k, j);
          }
        }
      }
      if (normal == 0) {  // avoid division by zero
        normal = 1;
      }
      if ((reciprocal == true && !R_IsNA(mat(i, j)) && !R_IsNA(mat(j, i)) && 
          mat(i, j) == mat(j, i)) || reciprocal == false) {
        recip = 1;
      } else {
        recip = 0;
      }
      int decayindex = -1;
      for (int k = 0; k < pathdist.size(); k++) {
        if (!R_IsNA(pdistmat(i, j)) && pathdist[k] == pdistmat(i, j)) {
          decayindex = k;
        }
      }
      if (decayindex > -1) {
        if (R_IsNA(y[j])) {
          result[i] = NA_REAL;
        } else {
          result[i] = result[i] + ((y[j] / normal) * decay[decayindex] * recip);
        }
      }
    }
  }
  return result;
}
