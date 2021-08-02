#ifndef RCLUSTALW_H
#define RCLUSTALW_H

#include <Rcpp.h>

 // [[Rcpp::export]]
SEXP RClustalW(SEXP rInputSeqs,
               SEXP rCluster,
               SEXP rGapOpen,
               SEXP rGapExtend,
               SEXP rMaxIters,
               SEXP rsubstitutionMatrix,
               SEXP rType,
               SEXP rVerbose,
               SEXP rParams);

struct ClustalWInput {
    std::vector<std::string> inputSeqs;
    std::vector<std::string> seqNames;
    Rcpp::NumericMatrix substitutionMatrix;
};

struct ClustalWOutput { //can be use for further result objects
    std::vector<std::string> msa; //multiple sequence alignment
};

#endif
