#ifndef _RMuscle_R_MUSCLE_H
#define _RMuscle_R_MUSCLE_H

#include "Muscle/muscle.h"
#include "Muscle/seqvect.h"

#include <Rcpp.h>

// [[Rcpp::export]]
SEXP RMuscle(SEXP rInputSeq,
             SEXP rCluster,
             SEXP rGapOpening,
             SEXP rGapExtension,
             SEXP rMaxiters,
             SEXP rSubstitutionMatrix,
             SEXP rType,
             SEXP rVerbose,
             SEXP rParams);

struct MuscleInput {
    SeqVect inputSeqs;
    std::vector<std::string> seqNames;
    std::vector<std::string> colNames;
    bool hasSubstitutionMatrix;
    float substitutionMatrix[32][32];
};

struct MuscleOutput { //can be used for further result objects
    std::vector<std::string> msa; //multiple sequence alignment
};

void DoMuscle(MuscleInput *msaInput, MuscleOutput *msaOutput);
void Run(MuscleInput *msaInput, MuscleOutput *msaOutput);

#endif
