#include <Rcpp.h>
#include <getopt.h>
#include <string>
#include <vector>
#include <fstream>
#include <time.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "ushuffle.h"
#ifdef __cplusplus
}
#endif

#ifndef USHUFFLER_CONF
#define USHUFFLER_CONF
#define USHUFFLER_VERSION "0.0.1"
#endif

using namespace Rcpp;

char** ou_ushuffler(char* s, int l, int k, int n){
  shuffle1(s, l, k);
  char** t = new char*[n];
  for(int i = 0; i<n; i++){
    if ((t[i] = (char*) malloc(l + 1)) == NULL) {
      stop("malloc failed\n");
    }
    t[i][l] = '\0';
    shuffle2(t[i]);
  }
  return(t);
}

// [[Rcpp::export]]
CharacterVector rushuffle(CharacterVector x, IntegerVector k, IntegerVector n) {
  int N = as<int>(n);
  int L = x.length();
  int K = as<int>(k);
  CharacterVector seq;
  for(int i=0; i<L; i++){
    std::string ins = as<std::string>(x[i]);
    char** res = ou_ushuffler((char*)ins.c_str(), ins.size(), K, N);
    for(int j=0; j<N; j++){
      seq.push_back(std::string(res[j]));
    }
  }
  return(seq);
}

