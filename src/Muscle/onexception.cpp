#include "muscle.h"
#include <stdio.h>
#include <Rcpp.h>

static char szOnExceptionMessage[] =
	{
	"\nFatal error, exception caught.\n"
	};

void OnException()
	{
	Rcpp::Rcerr << szOnExceptionMessage;
	Log("%s", szOnExceptionMessage);
	Log("Finished %s\n", GetTimeAsStr());
		throw EXIT_Except;
	}
