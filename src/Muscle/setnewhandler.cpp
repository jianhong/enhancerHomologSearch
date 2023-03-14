#include "muscle.h"
#include <stdio.h>
#include <new>
#include <Rcpp.h>

const int ONE_MB = 1024*1024;
const size_t RESERVE_BYTES = 8*ONE_MB;
static void *EmergencyReserve = 0;

void OnOutOfMemory()
	{
	free(EmergencyReserve);
	Rcpp::Rcerr << "\n*** OUT OF MEMORY ***\n";
	char msg[50];
	snprintf(msg, 50, "Memory allocated so far %g MB\n", GetMemUseMB());
	Rcpp::Rcerr << msg;
	extern MSA *ptrBestMSA;
	if (ptrBestMSA == 0)
	  Rcpp::Rcerr << "No alignment generated\n";
	else
		SaveCurrentAlignment();
	throw EXIT_FatalError;
	}

void SetNewHandler()
	{
	EmergencyReserve = malloc(RESERVE_BYTES);
	std::set_new_handler(OnOutOfMemory);
	}

void CleanupNewHandler() {
	free(EmergencyReserve);
}
