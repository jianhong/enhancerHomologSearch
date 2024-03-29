#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include <time.h>
#include <Rcpp.h>

MSA *ptrBestMSA;
static const char *pstrOutputFileName;

void SetOutputFileName(const char *out)
	{
	pstrOutputFileName = out;
	}

void SetCurrentAlignment(MSA &msa)
	{
	ptrBestMSA = &msa;
	}

void SaveCurrentAlignment()
	{
	static bool bCalled = false;
	if (bCalled)
		{
		Rcpp::Rcerr <<
		  "\nRecursive call to SaveCurrentAlignment, giving up attempt to save.\n";
		throw EXIT_FatalError;
		}

	if (0 == ptrBestMSA)
		{
	  Rcpp::Rcerr << "\nAlignment not completed, cannot save.\n";
		Log("Alignment not completed, cannot save.\n");
		throw EXIT_FatalError;
		}

	if (0 == pstrOutputFileName)
		{
	  Rcpp::Rcerr << "\nOutput file name not specified, cannot save.\n";
		throw EXIT_FatalError;
		}

	Rcpp::Rcerr << "\nSaving current alignment ...\n";

	TextFile fileOut(pstrOutputFileName, true);
	ptrBestMSA->ToFASTAFile(fileOut);

	Rcpp::Rcerr << "Current alignment saved to \""<< pstrOutputFileName <<"\".\n";
	Log("Current alignment saved to \"%s\".\n", pstrOutputFileName);
	}

void CheckMaxTime()
	{
	if (0 == g_ulMaxSecs)
		return;

	time_t Now = time(0);
	time_t ElapsedSecs = Now - GetStartTime();
	if (ElapsedSecs <= (time_t) g_ulMaxSecs)
		return;

	Log("Max time %s exceeded, elapsed seconds = %ul\n",
	  MaxSecsToStr(), ElapsedSecs);

	SaveCurrentAlignment();
	throw EXIT_Success;
	}
