#include "muscle.h"
#include <stdio.h>
#include <Rcpp.h>

void Credits()
	{
	static bool Displayed = false;
	if (Displayed)
		return;

	Rcpp::Rcerr << std::endl << MUSCLE_LONG_VERSION << std::endl << std::endl;
	Rcpp::Rcerr << "http://www.drive5.com/muscle\n";
	Rcpp::Rcerr << "This software is donated to the public domain.\n";
	Rcpp::Rcerr << "Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.\n\n";
	Displayed = true;
	}

void Usage()
	{
	Credits();
  Rcpp::Rcerr <<
"\n"
"Basic usage\n"
"\n"
"    muscle -in <inputfile> -out <outputfile>\n"
"\n"
"Common options (for a complete list please see the User Guide):\n"
"\n"
"    -in <inputfile>    Input file in FASTA format (default stdin)\n"
"    -out <outputfile>  Output alignment in FASTA format (default stdout)\n"
"    -diags             Find diagonals (faster for similar sequences)\n"
"    -maxiters <n>      Maximum number of iterations (integer, default 16)\n"
"    -maxhours <h>      Maximum time to iterate in hours (default no limit)\n"
"    -html              Write output in HTML format (default FASTA)\n"
"    -msf               Write output in GCG MSF format (default FASTA)\n"
"    -clw               Write output in CLUSTALW format (default FASTA)\n"
"    -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header\n"
"    -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)\n"
"    -quiet             Do not write progress messages to stderr\n"
"    -version           Display version information and exit\n"
"\n"
"Without refinement (very fast, avg accuracy similar to T-Coffee): -maxiters 2\n"
"Fastest possible (amino acids): -maxiters 1 -diags -sv -distance1 kbit20_3\n"
"Fastest possible (nucleotides): -maxiters 1 -diags\n";
	}
