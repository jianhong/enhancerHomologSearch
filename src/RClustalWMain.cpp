#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream>
#include "ClustalW/alignment/Alignment.h"
#include "ClustalW/alignment/Sequence.h"
#include "ClustalW/general/clustalw.h"
#include "ClustalW/general/UserParameters.h"
#include "ClustalW/substitutionMatrix/SubMatrix.h"
#include "ClustalW/general/Utility.h"
#include "ClustalW/fileInput/FileReader.h"
#include "ClustalW/interface/CommandLineParser.h"
#include "ClustalW/general/DebugLog.h"
#include "ClustalW/general/ClustalWResources.h"
#include "ClustalW/general/Stats.h"
#include <ctime>

#include <vector>
#include "RClustalWMain.h"
using namespace std;
using namespace clustalw;
using namespace Rcpp;

namespace clustalw
{
    UserParameters* userParameters;
    SubMatrix *subMatrix;
    Utility* utilityObject;
    DebugLog* logObject;
    Stats* statsObject;

RClustalWMain::RClustalWMain() {}

RClustalWMain::~RClustalWMain() {}

void RClustalWMain::run(UserArgs args, ClustalWInput *input, ClustalWOutput *output) {

	userParameters = new UserParameters(false);
	utilityObject = new Utility();
	subMatrix = new SubMatrix();
	statsObject = new Stats();
	ClustalWResources *resources = ClustalWResources::Instance();
	resources->setPathToExecutable(args.at(0));
	userParameters->setDisplayInfo(true);



	//userParameters->setDebug(5);
	#if DEBUGFULL
		if(DEBUGLOG) {
			Rcpp::Rcout << "debugging is on\n\n\n";
			logObject = new DebugLog("logfile.txt");
			logObject->logMsg("Loggin is on!");
		}
	#endif

	if (args.size() > 1) {
		//time_t start, end;
		//double dif;
		//start = time (NULL);
		//userParameters->setDisplayInfo(false);
		CommandLineParser cmdLineParser(&args, false);
		cmdLineParser.run(&args, false, input, output);
		//for (int i = 0; i < result.size(); i++) {
		//	utilityObject->info("[%s]", result[i].c_str());
		//}

		//if (statsObject->isEnabled())
		//	statsObject->logCmdLine(argc,argv);

		//end = time (NULL);
		//dif = difftime(end, start);
		//Rcpp::Rcout << "It took " << dif << " seconds\n";
	}

	delete userParameters;
	delete utilityObject;
	delete subMatrix;
	delete statsObject;

	if(logObject)
	{
		delete logObject;
	}
	return;
}
}
