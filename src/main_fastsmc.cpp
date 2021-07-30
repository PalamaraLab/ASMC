#include <iostream>
#include <string>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "HMM.hpp"
#include "Timer.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  // Parse input arguments
  DecodingParams params;
  if (!params.processCommandLineArgsFastSMC(argc, argv)) {
    cerr << "Error processing command line; exiting." << endl;
    exit(1);
  }

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  return 0;
}
