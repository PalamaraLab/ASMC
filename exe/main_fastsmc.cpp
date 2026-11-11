#include <exception>
#include <iostream>
#include <string>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "HMM.hpp"
#include "Timer.hpp"

using namespace std;

namespace
{

int runFastSmc(int argc, char* argv[])
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

} // namespace

// Library code reports invalid input by throwing; report it here rather than letting std::terminate abort.
int main(int argc, char* argv[])
{
  try {
    return runFastSmc(argc, argv);
  } catch (const std::exception& e) {
    cerr << e.what() << '\n';
    return 1;
  }
}
