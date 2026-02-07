#include "run.hpp"

// ==============================================================================
// Initializes and runs a simulation using the MFloe library.
//
// This function initializes a `simu` object from the `MFloe` class,
// checks for correct command-line arguments,
// loads either an integer value or directly uses a string filename as the
// configuration file, resets close lists related to collision detection,
// and then integrates the simulation.
// ==============================================================================
int main(int argc, char const *argv[]) {
  MFloe simu;
  simu.head();
  
  try {
    // Define the command line object
    TCLAP::CmdLine cmd("MapleFloe Simulation", ' ', "0.0");

    // Define a value argument and add it to the command line
    TCLAP::UnlabeledValueArg<std::string> inputConfArg("input-conf-file", "Input configuration file", true, "",
                                                       "string");
    cmd.add(inputConfArg);

    // Parse the command line arguments
    cmd.parse(argc, argv);

    // Get the value parsed by each argument
    std::string inputConf = inputConfArg.getValue();

    if (fileTool::containsOnlyDigits(inputConf.c_str())) {
      int num = std::atoi(inputConf.c_str());
      std::cout << MFLOE_INFO << "load conf" << num << std::endl;
      simu.loadConf(num);
    } else {
      std::cout << MFLOE_INFO << "load " << inputConf << std::endl;
      simu.loadConf(inputConf.c_str());
    }

    simu.updateNeighbors(simu.dVerlet);
    simu.integrate();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
