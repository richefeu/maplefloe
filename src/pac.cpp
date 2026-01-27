#include "pac.hpp"

int main(int argc, char const *argv[]) {  
  
  /*
  if (argc != 2) {
    std::cout << "usage: ./pac <pac-command.txt>" << std::endl;
    return 0;
  }
  
  MFloeAssembler Assemb;
  Assemb.read(argv[1]);
  Assemb.run(); 
  */
  
  
  try {
    // Define the command line object
    TCLAP::CmdLine cmd("MappleFloe Packing", ' ', "0.0");

    // Define a value argument and add it to the command line
    TCLAP::UnlabeledValueArg<std::string> inputFileArg("input-file", "Input file", true, "",
                                                       "string");
    cmd.add(inputFileArg);

    // Parse the command line arguments
    cmd.parse(argc, argv);

    // Get the value parsed by each argument
    std::string inputFile = inputFileArg.getValue();
    
    MFloeAssembler Assemb;
    Assemb.read(inputFile.c_str());
    Assemb.run();
    //Assemb.save();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}