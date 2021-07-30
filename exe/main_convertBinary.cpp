#include "BinaryDataReader.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{

  // make sure parameters are ok
  if (argc != 2) {
    std::cout << "Number of parameters is wrong." << std::endl;
    std::cout << "Only one parameter (name of binary file) is required." << std::endl;
    exit(1);
  }

  BinaryDataReader binaryDataReader(argv[1]);

  while (binaryDataReader.moreLinesInFile()) {
    std::cout << binaryDataReader.getNextLine().toString() << std::endl;
  }

  return 0;
}
