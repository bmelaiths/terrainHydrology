#include <stdio.h>

#include "hydrologyParameters.hpp"
#include "hydrologyFunctions.hpp"

/*
A list of data structures that will be used:
* A KDTree that will be used to query node locations
* A vector that will store actual node structs (this is where
  the node information will be, such as elevation, priority,
  etc). The tree will merely store an index in the vector
* Some data structure to represent the shore
* A graph
* A list of candidate nodes (this must be parallel)
*/

int main() {
  //gather inputs

  #define INPUT stdin
  // #define FILEINPUT
  #ifdef FILEINPUT
  FILE *input = fopen("./binaryFile", "rb");

  if (input == NULL)
  {
    printf("Unable to open file\n");
    exit(1);
  }
  #endif

  HydrologyParameters params = readParamsFromStream(INPUT);


  // perform computatons

  const uint8_t anotherNode = 0x2e, allDone = 0x21;

  while (params.candidates.size() > 0)
  {
    Primitive selectedCandidate = selectNode(params);
    alpha(selectedCandidate, params);
    fwrite(&anotherNode, sizeof(uint8_t), 1, stdout);
    fflush(stdout);
  }
  fwrite(&allDone, sizeof(uint8_t), 1, stdout);
  fflush(stdout);


  //export outputs

  writeBinary(params.hydrology, stdout);
  fflush(stdout);


  //free resources

  #ifdef FILEINPUT
  fclose(INPUT);
  #endif

  return 0;
}