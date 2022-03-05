#include <stdio.h>

#include <omp.h>

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
  #define INPUT input
  #define FILEINPUT
  #ifdef FILEINPUT
  FILE *input = fopen("./src/native-module/bin/binaryFile", "rb");

  if (input == NULL)
  {
    printf("Unable to open file\n");
    exit(1);
  }
  #endif

  HydrologyParameters params(INPUT);


  // perform computatons
  const uint8_t anotherNode = 0x2e, allDone = 0x21;
  #pragma omp parallel
  {
  // printf("Thread ID: %d\n", omp_get_thread_num());
  while (params.candidates.size() > 0)
  {
    params.lockCandidateVector();
    Primitive selectedCandidate = selectNode(params);
    params.unlockCandidateVector();

    alpha(selectedCandidate, params);

    #pragma omp critical
    {
    // write a byte to the calling program, signalling
    // that a candidate has been processed
    fwrite(&anotherNode, sizeof(uint8_t), 1, stdout);
    fflush(stdout);
    }
  }
  }
  // signall to the calling program that processing
  // is complete
  fwrite(&allDone, sizeof(uint8_t), 1, stdout);
  fflush(stdout);


  //export outputs
  params.hydrology.writeBinary(stdout);
  fflush(stdout);


  //free resources
  #ifdef FILEINPUT
  fclose(INPUT);
  #endif

  return 0;
}