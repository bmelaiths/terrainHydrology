#include <stdio.h>

#include "terrainPrimitives.hpp"
#include "terrainElevation.hpp"
#include "floatEndian.hpp"

void writeError(const char *fmt, ...)
{
  printf(fmt);
}

int main()
{
  //initialize the GEOS library
  initGEOS(NULL, &writeError);


  //gather inputs
  #define INPUT input
  #define FILEINPUT
  #ifdef FILEINPUT
  FILE *input = fopen("src/binaryFile", "rb");

  if (input == NULL)
  {
    printf("Unable to open file\n");
    exit(1);
  }
  #endif

  PrimitiveParameters params(input);


  //perform computations
  const uint8_t anotherNode = 0x2e, allDone = 0x21;
  for (size_t i = 0; i < params.ts.numTs(); i++)
  {
    T& t = params.ts.getT(i);
    t.setElevation(computePrimitiveElevation(
      t, params.hydrology, params.cells, params.ts,
      params.contour, params.resolution
    ));
    fwrite(&anotherNode, sizeof(uint8_t), 1, stdout);
    fflush(stdout);
  }
  fwrite(&allDone, sizeof(uint8_t), 1, stdout);
  fflush(stdout);

  //write results
  for (size_t i = 0; i < params.ts.numTs(); i++)
  {
    float elev = params.ts.getT(i).getElevation();
    elev = float_tobe(elev);
    fwrite(&elev, sizeof(float), 1, stdout);
    fflush(stdout);
  }
}