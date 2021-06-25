#include <stdio.h>

#include <omp.h>

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
  std::vector<GEOSContextHandle_t> geosContexts;
  for (int i = 0; i < omp_get_max_threads(); i++)
  {
    geosContexts.push_back(GEOS_init_r());
  }


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

  PrimitiveParameters params(INPUT, geosContexts[0]);


  //perform computations
  const uint8_t anotherNode = 0x2e, allDone = 0x21;
  #pragma omp parallel for
  for (size_t i = 0; i < params.ts.numTs(); i++)
  {
    T& t = params.ts.getT(i);
    t.setElevation(computePrimitiveElevation(
      t, params.hydrology, params.cells, params.ts, params.contour,
      params.resolution, geosContexts[omp_get_thread_num()]
    ));
    fwrite(&anotherNode, sizeof(uint8_t), 1, stdout);
    fflush(stdout);
  }
  fwrite(&allDone, sizeof(uint8_t), 1, stdout);
  fflush(stdout);

  // write results
  for (size_t i = 0; i < params.ts.numTs(); i++)
  {
    float elev = params.ts.getT(i).getElevation();
    elev = float_tobe(elev);
    fwrite(&elev, sizeof(float), 1, stdout);
    fflush(stdout);
  }

  for (GEOSContextHandle_t geosContext : geosContexts)
  {
    GEOS_finish_r(geosContext);
  }
  
  #ifdef FILEINPUT
  fclose(input);
  #endif
}