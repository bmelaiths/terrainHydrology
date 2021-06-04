#include "ts.hpp"

Terrain::~Terrain()
{
  for (T *t : allTs)
  {
    delete t;
  }
  
}

void Terrain::dumpT(Point loc, size_t cellID)
{
  T *t = new T(loc, cellID);
  allTs.push_back(t);
}

size_t Terrain::numTs() const
{
  return allTs.size();
}

T& Terrain::getT(size_t id)
{
  return *allTs[id];
}