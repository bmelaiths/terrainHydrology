#include "terrainHoneycomb.hpp"

TerrainHoneycomb::~TerrainHoneycomb()
{
  for (Q *q : allQs)
  {
    if (q != NULL)
    {
      delete q;
    }
  }
}

void TerrainHoneycomb::dumpQ
(
  Point position, float elevation, size_t vorIndex,
  std::vector<size_t> neighbors
)
{
  Q *q = new Q(position, elevation, vorIndex, neighbors);

  allQs.push_back(q);
}

void TerrainHoneycomb::dumpNull()
{
  allQs.push_back(NULL);
}

void TerrainHoneycomb::dumpCellRidge(size_t cellID, Ridge ridge)
{
  cellRidges[cellID].push_back(ridge);
}

Q* TerrainHoneycomb::getQ(size_t idx)
{
  return allQs[idx];
}

std::vector<Ridge> TerrainHoneycomb::getCellRidges(size_t nodeID)
{
  return cellRidges[nodeID];
}