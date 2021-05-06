#ifndef HYDROFUNC_H
#define HYDROFUNC_H

#include "hydrology.hpp"
#include "hydrologyParameters.hpp"

class LockedPoint
{
  private:
  Point point;
  Lock lock;

  public:
  LockedPoint(Point point, Lock lock);
  void release();
  Point getLoc();
};

Primitive selectNode(HydrologyParameters& params);

bool isAcceptablePosition(Point testLoc, float radius, size_t parentID, HydrologyParameters& params);

float coastNormal(Primitive candidate, HydrologyParameters& params);

LockedPoint pickNewNodeLoc(Primitive candidate, HydrologyParameters& params);

void tao(Primitive node, HydrologyParameters& params);

void beta
(Primitive node, int priority, HydrologyParameters& params);

void ruleBase(Primitive candidate, HydrologyParameters& params);

void pa(Primitive candidate, HydrologyParameters& params);

void pc(Primitive candidate, HydrologyParameters& params);

void ps(Primitive candidate, HydrologyParameters& params);

void alpha(Primitive candidate, HydrologyParameters& params);

#endif