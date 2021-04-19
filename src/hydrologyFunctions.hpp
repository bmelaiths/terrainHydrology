#ifndef HYDROFUNC_H
#define HYDROFUNC_H

#include "hydrology.hpp"
#include "hydrologyParameters.hpp"

Primitive selectNode(HydrologyParameters& params);

bool isAcceptablePosition(Point testLoc, HydrologyParameters& params);

float coastNormal(Primitive candidate, HydrologyParameters& params);

Point pickNewNodeLoc(Primitive candidate, HydrologyParameters& params);

void tao(Primitive node, HydrologyParameters& params);

void beta
(Primitive node, int priority, HydrologyParameters& params);

void ruleBase(Primitive candidate, HydrologyParameters& params);

void pa(Primitive candidate, HydrologyParameters& params);

void pc(Primitive candidate, HydrologyParameters& params);

void ps(Primitive candidate, HydrologyParameters& params);

void alpha(Primitive candidate, HydrologyParameters& params);

#endif