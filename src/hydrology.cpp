#include "hydrology.hpp"

#include <string.h>
#include <stdlib.h>
#include <endian.h>

Primitive::Primitive()
:
id(0), parent(NULL), loc(Point(0,0)),
elevation(0.0f), priority(0), contourIndex(0)
{

}

Primitive::Primitive
(
  size_t id, Point loc, float elevation,
  int priority, int contourIndex
)
:
id(id), parent(NULL), loc(loc), elevation(elevation),
priority(priority), contourIndex(contourIndex)
{

}

Primitive::Primitive
(
  size_t id, Primitive *parentID, Point loc,
  float elevation, int priority
)
:
id(id), parent(parentID), loc(loc),
elevation(elevation), priority(priority)
{

}

size_t Primitive::binarySize()
{
  return
  (
    sizeof(size_t) * (2 + children.size()) +
    sizeof(uint8_t) * 1 +
    sizeof(float) * 3
  );
}

float float_tobe(float value) {
    union v {
        float f;
        uint32_t i;
    };

    union v val;

    val.f = value;
    val.i = htobe32(val.i);

    return val.f;
}

void Primitive::toBinary(uint8_t *buffer)
{
  size_t idx = 0;

  uint64_t ID = htobe64((uint64_t)id);
  memcpy(buffer + idx, &ID, sizeof(uint64_t));
  idx += sizeof(size_t);

  uint64_t PARENT;
  if (parent == NULL)
  {
    PARENT = ID;
  }
  else
  {
    PARENT = htobe64((uint64_t)parent->id);
  }
  
  memcpy(buffer + idx, &PARENT, sizeof(uint64_t));
  idx += sizeof(size_t);

  uint8_t numChildren = (uint8_t) children.size();
  memcpy(buffer + idx, &numChildren, sizeof(uint8_t));
  idx += sizeof(uint8_t);

  for (size_t child = 0; child < children.size(); child++)
  {
    uint64_t childID = htobe64((uint64_t)children[child]->id);
    memcpy(buffer + idx, &childID, sizeof(uint64_t));
    idx += sizeof(uint64_t);
  }

  float locX = float_tobe(loc.x());
  memcpy(buffer + idx, &locX, sizeof(float));
  idx += sizeof(float);

  float locY = float_tobe(loc.y());
  memcpy(buffer + idx, &locY, sizeof(float));
  idx += sizeof(float);

  float ELEV = float_tobe(elevation);
  memcpy(buffer + idx, &ELEV, sizeof(float));
  idx += sizeof(float);
}

size_t Primitive::getID() {return id;}

Primitive* Primitive::getParent() {return parent;}

bool Primitive::hasParent() {return parent != NULL;}

std::vector<Primitive*> Primitive::getChildren() {return children;}

size_t Primitive::numChildren() {return children.size();}

Point Primitive::getLoc() {return loc;}

float Primitive::getElevation() {return elevation;}

int Primitive::getPriority() {return priority;}

int Primitive::getContourIndex() {return contourIndex;}

Edge::Edge(Primitive node0, Primitive node1)
: node0(node0), node1(node1)
{
}

Primitive Hydrology::addMouthNode
(Point loc, float elevation, int priority, int contourIndex)
{
  primitiveStorage.push_back(
    Primitive(indexedNodes.size(), loc, elevation, priority, contourIndex)
  );

  tree.insert(loc, &primitiveStorage.back());
  indexedNodes.push_back(&primitiveStorage.back());

  return primitiveStorage.back();
}

Primitive Hydrology::addRegularNode
(Point loc, float elevation, int priority, size_t parent)
{
  primitiveStorage.push_back(
    Primitive(indexedNodes.size(), indexedNodes[parent], loc, elevation, priority)
  );
  
  Primitive* node = &primitiveStorage.back();

  tree.insert(loc, node);
  indexedNodes.push_back(node);

  node->getParent()->children.push_back(node);

  // Classify new leaf
  node->priority = 1;

  //loop invariants:
  //classifyNode is a pointer to a node whose parent had to be changed
  Primitive* classifyNode = node->getParent();
  while (true)
  {
    int maxPriority = classifyNode->getChildren()[0]->getPriority();
    int numMax = 1;
    for (size_t i = 0; i < classifyNode->numChildren(); i++)
    {
      if (classifyNode->getChildren()[i]->getPriority() == maxPriority)
      {
        numMax++;
        continue;
      }
      if (classifyNode->getChildren()[i]->getPriority() > maxPriority)
      {
        maxPriority = classifyNode->getChildren()[i]->getPriority();
        numMax = 1;
      }
    }

    if (numMax > 1 && classifyNode->getPriority() < numMax + 1)
    {
      classifyNode->priority = numMax + 1;
    }
    else if (classifyNode->getPriority() < numMax)
    {
      classifyNode->priority = numMax;
    }
    else
    {
      break;
    }
    
    if (!classifyNode->hasParent())
    {
      break;
    }
    
    classifyNode = classifyNode->getParent();
  }

  return *node;
}

//note: this method may double-count edges that have both ends in the area
std::vector<Edge> Hydrology::edgesWithinRadius(Point loc, float radius)
{
  std::vector<Primitive*> closeIdxes = tree.rangeSearch(loc, radius);

  std::vector<Edge> edges;
  for (Primitive* closeIdx : closeIdxes)
  {
    // Primitive idxNode = indexedNodes[closeIdx];
    if (closeIdx->hasParent())
    {
      edges.push_back(Edge(*closeIdx, *closeIdx->getParent()));
    }
    for (size_t i = 0; i < closeIdx->numChildren(); i++)
    {
      edges.push_back(Edge(*closeIdx, *closeIdx->getChildren()[i]));
    }
  }
  return edges;
}

Primitive Hydrology::getNode(size_t idx) {
  return *indexedNodes[idx];
}

void writeBinary(Hydrology hydrology, FILE *stream)
{
  uint64_t numPrimitives = htobe64((uint64_t) hydrology.indexedNodes.size());
  fwrite(&numPrimitives, sizeof(uint64_t), 1, stream);
  fflush(stream);

  uint8_t *buffer;
  for (size_t i = 0; i < hydrology.indexedNodes.size(); i++)
  {
    size_t size = hydrology.indexedNodes[i]->binarySize();
    buffer = new uint8_t[size];

    hydrology.indexedNodes[i]->toBinary(buffer);

    fwrite(buffer, size, 1, stream);

    delete buffer;

    fflush(stream);
  }
}