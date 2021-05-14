#include "hydrology.hpp"

#include "floatEndian.hpp"

#include <string.h>
#include <stdlib.h>

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

Edge::Edge(Primitive *node0, Primitive *node1)
: node0(node0), node1(node1)
{
}

Hydrology::Hydrology(Point lowerLeft, Point upperRight)
{
  float dimension;
  if (upperRight.x() - lowerLeft.x() > upperRight.y() - lowerLeft.y())
  {
    dimension = (upperRight.x() - lowerLeft.x()) / 10;
  }
  else
  {
    dimension = (upperRight.y() - lowerLeft.y()) / 10;
  }
  trees = Forest<Primitive*>(lowerLeft, upperRight, dimension);
}

Hydrology::~Hydrology()
{
  for (Primitive *node : indexedNodes)
  {
    delete node;
  }
}

Hydrology::Hydrology(const Hydrology& other)
: indexedNodes(other.indexedNodes), trees(other.trees)
{}

Hydrology::Hydrology(Hydrology&& other)
: indexedNodes(std::move(other.indexedNodes)), trees(std::move(other.trees))
{
  other.indexedNodes = std::vector<Primitive*>();
  other.trees = Forest<Primitive*>();
}

Hydrology& Hydrology::operator=(const Hydrology& other)
{
  if (this == &other)
  {
    return *this;
  }

  indexedNodes = other.indexedNodes;
  trees = other.trees;

  return *this;
}

Hydrology& Hydrology::operator=(Hydrology&& other)
{
  if (this == &other)
  {
    return *this;
  }

  indexedNodes = std::move(other.indexedNodes);
  trees = std::move(other.trees);

  other.indexedNodes = std::vector<Primitive*>();
  other.trees = Forest<Primitive*>();

  return *this;
}

Primitive* Hydrology::addMouthNode
(Point loc, float elevation, int priority, int contourIndex)
{
  Primitive *node = new Primitive(indexedNodes.size(), loc, elevation, priority, contourIndex);

  trees.insert(loc, node);
  indexedNodes.push_back(node);

  return node;
}

Primitive* Hydrology::addRegularNode
(Point loc, float elevation, int priority, size_t parent)
{
  Primitive *node = new Primitive(indexedNodes.size(), indexedNodes[parent], loc, elevation, priority);

  trees.insert(loc, node);
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

  return node;
}

AreaLock Hydrology::lockArea(Point loc, float radius)
{
  return trees.lockArea(loc, radius);
}

//note: this method may double-count edges that have both ends in the area
std::vector<Edge> Hydrology::queryArea(Point loc, float radius)
{
  std::vector<Primitive*> closeIdxes = trees.searchRange(loc, radius);

  std::vector<Edge> edges;
  for (Primitive *closeIdx : closeIdxes)
  {
    // Primitive idxNode = indexedNodes[closeIdx];
    if (closeIdx->hasParent())
    {
      edges.push_back(Edge(closeIdx, closeIdx->getParent()));
    }
    for (size_t i = 0; i < closeIdx->numChildren(); i++)
    {
      edges.push_back(Edge(closeIdx, closeIdx->getChildren()[i]));
    }
  }
  return edges;
}

Primitive Hydrology::getNode(size_t idx) {
  return *indexedNodes[idx];
}

size_t Hydrology::numNodes()
{
  return indexedNodes.size();
}

void Hydrology::writeBinary(FILE *stream)
{
  uint64_t numPrimitives = htobe64((uint64_t) indexedNodes.size());
  fwrite(&numPrimitives, sizeof(uint64_t), 1, stream);
  fflush(stream);

  uint8_t *buffer;
  for (size_t i = 0; i < indexedNodes.size(); i++)
  {
    size_t size = indexedNodes[i]->binarySize();
    buffer = new uint8_t[size];

    indexedNodes[i]->toBinary(buffer);

    fwrite(buffer, size, 1, stream);

    delete buffer;

    fflush(stream);
  }
}