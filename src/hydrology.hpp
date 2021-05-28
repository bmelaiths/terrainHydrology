#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include <vector>
#include <list>
#include <stdio.h>
#include <stdint.h>

#include "point.hpp"
#include "kdtree.hpp"
#include "forest.hpp"

class Hydrology;

/**
 * @brief A Hydrology primitive. Represents a certain stretch of river
 * 
 */
class Primitive
{
  private:
  size_t id;
  Primitive *parent;
  std::vector<Primitive*> children;
  Point loc;
  float elevation;
  int priority;
  int contourIndex;

  public:
  friend Hydrology;

  public:
  Primitive();
  /**
   * @brief This constructor is for mouth nodes
   * 
   * @param id The ID of this node. See `Hydrology` for this value’s significance
   * @param loc The location of the node in meters
   * @param elevation The node’s elevation in meters
   * @param priority The priority of the node. See `Hydrology` for this value’s significance
   * @param contourIndex This is the index in the contour vector that is closest to this node
   */
  Primitive
  ( //for mouth nodes
    size_t id, Point loc, float elevation, int priority, int contourIndex
  );
  /**
   * @brief This constructor is for regular nodes
   * 
   * @param id The ID of this node. See `Hydrology` for this value’s significance
   * @param parent A pointer to the node's parent
   * @param loc The location of the node in meters
   * @param elevation The node’s elevation in meters
   * @param priority The priority of the node. See `Hydrology` for this value’s significance
   */
  Primitive
  ( //for regular nodes
    size_t id, Primitive *parent, Point loc, float elevation, int priority
  );

  /**
   * @brief The size of this node's binary representation
   * 
   * @return size_t The size (in bytes?)
   */
  size_t binarySize();
  /**
   * @brief Writes a binary representation of this primitive to the buffer
   * 
   * @param buffer This buffer should be of the appropriate size (get size with binarySize())
   */
  void toBinary(uint8_t *buffer);

  /**
   * @brief Get the ID of this node
   * 
   * @return size_t 
   */
  size_t getID();
  /**
   * @brief Get a pointer to this primitive's parent
   * 
   * @return Primitive* A pointer to this primitive's parent
   */
  Primitive* getParent();
  /**
   * @brief Check whether or not this primitive has a parent or not
   * 
   * @return true This is a regular node
   * @return false This is a mouth node, without a parent
   */
  bool hasParent();
  /**
   * @brief Gets pointers to this node's children
   * 
   * @return std::vector<Primitive*> A vector containing pointers to this node's child nodes
   */
  std::vector<Primitive*> getChildren();
  /**
   * @brief The number of children that this node has
   * 
   * @return size_t 
   */
  size_t numChildren();
  /**
   * @brief Gets the location of this node
   * 
   * @return Point 
   */
  Point getLoc();
  /**
   * @brief Gets the elevation of this node
   * 
   * @return float The elevation, in meters
   */
  float getElevation();
  /**
   * @brief Get the priority
   * 
   * @return int 
   */
  int getPriority();
  /**
   * @brief If this node is on the coast, this is the index in `Shore` that is closest to this node
   * 
   * @return int 
   */
  int getContourIndex();
};

/**
 * @brief Represents an edge, containing the two ends of which it consists
 * 
 */
class Edge
{
  public:
  Primitive *node0, *node1;

  public:
  Edge(Primitive *node0, Primitive *node1);
};

/**
 * @brief This class represents the network of rivers that flow over the land
 * 
 * A HydrologyNetwork is basically a forest of trees, with each tree representing a river that merges and drains into the ocean through a single mouth node.
 * 
 * A HydrologyNetwork is empty when it is instantiated. The network is built incrementally using the addNode() methods. Edges connect all nodes.
 * 
 * It should be noted that each node is associated with an integer ID. This ID is strictly the order that each node was added in, starting at 0. Nodes cannot be removed, thus range(len(hydrology)) will iterate over all the nodes in the network.
 * 
 */
class Hydrology
{
  private:
  std::vector<Primitive*> indexedNodes;
  Forest<Primitive*> trees;

  public:
  Hydrology() = default;
  /**
   * @brief Constructs the hydrology network with the specified area
   * 
   * @param lowerLeft The lower left corner of the area for the hydrology network
   * @param upperRight The upper right corner of the area for the hydrology network
   * @param edgeLength The edge length parameter of the simulation
   */
  Hydrology(Point lowerLeft, Point upperRight, float edgeLength);
  ~Hydrology();
  Hydrology(const Hydrology& other);
  Hydrology(Hydrology&& other);

  Hydrology& operator=(const Hydrology& other);
  Hydrology& operator=(Hydrology&& other);

/**
 * @brief Add a mouth node to the network
 * 
 * @param loc The location of the new node
 * @param elevation The elevation of the new node
 * @param priority The priority of the new node
 * @param contourIndex The index of the node's location on the shore
 * @return Primitive* The node created
 */
  Primitive* addMouthNode(
    Point loc, float elevation, int priority, int contourIndex
  );
  /**
   * @brief 
   * 
   * @param loc The location of the new node
   * @param elevation The elevation of the new node
   * @param priority The priority of the new node
   * @param parent The ID of the parent node
   * @return Primitive* 
   */
  Primitive* addRegularNode(
    Point loc, float elevation, int priority, size_t parent
  );

  /**
   * @brief Acquires a lock on the specified area
   * 
   * @param loc The center of the area to search
   * @param radius The width and height of the area to search
   * @return AreaLock
   */
  AreaLock lockArea(Point loc, float radius);
  /**
   * @brief Returns all edges with one or both nodes within the specified area
   * 
   * NOTE: This method double-counts edges where both nodes are within the search range
   * 
   * @param loc The center of the area to search
   * @param radius The width and height of the area to search
   * @return std::vector<Edge> 
   */
  std::vector<Edge> queryArea(Point loc, float radius);

  /**
   * @brief Get the Node with the specified index
   * 
   * @param idx The index of the node to find
   * @return Primitive The node at that index
   */
  Primitive getNode(size_t idx);
  /**
   * @brief Returns the total number of nodes in the network
   * 
   * @return size_t 
   */
  size_t numNodes();

  /**
   * @brief Writes a binary representation of this network to the specified stream
   * 
   * @param stream 
   */
  void writeBinary(FILE *stream);
};

#endif