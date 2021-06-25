#ifndef HYDROFUNC_H
#define HYDROFUNC_H

#include "hydrology.hpp"
#include "hydrologyParameters.hpp"

/**
 * @brief A simple struct that pairs a point in a locked area, and the lock for that area
 * 
 */
class LockedPoint
{
  private:
  Point point;
  AreaLock lock;

  public:
  /**
   * @brief Construct a new Locked Point object
   * 
   * @param point A point within the locked area
   * @param lock The lock on that area
   */
  LockedPoint(Point point, AreaLock lock);
  /**
   * @brief Releases the lock on the area
   * 
   */
  void release();
  /**
   * @brief Get the point that is within the area
   * 
   * @return Point 
   */
  Point getLoc();
};

/**
 * @brief Given a list of candidate nodes, this function selects the next node to expand
 * 
 * The selection is based on Genevaux et al §4.2.1
 * 
 * @param params All the information needed to select a node
 * @return Primitive 
 */
Primitive selectNode(HydrologyParameters& params);

/**
 * @brief Determines whether or not a point is acceptable as a new node’s location
 * 
 * @param testLoc The point in question
 * @param radius The radius to search
 * @param parentID The ID of the potential node's parent
 * @param params All the other information needed to make the determination
 * @return true The position is acceptable
 * @return false The position fails to meet the criteria for a new node
 */
bool isAcceptablePosition(Point testLoc, float radius, size_t parentID, HydrologyParameters& params);

/**
 * @brief Gets an angle that is (approximately) perpendicular to the coast
 * 
 * @param candidate The node to query
 * @param params All the information needed
 * @return float The angle in radians
 */
float coastNormal(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Tries to pick a new position for a candidate node
 * 
 * Tries to avoid the coast and crowding other nodes
 * 
 * @param candidate The node to start from
 * @param params The struct holding the relevant parameters
 * @return LockedPoint The location in this object is the new position. Don't forget to release the lock on the area
 */
LockedPoint pickNewNodeLoc(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Removes a node from the set of candidates (rule 3.2 from Table 1 of Genevaux et al)
 * 
 * @param node The node to remove
 * @param params The set (list) of candidates
 */
void tao(Primitive node, HydrologyParameters& params);

/**
 * @brief Instantiates a new node with a given priority
 * 
 * @param node The node to expand
 * @param priority The priority for the new node
 * @param params The parameter struct
 */
void beta
(Primitive node, int priority, HydrologyParameters& params);

/**
 * @brief Rule 1 from Table 1 in Genevaux et al
 * 
 * @param candidate The node to expand
 * @param params The parameter struct
 */
void ruleBase(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Rule 2.3 from Table 1 in Genevaux et al
 * 
 * @param candidate The node to expand
 * @param params The parameter struct
 */
void pa(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Rule 2.1 from Table 1 in Genevaux et al
 * 
 * @param candidate The node to expand
 * @param params The parameter struct
 */
void pc(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Rule 2.2 from Table 1 in Genevaux et al
 * 
 * @param candidate The node to expand
 * @param params The parameter struct
 */
void ps(Primitive candidate, HydrologyParameters& params);

/**
 * @brief Alpha node expansion rule
 * 
 * @param candidate The node to expand
 * @param params The parameter struct
 */
void alpha(Primitive candidate, HydrologyParameters& params);

#endif