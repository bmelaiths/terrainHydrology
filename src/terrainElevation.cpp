#include "terrainElevation.hpp"

#define FLOAT_TOLERANCE 0.001

float distance(const Point& p0, const Point& p1)
{
  // traditional Euclidean distance calculation
  return sqrtf32(
    powf32(p1.x()-p0.x(),2)+
    powf32(p1.y()-p0.y(),2)
  );
}

endpointAndDistance point_segment_distance
(
  const Point& tLoc, const Point& q0Loc, const Point& q1Loc
)
{
  float dx = q1Loc.x() - q0Loc.x();
  float dy = q1Loc.y() - q0Loc.y();
  if (abs(dx) < FLOAT_TOLERANCE && abs(dy) < FLOAT_TOLERANCE)
  { // the segment is just a point
    return {
      .dist = hypotf32(tLoc.x() - q0Loc.x(), tLoc.y() - q0Loc.y()),
      .isEndpoint = true
    };
  }

  // calculate the t that minimizes the distance
  float t =
    (((tLoc.x()-q0Loc.x()) * dx) + ((tLoc.y() - q0Loc.y()) * dy)) /
    (dx * dx + dy * dy)
  ;

  bool isEndpoint = true;

  // See if this represents one of the segment's endpoints,
  // or a point in the middle
  if (t < 0)
  {
    dx = tLoc.x() - q0Loc.x();
    dy = tLoc.y() - q0Loc.y();
  }
  else if (t > 1)
  {
    dx = tLoc.x() - q1Loc.x();
    dy = tLoc.y() - q1Loc.y();
  }
  else
  {
    float near_x = q0Loc.x() + t * dx;
    float near_y = q0Loc.y() + t * dy;
    dx = tLoc.x() - near_x;
    dy = tLoc.y() - near_y;
    isEndpoint = false;
  }
  
  return {.dist=hypotf32(dx,dy), .isEndpoint=isEndpoint};
}

float lerpRidge(Q *q0, Q *q1, T& t, float dist)
{
  float result =
    q0->getElevation() +
    (
      ( // this gets the distance along the ridge that the t corresponds to
        sqrtf32(
          pow(distance(q0->getPosition(),t.getLoc()),2) -
          pow(dist,2)
        ) /
        distance(q0->getPosition(),q1->getPosition())
      ) *
      (q1->getElevation() - q0->getElevation()) // the slope
    )
  ;
  // If the t is really close to the first endpoint, the calculation
  // may fail. In that case, return the elevation of that primitive
  return isnan(result) ? q0->getElevation() : result;
}

float computePrimitiveElevation
(
  T& t, Hydrology& hydrology, TerrainHoneycomb& cells, Terrain& ts,
  std::vector<cv::Point>& contour, float resolution,
  GEOSContextHandle_t geosContext
)
{
  /*
    The point of this loop is to find the ridge that the primitive
    is closest to, and to interpolate the elevation of the point
    along the ridge that the primitive corresponds to.
  */
  std::vector<Ridge> ridges = cells.getCellRidges(t.getCellID());
  float closestRidgeDist = -1.0;
  float ridgeElevation = -1.0;
  for (Ridge& ridge : ridges)
  {
    if (ridge.getSize() < 2)
    {
      Q *q = ridge.getPoint0();
      float dist = distance(q->getPosition(), t.getLoc());
      if (closestRidgeDist < 0 || dist < closestRidgeDist)
      {
        closestRidgeDist = dist;
        ridgeElevation = ridge.getPoint0()->getElevation();
      }
      continue;
    }

    Q *q0, *q1;
    q0 = ridge.getPoint0();
    q1 = ridge.getPoint1();
    endpointAndDistance res = point_segment_distance(
      t.getLoc(),q0->getPosition(),q1->getPosition()
    );
    float dist = res.dist;
    bool isEndpoint = res.isEndpoint;

    if (closestRidgeDist > 0 && closestRidgeDist < dist)
    {
      continue;
    }

    if (isEndpoint)
    {
      if
      (
        distance(q0->getPosition(),t.getLoc()) <
        distance(q1->getPosition(),t.getLoc())
      )
      {
        closestRidgeDist = distance(q0->getPosition(),t.getLoc());
        ridgeElevation = q0->getElevation();
      }
      else
      {
        closestRidgeDist = distance(q1->getPosition(),t.getLoc());
        ridgeElevation = q1->getElevation();
      }
    }
    else
    {
      closestRidgeDist = dist;
      ridgeElevation = lerpRidge(q0, q1, t, dist);
    }
  }

  /*
    If the shore is closer than the closest ridge, then use that
    instead.
  */
  float distToGamma = resolution * (float) cv::pointPolygonTest(
    contour,
    cv::Point2f(
      (float) t.getLoc().x() / resolution,
      (float) t.getLoc().y() / resolution
    ),
    true
  );
  if (closestRidgeDist < 0 || distToGamma < closestRidgeDist)
  {
    closestRidgeDist = distToGamma;
    ridgeElevation = 0;
  }
  

  double closestRiverDist;
  GEOSGeometry *projectedPoint;

  GEOSGeometry *point = GEOSGeom_createPointFromXY_r(
    geosContext, t.getLoc().x(),t.getLoc().y()
  );
  Primitive* node = hydrology.getNodeP(t.getCellID());

  if (node->numRivers() > 0)
  {
    GEOSGeometry *closestRiver = NULL;
    for (GEOSGeometry *river : node->getRivers())
    {
      double riverDistance;
      GEOSDistance_r(geosContext, point, river, &riverDistance);
      if (closestRiver == NULL || riverDistance < closestRiverDist)
      {
        closestRiverDist = riverDistance;
        closestRiver = river;
      }
    }

    double projectedDistance = GEOSProject_r(geosContext, closestRiver, point);
    projectedPoint = GEOSInterpolate_r(
      geosContext, closestRiver, projectedDistance
    );
    //this line is probably redundant
    GEOSDistance_r(geosContext, point, projectedPoint, &closestRiverDist);
  }
  else
  {
    GEOSCoordSequence *pointSeq = GEOSCoordSeq_create_r(geosContext, 1,3);
    GEOSCoordSeq_setXYZ_r(
      geosContext, pointSeq, 0,
      node->getLoc().x(), 
      node->getLoc().y(),
      node->getElevation()
    );
    projectedPoint = GEOSGeom_createPoint_r(geosContext, pointSeq);
    GEOSDistance_r(geosContext, point, projectedPoint, &closestRiverDist);
  }

  if
  (
    abs(closestRiverDist) < FLOAT_TOLERANCE &&
    abs(closestRidgeDist) < FLOAT_TOLERANCE
  )
  {
    closestRiverDist = 1;
  }

  // return lerped elevation
  double projectedZ;
  GEOSGeomGetZ_r(geosContext, projectedPoint, &projectedZ);
  return
    projectedZ * (closestRidgeDist/(closestRidgeDist+closestRiverDist)) +
    ridgeElevation*(closestRiverDist/(closestRidgeDist+closestRiverDist))
  ;
}