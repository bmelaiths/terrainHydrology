#include "hydrologyFunctions.hpp"
#include "point.hpp"

#include <random>

#include <opencv2/imgproc.hpp>

#define FLOAT_THRESH 0.1f

class ComparePrimitives
{
  private:
  HydrologyParameters& params;

  public:
  ComparePrimitives(HydrologyParameters& params)
  : params(params)
  {}
  ~ComparePrimitives() {}
  bool operator() (Primitive *a, Primitive *b) {
    return (a->getPriority() > b->getPriority());
  }
};

Primitive selectNode(HydrologyParameters& params) {
  float lowestCandidateZ = params.candidates[0]->getElevation();
  for (size_t i = 0; i < params.candidates.size(); i++)
  {
    if
    (
      params.candidates[i]->getElevation()
      <
      lowestCandidateZ
    )
    {
      lowestCandidateZ = params.candidates[i]->getElevation();
    }
  }

  std::vector<Primitive*> subselection;
  for (size_t i = 0; i < params.candidates.size(); i++)
  {
    if
    (
      params.candidates[i]->getElevation()
      <
      (lowestCandidateZ + params.zeta)
    )
    {
      subselection.push_back(params.candidates[i]);
    }
  }

  std::sort(
    subselection.begin(),
    subselection.end(),
    ComparePrimitives(params)
  );

  size_t idx = 0;
  Primitive lowestElevation = *subselection[idx];
  while
  (
    idx < subselection.size() &&
    subselection[idx]->getPriority() == subselection[0]->getPriority()
  )
  {
    if
    (
      subselection[idx]->getElevation()
      <
      lowestElevation.getElevation()
    )
    {
      lowestElevation = *subselection[idx];
    }
    idx++;
  }

  return lowestElevation;
}

float point_segment_distance(float px, float py, float x1, float y1, float x2, float y2) {
  float dx = x2 - x1;
  float dy = y2 - y1;

  if (abs(dx) < FLOAT_THRESH && abs(dy) < 0)  // the segment's just a point
  {
    return hypotf32(px - x1, py - y1);
  }

  // Calculate the t that minimizes the distance.
  float t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

  // See if this represents one of the segment's
  // end points or a point in the middle.
  if (t < 0)
  {
    dx = px - x1;
    dy = py - y1;
  }
  else if (t > 1)
  {
    dx = px - x2;
    dy = py - y2;
  }
  else
  {
    float near_x = x1 + t * dx;
    float near_y = y1 + t * dy;
    dx = px - near_x;
    dy = py - near_y;
  }

  return hypotf32(dx, dy);
}

bool isAcceptablePosition(Point testLoc, float radius, size_t parentID, HydrologyParameters& params) {
  if (testLoc.x() < 0)
  {
    return false;
  }
  
  float distToGamma = params.resolution * (float) cv::pointPolygonTest(
    params.contour,
    cv::Point2f(
      (float) testLoc.x() / params.resolution,
      (float) testLoc.y() / params.resolution
    ),
    true
  );
  if (distToGamma < params.edgeLength * params.eta)
  {
    return false;
  }
  
  std::vector<Edge> edges = params.hydrology.queryArea(
    testLoc, radius
  );
  for (Edge edge : edges)
  {
    float dist = point_segment_distance(
      testLoc.x(), testLoc.y(),
      edge.node0.getLoc().x(), edge.node0.getLoc().y(),
      edge.node1.getLoc().x(), edge.node1.getLoc().y()
    );
    if (dist < params.sigma * params.edgeLength)
    {
      return false;
    }
  }
  return true;
}

float coastNormal(Primitive candidate, HydrologyParameters& params) {
  Point p1(
    params.contour[candidate.getContourIndex()+3].x,
    params.contour[candidate.getContourIndex()+3].y
  );
  Point p2(
    params.contour[candidate.getContourIndex()-3].x,
    params.contour[candidate.getContourIndex()-3].y
  );
  float theta = atan2(p2.y() - p1.y(), p2.x() - p1.x());
  return theta + M_PI_2f32;
}

LockedPoint::LockedPoint(Point point, Lock lock)
: point(point), lock(lock)
{}

void LockedPoint::release() {
  lock.release();
}

Point LockedPoint::getLoc() {
  return point;
}

LockedPoint pickNewNodeLoc(Primitive candidate, HydrologyParameters& params) {
  float angle;
  // if candidate has no parent
  if (!candidate.hasParent())
  {
    angle = coastNormal(candidate, params);
  }
  else
  {
    // else 'angle will be the direction of the river
    angle = atan2(
      candidate.getLoc().y() - candidate.getParent()->getLoc().y(),
      candidate.getLoc().x() - candidate.getParent()->getLoc().x()
    );
  }

  // std::default_random_engine generator;
  // std::normal_distribution<float> distribution(0.0, params.riverAngleDev);

  float newAngle;
  Point newLoc(-1,-1);
  Lock lock;
  for (size_t i = 0; i < params.maxTries; i++)
  {
    newAngle = angle + params.distribution(params.generator);
    newLoc = Point(
      candidate.getLoc().x() + cosf32(newAngle) * params.edgeLength,
      candidate.getLoc().y() + sinf32(newAngle) * params.edgeLength
    );
    lock = params.hydrology.lockArea(newLoc, 2 * params.edgeLength);
    if (isAcceptablePosition(newLoc, 2 * params.edgeLength, candidate.getID(), params))
    {
      break;
    }
    else
    {
      lock.release();
      newLoc = Point(-1,-1);
    }
  }

  return LockedPoint(newLoc, lock);
}

void tao(Primitive node, HydrologyParameters& params) {
  #pragma omp critical
  {
  size_t candidateIdx = 0;
  for (size_t i = 0; i < params.candidates.size(); i++)
  {
    if (params.candidates[i]->getID() == node.getID())
    {
      candidateIdx = i;
      params.candidates.erase(params.candidates.begin() + candidateIdx);
      break;
    }
  }
  }
}

void beta
(Primitive node, int priority, HydrologyParameters& params)
{
  LockedPoint lockedPoint = pickNewNodeLoc(node, params);
  Point newLoc = lockedPoint.getLoc();
  if (newLoc.x() != -1)
  {
    float slope =
      params.slopeRate *
      params.riverSlope.get(node.getLoc().x(), node.getLoc().y()) / 255
    ;
    float newZ = node.getElevation() + slope * params.edgeLength;
    #pragma omp critical
    {
      params.candidates.push_back(
        params.hydrology.addRegularNode(
          newLoc, newZ, priority, node.getID()
        )
      );
    }
    lockedPoint.release();
  }
  else
  {
    tao(node, params);
  }
}

void ruleBase(Primitive candidate, HydrologyParameters& params) {
  //int numBranches = (int) (rand() * 4) + 1;
  const int numBranches = 5;
  for (int i = 0; i < numBranches; i++)
  {
    beta(candidate, candidate.getPriority(), params);
  }
}

void pa(Primitive candidate, HydrologyParameters& params) {
  beta(candidate, candidate.getPriority(), params);
  beta(candidate, candidate.getPriority() - 1, params);
}

void pc(Primitive candidate, HydrologyParameters& params) {
  beta(candidate, candidate.getPriority(), params);
}

void ps(Primitive candidate, HydrologyParameters& params) {
  beta(candidate, candidate.getPriority() - 1, params);
  beta(candidate, candidate.getPriority() - 1, params);
}

void alpha(Primitive candidate, HydrologyParameters& params) {
  if (candidate.getPriority() == 1)
  {
    ruleBase(candidate, params);
  }
  else
  {
    float pval = rand();
    if (pval <= params.Pa)
    {
      pa(candidate, params);
    }
    else if (pval <= (params.Pa + params.Pc))
    {
      pc(candidate, params);
    }
    else
    {
      ps(candidate, params);
    }
  }
}