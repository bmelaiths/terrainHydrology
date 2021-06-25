import numpy as np
import math

import typing

def pointInConvexPolygon(point: typing.Tuple[float,float], vertices: np.ndarray, pivotPoint: typing.Tuple[float,float]) -> bool:
  """Determine if a point is within a convex polygon

  Algorithm derived from Wolfram: _An Efficient Test for a Point to Be in a Convex Polygon: http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/

  .. todo::
    This function must be reworked. The fact that it requires a known pivot point is
    a hack. There is certainly a better way to do this.

  :param point: The point to test
  :param point: tuple[float,float]
  :param vertices: The polygon to test
  :type vertices: numpy.ndarray(2,n)
  :param pivotPoint: A point that is known to be within the polygon
  :type pivotPoint: tuple[float,float]
  :return: True if the point is within the polygon, False otherwise
  :rtype: bool
  """
  previousSignPositive = None

  # make the test point the new origin
  vertices = np.subtract(vertices, point)

  for i in range(len(vertices)):
    v0, v1 = vertices[i], vertices[(i+1)%len(vertices)]
    # a_i = X_{i+1} Y_{i} - X_{i} Y_{i+1}
    ai = v1[0]*v0[1] - v0[0]*v1[1]
    if previousSignPositive is None:
      previousSignPositive = (ai > 0)
      continue
    if previousSignPositive != (ai > 0):
      # ai must always have the same sign
      return False
  return True

def convexPolygonArea(pivotPoint: typing.Tuple[float,float], vertices: np.ndarray) -> bool:
  """Determine the area of a convex polygon

  .. todo::
    This function must be reworked. The fact that it requires a known pivot point is
    a hack. There is certainly a better way to do this.

  :param pivotPoint: A point that is known to be within the polygon
  :type pivotPoint: tuple[float,float]
  :param vertices: An array of points that defines the shape
  :type vertices: np.ndarray[[float,float]]
  :return: The area of the polygon
  :rtype: float
  """
  vertices = np.subtract(vertices, pivotPoint)
  vertices = list(vertices)
  vertices.sort(key = lambda coord: math.atan2(coord[1],coord[0]))
  vertices = np.array(vertices)

  x = vertices[:,0]
  y = vertices[:,1]

  # This one-line implementation of Gauss's Shoelace Formula
  # written by stackoverflow user Mahdi
  return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

# Borrowed , all of it
def segments_distance(a: typing.Tuple[float,float], b: typing.Tuple[float,float], c: typing.Tuple[float,float], d: typing.Tuple[float,float]) -> float:
  """Distance between two line segments

      One segment is a to b, and the other is c to d

      :param a: Point A of the first segment
      :type a: tuple[float,float]
      :param b: Point B of the first segment
      :type b: tuple[float,float]
      :param c: Point C of the second segment
      :type c: tuple[float,float]
      :param d: Point D of the second segment
      :type d: tuple[float,float]
      :return: The distance between the two segments
      :rtype: float
  """
  #print(a[0],a[1],b[0],b[1],c[0],c[1],d[1],d[0])
  #print(segments_distance_internal(a[0],a[1],b[0],b[1],c[0],c[1],d[1],d[0]))
  return segments_distance_internal(a[0],a[1],b[0],b[1],c[0],c[1],d[0],d[1])

def segments_distance_internal(x11: float, y11: float, x12: float, y12: float, x21: float, y21: float, x22: float, y22: float) -> float:
  """Distance between two segments in a plane, used by :func:`segments_distance`

      One segment is (x11, y11) to (x12, y12), and
      the other is   (x21, y21) to (x22, y22)
  """
  if segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22): return 0
  # try each of the 4 vertices w/the other segment
  distances = []
  distances.append(point_segment_distance(x11, y11, x21, y21, x22, y22))
  distances.append(point_segment_distance(x12, y12, x21, y21, x22, y22))
  distances.append(point_segment_distance(x21, y21, x11, y11, x12, y12))
  distances.append(point_segment_distance(x22, y22, x11, y11, x12, y12))
  return min(distances)

def segments_intersect(x11: float, y11: float, x12: float, y12: float, x21: float, y21: float, x22: float, y22: float) -> bool:
  """ Whether two segments in the plane intersect

      One segment is (x11, y11) to (x12, y12), and the
      other is (x21, y21) to (x22, y22).

      :return: True if the line segments intersect, False if otherwise
      :rtype: bool
  """
  dx1 = x12 - x11
  dy1 = y12 - y11
  dx2 = x22 - x21
  dy2 = y22 - y21
  delta = dx2 * dy1 - dy2 * dx1
  if delta == 0: return False  # parallel segments
  s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta
  t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / (-delta)
  return (0 <= s <= 1) and (0 <= t <= 1)

# I think this finds the distance between a point and a line segment
def point_segment_distance(px: float, py: float, x1: float, y1: float, x2: float, y2: float):
  """Returns the distance between a point and a line segment

      :param px: The x location of the coordinate
      :type px: float
      :param py: The y location of the coordinate
      :type py: float
      :param x1: The x location of the first end of the segment
      :type x1: float
      :param y1: The y location of the first end of the segment
      :type y1: float
      :param x2: The x location of the second end of the segment
      :type x2: float
      :param y2: The y location of the second end of the segment
      :type y2: float
      :return: The distance
      :rtype: float
  """
  dx = x2 - x1
  dy = y2 - y1
  if dx == dy == 0:  # the segment's just a point
    return math.hypot(px - x1, py - y1)

  # Calculate the t that minimizes the distance.
  t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  # See if this represents one of the segment's
  # end points or a point in the middle.
  if t < 0:
    dx = px - x1
    dy = py - y1
  elif t > 1:
    dx = px - x2
    dy = py - y2
  else:
    near_x = x1 + t * dx
    near_y = y1 + t * dy
    dx = px - near_x
    dy = py - near_y

  return math.hypot(dx, dy)

def point_segment_distance_is_endpoint(px: float, py: float, x1: float, y1: float, x2: float, y2: float) -> typing.Tuple[float,bool]:
  """Same as :func:`point_segment_distance`, but also returns a bool that indicates whether or not the nearest point on the segment is just an endpoint

      :param px: The x location of the coordinate
      :type px: float
      :param py: The y location of the coordinate
      :type py: float
      :param x1: The x location of the first end of the segment
      :type x1: float
      :param y1: The y location of the first end of the segment
      :type y1: float
      :param x2: The x location of the second end of the segment
      :type x2: float
      :param y2: The y location of the second end of the segment
      :type y2: float
      :return: A tuple with the distance, and True if the nearest part of the segment is an endpoint, False if otherwise
      :rtype: tuple[float,bool]
  """
  dx = x2 - x1
  dy = y2 - y1
  if dx == dy == 0:  # the segment's just a point
    return (math.hypot(px - x1, py - y1), True)

  # Calculate the t that minimizes the distance.
  t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  isEndpoint = True

  # See if this represents one of the segment's
  # end points or a point in the middle.
  if t < 0:
    dx = px - x1
    dy = py - y1
  elif t > 1:
    dx = px - x2
    dy = py - y2
  else:
    near_x = x1 + t * dx
    near_y = y1 + t * dy
    dx = px - near_x
    dy = py - near_y
    isEndpoint = False

  return (math.hypot(dx, dy), isEndpoint)

def point_segment_distance_tuple(p: typing.Tuple[float,float], a: typing.Tuple[float,float], b: typing.Tuple[float,float]) -> float:
    """Just a front for :func:`point_segment_distance`, but allows tuple arguments

      :param p: The point
      :type p: tuple[float,float]
      :param a: The first endpoint of the segment
      :type a: tuple[float,float]
      :param b: The second endpoint of the segment
      :type b: tuple[float,float]
      :return: The distance
      :rtype: float
    """
    return point_segment_distance(p[0],p[1],a[0],a[1],b[0],b[1])
    
def segments_intersect_tuple(a1: typing.Tuple[float,float], a2: typing.Tuple[float,float], b1: typing.Tuple[float,float], b2: typing.Tuple[float,float]) -> bool:
    """A front for :func:`segments_intersect`, but with tuple arguments
    
    :param a1: Point 1 of the first segment
    :type a1: tuple[float,float]
    :param a2: Point 2 of the first segment
    :type a2: tuple[float,float]
    :param b1: Point 1 of the second segment
    :type b1: tuple[float,float]
    :param b2: Point 2 of the second segment
    :type b2: tuple[float,float]
    :return: True if the segments intersect, False if otherwise
    :rtype: bool
    """
    return segments_intersect(a1[0],a1[1],a2[0],a2[1],b1[0],b1[1],b2[0],b2[1])

def window(seq, n=2): ##Borrowed as is
    """Returns a sliding window (of width n) over data from the iterable

       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

       Comments indicate that this function is "borrowed as is", but I
       don't know where from.

       .. todo::
          This function is no longer used. The algorithm that used it
          should either be brought back in a future version, or this
          function should be deleted.
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
def clean_asin(sinAngle): ## Borrowed but modified
    return math.asin(min(1,max(sinAngle,-1)))

# This is the least intuitive way to calculate Euclidean distance, but I respect it
def distance(a: np.ndarray, b: np.ndarray) -> float:
   """Euclidean distance between two points

      :param a: Point A
      :type a: np.ndarray[float,float]
      :param b: Point B
      :type b: np.ndarray[float,float]
      :return: The distance between A and B
      :rtype: float
   """
   return np.linalg.norm( np.subtract(a , b))

def projection(p: np.ndarray, u: np.ndarray, v: np.ndarray):
    """I don't actually know what this function does 🙁

    .. todo::
       This function does not appear to be used anywhere. If it is
       indeed useless, then it should be removed.
    """
    n = np.subtract(u,v)
    nnorm = np.linalg.norm(n, 2)
    n = n/nnorm
    ret = np.dot(np.subtract(p,v), n)
    proj = ret/nnorm
    if proj >1 :
        proj=1
    if proj <0 :
        proj =0
    return proj*n
