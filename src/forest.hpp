#ifndef FOREST_H
#define FOREST_H

#include "kdtree.hpp"
#include <omp.h>

class AreaLock
{
    private:
    std::vector<omp_lock_t*> locks;

    public:
    AreaLock() {}
    AreaLock(std::vector<omp_lock_t*> locks)
    : locks(locks)
    {
      for (omp_lock_t *lock : locks)
      {
        omp_set_lock(lock);
      }
      
    }
    void release()
    {
      for (omp_lock_t *lock : locks)
      {
          omp_unset_lock(lock);
      }
    }
};

template <typename T>
class Forest
{
  private:
  Point lowerLeft;
  Point upperRight;
  float dimensionLength;

  size_t xDimension, yDimension;
  KDTree<T>*  *forest;
  omp_lock_t*  forestLocks;

  size_t getTileX(float x);
  size_t getTileY(float y);
  size_t getIndex(size_t x, size_t y);

  public:
  Forest();
  ~Forest();

  void set(Point lowerLeft, Point upperRight, float dimension);

  void insert(Point loc, T idx);
  AreaLock lockArea(Point loc, float radius);
  std::vector<T> searchRange(Point loc, float radius);
};

template <typename T>
Forest<T>::Forest()
: lowerLeft(Point(0,0)), upperRight(Point(0,0)), forest(NULL), forestLocks(NULL)
{}

template <typename T>
void Forest<T>::set(Point lowerLeft, Point upperRight, float dimension)
{
  this->lowerLeft = lowerLeft;
  this->upperRight = upperRight;
  this->dimensionLength = dimension;

  xDimension = ceill((upperRight.x()-lowerLeft.x())/dimensionLength);
  yDimension = ceill((upperRight.y()-lowerLeft.y())/dimensionLength);

  forest =      new KDTree<T>*[yDimension * xDimension];
  forestLocks = new omp_lock_t[yDimension * xDimension];

  for (size_t y = 0; y < yDimension; y++)
  {
    for (size_t x = 0; x < xDimension; x++)
    {
      forest[getIndex(x, y)] = new KDTree<T>();
    }
  }
  for (size_t y = 0; y < yDimension; y++)
  {
    for (size_t x = 0; x < xDimension; x++)
    {
      omp_init_lock(forestLocks + (y * xDimension + x));
    }
  }
}

template <typename T>
Forest<T>::~Forest()
{
  if (forest != NULL)
  {
    for (size_t y = 0; y < yDimension; y++)
    {
      for (size_t x = 0; x < xDimension; x++)
      {
        delete forest[getIndex(x,y)];
      }
    }
    delete forest;
  }
  if (forestLocks != NULL)
  {
    for (size_t y = 0; y < yDimension; y++)
    {
      for (size_t x = 0; x < xDimension; x++)
      {
        omp_destroy_lock(forestLocks + (y * xDimension + x));
      }
    }
    delete forestLocks;
  }
}

template <typename T>
size_t Forest<T>::getTileX(float x)
{
  if (x < lowerLeft.x())
  {
    return 0;
  }
  if (upperRight.x() <= x)
  {
    return xDimension - 1;
  }
  return (x - lowerLeft.x()) / dimensionLength;
}

template <typename T>
size_t Forest<T>::getTileY(float y)
{
  if (y < lowerLeft.y())
  {
    return 0;
  }
  if (upperRight.y() <= y)
  {
    return yDimension - 1;
  }
  return (y - lowerLeft.y()) / dimensionLength;
}

template <typename T>
size_t Forest<T>::getIndex(size_t x, size_t y)
{
  return y * xDimension + x;
}

template <typename T>
void Forest<T>::insert(Point loc, T idx)
{
  forest[getIndex(getTileX(loc.x()),getTileY(loc.y()))]->insert(loc, idx);
}

template <typename T>
AreaLock Forest<T>::lockArea(Point loc, float radius)
{
  std::vector<omp_lock_t*> locks;

  size_t minXTile = getTileX(loc.x()-radius);
  size_t minYTile = getTileY(loc.y()-radius);
  size_t maxXTile = getTileX(loc.x()+radius);
  size_t maxYTile = getTileY(loc.y()+radius);
  for (size_t yTile = minYTile; yTile <= maxYTile; yTile++)
  {
    for (size_t xTile = minXTile; xTile <= maxXTile; xTile++)
    {
      // treesToQuery.push_back(getIndex(xTile, yTile));
      locks.push_back(&forestLocks[getIndex(xTile,yTile)]);
    }
  }

  return AreaLock(locks);
}

template <typename T>
std::vector<T> Forest<T>::searchRange(Point loc, float radius)
{
  std::vector<size_t> treesToQuery;

  size_t minXTile = getTileX(loc.x()-radius);
  size_t minYTile = getTileY(loc.y()-radius);
  size_t maxXTile = getTileX(loc.x()+radius);
  size_t maxYTile = getTileY(loc.y()+radius);
  for (size_t yTile = minYTile; yTile <= maxYTile; yTile++)
  {
    for (size_t xTile = minXTile; xTile <= maxXTile; xTile++)
    {
      treesToQuery.push_back(getIndex(xTile, yTile));
    }
  }

  std::vector<T> foundItems;
  for (size_t treeIdx : treesToQuery)
  {
    KDTree<T>* tree = forest[treeIdx];
    std::vector<T> foundSome = tree->searchRange(loc, radius);
    foundItems.insert(foundItems.end(), foundSome.begin(), foundSome.end());
  }

  return foundItems;
}

#endif