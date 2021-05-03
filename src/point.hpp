#ifndef POINT_H
#define POINT_H

struct Point
{
    private:
    float xLoc,yLoc;

    public:
    Point(float x, float y)
    : xLoc(x), yLoc(y)
    {};
    float x() const {return xLoc;};
    float y() const {return yLoc;};
};

#endif