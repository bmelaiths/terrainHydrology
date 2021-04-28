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
    float x() {return xLoc;};
    float y() {return yLoc;};
};

#endif