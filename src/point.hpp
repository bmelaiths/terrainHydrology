#ifndef POINT_H
#define POINT_H

/**
 * @brief Encapsulates a 2D point
 * 
 */
struct Point
{
    private:
    float xLoc,yLoc;

    public:
    /**
     * @brief Construct a new Point object
     * 
     * @param x 
     * @param y 
     */
    Point(float x, float y)
    : xLoc(x), yLoc(y)
    {};
    /**
     * @brief Gets the x value
     * 
     * @return float 
     */
    float x() const {return xLoc;};
    /**
     * @brief Gets the y value
     * 
     * @return float 
     */
    float y() const {return yLoc;};
};

#endif