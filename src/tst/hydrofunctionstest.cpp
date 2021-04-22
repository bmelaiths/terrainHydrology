#include "gtest/gtest.h"

#include "../hydrologyFunctions.hpp"

namespace
{
    TEST(HydrologyFunctionTests, SelectNodeTest)
    {
        HydrologyParameters testParams;
        testParams.zeta = 14.0f;
        testParams.candidates.push_back(
            Primitive(0, Point(0,0), 4.0, 1, 0)
        );
        testParams.candidates.push_back(
            Primitive(1, Point(0,0), 6.0, 2, 0)
        );
        testParams.candidates.push_back(
            Primitive(2, Point(0,0), 14.0, 3, 0)
        );
        testParams.candidates.push_back(
            Primitive(3, Point(0,0), 8.0, 3, 0)
        );
        testParams.candidates.push_back(
            Primitive(4, Point(0,0), 24.0, 1, 0)
        );
        testParams.candidates.push_back(
            Primitive(5, Point(0,0), 23.0, 4, 0)
        );

        Primitive selected = selectNode(testParams);

        ASSERT_EQ(selected.id, 3);
    }
}