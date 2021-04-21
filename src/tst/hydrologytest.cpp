#include "gtest/gtest.h"
#include "../hydrology.hpp"

#include <vector>

namespace
{
    TEST(HydrologyTest, IDTest)
    {
        Hydrology hydrology;

        Primitive node0 = hydrology.addMouthNode  (
            Point(6.0f,11.0f), 0.0f, 0, 0
        );
        Primitive node1 = hydrology.addMouthNode  (
            Point(6.0f, 7.0f), 0.0f, 0, 0
        );
        Primitive node2 = hydrology.addMouthNode  (
            Point(4.0f, 4.0f), 0.0f, 0, 0
        );
        Primitive node3 = hydrology.addMouthNode  (
            Point(9.0f, 2.0f), 0.0f, 0, 0
        );
        Primitive node4 = hydrology.addRegularNode(
            Point(4.0f, 7.0f), 0.0f, 0, node0.id
        );
        Primitive node5 = hydrology.addRegularNode(
            Point(10.0f,9.0f), 0.0f, 0, node1.id
        );
        Primitive node6 = hydrology.addRegularNode(
            Point(7.0f, 4.0f), 0.0f, 0, node2.id
        );
        Primitive node7 = hydrology.addRegularNode(
            Point(9.0f, 6.0f), 0.0f, 0, node3.id
        );

        ASSERT_EQ(node0.id, 0);
        ASSERT_EQ(node1.id, 1);
        ASSERT_EQ(node2.id, 2);
        ASSERT_EQ(node3.id, 3);
        ASSERT_EQ(node4.id, 4);
        ASSERT_EQ(node5.id, 5);
        ASSERT_EQ(node6.id, 6);
        ASSERT_EQ(node7.id, 7);
    }

    TEST(HydrologyTest, BallPointSearch)
    {
        Hydrology hydrology;

        Primitive node0 = hydrology.addMouthNode  (
            Point(6.0f,11.0f), 0.0f, 0, 0
        );
        Primitive node1 = hydrology.addMouthNode  (
            Point(6.0f, 7.0f), 0.0f, 0, 0
        );
        Primitive node2 = hydrology.addMouthNode  (
            Point(4.0f, 4.0f), 0.0f, 0, 0
        );
        Primitive node3 = hydrology.addMouthNode  (
            Point(9.0f, 2.0f), 0.0f, 0, 0
        );
        /*Primitive node4 = */hydrology.addRegularNode(
            Point(4.0f, 7.0f), 0.0f, 0, node0.id
        );
        /*Primitive node5 = */hydrology.addRegularNode(
            Point(10.0f,9.0f), 0.0f, 0, node1.id
        );
        /*Primitive node6 = */hydrology.addRegularNode(
            Point(7.0f, 4.0f), 0.0f, 0, node2.id
        );
        /*Primitive node7 = */hydrology.addRegularNode(
            Point(9.0f, 6.0f), 0.0f, 0, node3.id
        );

        std::vector<Edge> edges = hydrology.edgesWithinRadius(
            Point(5.0f,5.0f), 3.0f
        );

        ASSERT_EQ(edges.size(), 4);
    }

    TEST(HydrologyTest, BallPointSearchII)
    {
        Hydrology hydrology;

        Primitive node0 = hydrology.addMouthNode  (
            Point(3.0f, 7.0f), 0.0f, 0, 0
        );
        /*Primitive node1 = */hydrology.addRegularNode(
            Point(2.0f, 10.0f), 0.0f, 0, node0.id
        );
        Primitive node2 = hydrology.addMouthNode  (
            Point(7.0f, 7.0f), 0.0f, 0, 0
        );
        /*Primitive node3 = */hydrology.addRegularNode(
            Point(8.0f, 10.0f), 0.0f, 0, node2.id
        );
        Primitive node4 = hydrology.addMouthNode  (
            Point(3.0f, 3.0f), 0.0f, 0, 0
        );
        /*Primitive node5 = */hydrology.addRegularNode(
            Point(0.0f, 0.0f), 0.0f, 0, node4.id
        );
        Primitive node6 = hydrology.addMouthNode  (
            Point(7.0f, 3.0f), 0.0f, 0, 0
        );
        /*Primitive node7 = */hydrology.addRegularNode(
            Point(10.0f, 0.0f), 0.0f, 0, node6.id
        );

        std::vector<Edge> edges = hydrology.edgesWithinRadius(
            Point(5.0f,5.0f), 3.0f
        );

        ASSERT_EQ(edges.size(), 4);

        edges = hydrology.edgesWithinRadius(
            Point(5.0f,5.0f), 6.0f
        );

        ASSERT_EQ(edges.size(), 8);
    }
}