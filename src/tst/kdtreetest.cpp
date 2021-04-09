#include "../kdtree.hpp"
#include "gtest/gtest.h"

namespace
{
    TEST(KDTreeTest, CreationTest)
    {
        std::vector<Point>  initPoints;
        std::vector<size_t>  initIdxes;

        initPoints.push_back(Point(3.0f, 6.0f));
        initIdxes.push_back(0);

        KDTree tree(initPoints, initIdxes);

        std::vector<size_t> structure = tree.breadthFirstSearch();

        ASSERT_EQ(structure.size(), 1);
    }

    TEST(KDTreeTest, StructureTest) {
        std::vector<Point>  initPoints;
        std::vector<size_t>  initIdxes;

        initPoints.push_back(Point(3.0f, 6.0f));
        initIdxes.push_back(0);
        initPoints.push_back(Point(17.0f, 15.0f));
        initIdxes.push_back(2);
        initPoints.push_back(Point(13.0f, 15.0f));
        initIdxes.push_back(4);
        initPoints.push_back(Point(6.0f, 12.0f));
        initIdxes.push_back(3);
        initPoints.push_back(Point(9.0f, 1.0f));
        initIdxes.push_back(5);
        initPoints.push_back(Point(2.0f, 7.0f));
        initIdxes.push_back(1);
        initPoints.push_back(Point(10.0f, 19.0f));
        initIdxes.push_back(6);

        KDTree tree(initPoints, initIdxes);
        std::vector<size_t> structure = tree.breadthFirstSearch();

        for (size_t i = 0; i < 7; i++)
        {
            ASSERT_EQ(structure.at(i), i);
        }
    }

    TEST(KDTreeTest, RangeSearchTest) {
        std::vector<Point>  initPoints;
        std::vector<size_t>  initIdxes;

        initPoints.push_back(Point(3.0f, 6.0f));
        initIdxes.push_back(0);
        initPoints.push_back(Point(17.0f, 15.0f));
        initIdxes.push_back(2);
        initPoints.push_back(Point(13.0f, 15.0f));
        initIdxes.push_back(4);
        initPoints.push_back(Point(6.0f, 12.0f));
        initIdxes.push_back(3);
        initPoints.push_back(Point(9.0f, 1.0f));
        initIdxes.push_back(5);
        initPoints.push_back(Point(2.0f, 7.0f));
        initIdxes.push_back(1);
        initPoints.push_back(Point(10.0f, 19.0f));
        initIdxes.push_back(6);

        KDTree tree(initPoints, initIdxes);
        std::vector<size_t> searchResults = tree.rangeSearch(Point(2.0f,6.0f), 2.0f);

        ASSERT_EQ(searchResults.size(), 2);
        ASSERT_TRUE(searchResults[0] == 0 || searchResults[1] == 0);
        ASSERT_TRUE(searchResults[0] == 1 || searchResults[1] == 1);
    }

    TEST(KDTreeTest, RangeSearchTestII) {
        std::vector<Point>  initPoints;
        std::vector<size_t>  initIdxes;

        initPoints.push_back(Point(7,5));
        initIdxes.push_back(0);
        initPoints.push_back(Point(7,3));
        initIdxes.push_back(1);
        initPoints.push_back(Point(2,3));
        initIdxes.push_back(2);
        initPoints.push_back(Point(7,10));
        initIdxes.push_back(3);
        initPoints.push_back(Point(9,8));
        initIdxes.push_back(4);
        initPoints.push_back(Point(4,8));
        initIdxes.push_back(5);
        initPoints.push_back(Point(5,3));
        initIdxes.push_back(6);
        initPoints.push_back(Point(8,3));
        initIdxes.push_back(7);
        initPoints.push_back(Point(3,1));
        initIdxes.push_back(8);
        initPoints.push_back(Point(7,9));
        initIdxes.push_back(9);
        initPoints.push_back(Point(3,6));
        initIdxes.push_back(10);
        initPoints.push_back(Point(2,5));
        initIdxes.push_back(11);
        initPoints.push_back(Point(3,10));
        initIdxes.push_back(12);
        initPoints.push_back(Point(0,4));
        initIdxes.push_back(13);
        initPoints.push_back(Point(5,6));
        initIdxes.push_back(14);
        initPoints.push_back(Point(1,6));
        initIdxes.push_back(15);
        initPoints.push_back(Point(10,5));
        initIdxes.push_back(16);
        initPoints.push_back(Point(0,2));
        initIdxes.push_back(17);

        KDTree tree(initPoints, initIdxes);
        std::vector<size_t> searchResults = tree.rangeSearch(Point(2,5), 1.5);

        ASSERT_EQ(searchResults.size(), 3);
        ASSERT_TRUE((searchResults[0] == 11) || (searchResults[1] == 11) || (searchResults[2] == 11));
        ASSERT_TRUE((searchResults[0] == 10) || (searchResults[1] == 10) || (searchResults[2] == 10));
        ASSERT_TRUE((searchResults[0] == 15) || (searchResults[1] == 15) || (searchResults[2] == 15));
    }
} // namespace