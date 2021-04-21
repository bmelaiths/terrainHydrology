#include "../kdtree.hpp"
#include "gtest/gtest.h"

namespace
{
    TEST(KDTreeTest, CreationTest)
    {
        KDTree tree;

        tree.insert(Point(3.0f, 6.0f), 0);

        std::vector<size_t> structure = tree.breadthFirstSearch();

        ASSERT_EQ(structure.size(), 1);
    }

    TEST(KDTreeTest, StructureTest) {
        KDTree tree;

        tree.insert(Point(3.0f, 6.0f),0);
        tree.insert(Point(17.0f, 15.0f),2);
        tree.insert(Point(13.0f, 15.0f),4);
        tree.insert(Point(6.0f, 12.0f),3);
        tree.insert(Point(9.0f, 1.0f),5);
        tree.insert(Point(2.0f, 7.0f),1);
        tree.insert(Point(10.0f, 19.0f),6);

        std::vector<size_t> structure = tree.breadthFirstSearch();

        for (size_t i = 0; i < 7; i++)
        {
            ASSERT_EQ(structure.at(i), i);
        }
    }

    TEST(KDTreeTest, RangeSearchTest) {
        KDTree tree;

        tree.insert(Point(3.0f, 6.0f),0);
        tree.insert(Point(17.0f, 15.0f),2);
        tree.insert(Point(13.0f, 15.0f),4);
        tree.insert(Point(6.0f, 12.0f),3);
        tree.insert(Point(9.0f, 1.0f),5);
        tree.insert(Point(2.0f, 7.0f),1);
        tree.insert(Point(10.0f, 19.0f),6);

        std::vector<size_t> searchResults = tree.rangeSearch(Point(2.0f,6.0f), 2.0f);

        ASSERT_EQ(searchResults.size(), 2);
        ASSERT_TRUE(searchResults[0] == 0 || searchResults[1] == 0);
        ASSERT_TRUE(searchResults[0] == 1 || searchResults[1] == 1);
    }

    TEST(KDTreeTest, RangeSearchTestII) {
        KDTree tree;

        tree.insert(Point(7,5),0);
        tree.insert(Point(7,3),1);
        tree.insert(Point(2,3),2);
        tree.insert(Point(7,10),3);
        tree.insert(Point(9,8),4);
        tree.insert(Point(4,8),5);
        tree.insert(Point(5,3),6);
        tree.insert(Point(8,3),7);
        tree.insert(Point(3,1),8);
        tree.insert(Point(7,9),9);
        tree.insert(Point(3,6),10);
        tree.insert(Point(2,5),11);
        tree.insert(Point(3,10),12);
        tree.insert(Point(0,4),13);
        tree.insert(Point(5,6),14);
        tree.insert(Point(1,6),15);
        tree.insert(Point(10,5),16);
        tree.insert(Point(0,2),17);

        std::vector<size_t> searchResults = tree.rangeSearch(Point(2,5), 1.5);

        ASSERT_EQ(searchResults.size(), 3);
        ASSERT_TRUE((searchResults[0] == 11) || (searchResults[1] == 11) || (searchResults[2] == 11));
        ASSERT_TRUE((searchResults[0] == 10) || (searchResults[1] == 10) || (searchResults[2] == 10));
        ASSERT_TRUE((searchResults[0] == 15) || (searchResults[1] == 15) || (searchResults[2] == 15));
    }
} // namespace