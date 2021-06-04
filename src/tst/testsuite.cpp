#include "gtest/gtest.h"

#include <vector>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "../hydrologyFunctions.hpp"
#include "../kdtree.hpp"
#include "../hydrology.hpp"
#include "../forest.hpp"
#include "../terrainElevation.hpp"

namespace
{
    TEST(KDTreeTest, CreationTest)
    {
        KDTree<size_t> tree;

        tree.insert(Point(3.0f, 6.0f), 0);

        std::vector<size_t> structure = tree.breadthFirstSearch();

        ASSERT_EQ(structure.size(), 1);
    }

    TEST(KDTreeTest, StructureTest) {
        KDTree<size_t> tree;

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
        KDTree<size_t> tree;

        tree.insert(Point(3.0f, 6.0f),0);
        tree.insert(Point(17.0f, 15.0f),2);
        tree.insert(Point(13.0f, 15.0f),4);
        tree.insert(Point(6.0f, 12.0f),3);
        tree.insert(Point(9.0f, 1.0f),5);
        tree.insert(Point(2.0f, 7.0f),1);
        tree.insert(Point(10.0f, 19.0f),6);

        std::vector<size_t> searchResults = tree.searchRange(Point(2.0f,6.0f), 2.0f);

        ASSERT_EQ(searchResults.size(), 2);
        ASSERT_TRUE(searchResults[0] == 0 || searchResults[1] == 0);
        ASSERT_TRUE(searchResults[0] == 1 || searchResults[1] == 1);
    }

    TEST(KDTreeTest, RangeSearchTestII) {
        KDTree<size_t> tree;

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

        std::vector<size_t> searchResults = tree.searchRange(Point(2,5), 1.5);

        ASSERT_EQ(searchResults.size(), 3);
        ASSERT_TRUE((searchResults[0] == 11) || (searchResults[1] == 11) || (searchResults[2] == 11));
        ASSERT_TRUE((searchResults[0] == 10) || (searchResults[1] == 10) || (searchResults[2] == 10));
        ASSERT_TRUE((searchResults[0] == 15) || (searchResults[1] == 15) || (searchResults[2] == 15));
    }
    TEST(KDTreeTest, RangeSearchReconstructTest) {
        KDTree<size_t> tree;

        tree.insert(Point(3.0f, 6.0f),0);
        tree.insert(Point(17.0f, 15.0f),2);
        tree.insert(Point(13.0f, 15.0f),4);
        tree.insert(Point(6.0f, 12.0f),3);
        tree.insert(Point(9.0f, 1.0f),5);
        tree.insert(Point(2.0f, 7.0f),1);
        tree.insert(Point(10.0f, 19.0f),6);

        tree.reconstruct();

        std::vector<size_t> searchResults = tree.searchRange(Point(2.0f,6.0f), 2.0f);

        ASSERT_EQ(searchResults.size(), 2);
        ASSERT_TRUE(searchResults[0] == 0 || searchResults[1] == 0);
        ASSERT_TRUE(searchResults[0] == 1 || searchResults[1] == 1);
    }
    TEST(KDTreeTest, RangeSearchReconstructTestII) {
        KDTree<size_t> tree;

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

        tree.reconstruct();

        std::vector<size_t> searchResults = tree.searchRange(Point(2,5), 1.5);

        ASSERT_EQ(searchResults.size(), 3);
        ASSERT_TRUE((searchResults[0] == 11) || (searchResults[1] == 11) || (searchResults[2] == 11));
        ASSERT_TRUE((searchResults[0] == 10) || (searchResults[1] == 10) || (searchResults[2] == 10));
        ASSERT_TRUE((searchResults[0] == 15) || (searchResults[1] == 15) || (searchResults[2] == 15));
    }
    TEST(ForestTest, CreationTest)
    {
        Forest<size_t> trees(Point(-2.0,-2.0), Point(4.0, 7.0), 2.0);

        trees.insert(Point(2.0, 6.0), 0);

        std::vector<size_t> getPointBack = trees.searchRange(Point(3.0, 6.0), 1.0);

        ASSERT_EQ(getPointBack.size(), 1);
    }
    TEST(ForestTest, RangeSearchTest) {
        Forest<size_t> trees(Point(0,-5), Point(20,25), 5.0);

        trees.insert(Point(3.0f, 6.0f),0);
        trees.insert(Point(17.0f, 15.0f),2);
        trees.insert(Point(13.0f, 15.0f),4);
        trees.insert(Point(6.0f, 12.0f),3);
        trees.insert(Point(9.0f, 1.0f),5);
        trees.insert(Point(2.0f, 7.0f),1);
        trees.insert(Point(10.0f, 19.0f),6);

        std::vector<size_t> searchResults = trees.searchRange(Point(2.0f,6.0f), 2.0f);

        ASSERT_EQ(searchResults.size(), 2);
        ASSERT_TRUE(searchResults[0] == 0 || searchResults[1] == 0);
        ASSERT_TRUE(searchResults[0] == 1 || searchResults[1] == 1);
    }
    TEST(ForestTest, RangeSearchTestII) {
        Forest<size_t> trees(Point(0,0), Point(10,10), 5);

        trees.insert(Point(7,5),0);
        trees.insert(Point(7,3),1);
        trees.insert(Point(2,3),2);
        trees.insert(Point(7,10),3);
        trees.insert(Point(9,8),4);
        trees.insert(Point(4,8),5);
        trees.insert(Point(5,3),6);
        trees.insert(Point(8,3),7);
        trees.insert(Point(3,1),8);
        trees.insert(Point(7,9),9);
        trees.insert(Point(3,6),10);
        trees.insert(Point(2,5),11);
        trees.insert(Point(3,10),12);
        trees.insert(Point(0,4),13);
        trees.insert(Point(5,6),14);
        trees.insert(Point(1,6),15);
        trees.insert(Point(10,5),16);
        trees.insert(Point(0,2),17);

        std::vector<size_t> searchResults = trees.searchRange(Point(2,5), 1.5);

        ASSERT_EQ(searchResults.size(), 3);
        ASSERT_TRUE((searchResults[0] == 11) || (searchResults[1] == 11) || (searchResults[2] == 11));
        ASSERT_TRUE((searchResults[0] == 10) || (searchResults[1] == 10) || (searchResults[2] == 10));
        ASSERT_TRUE((searchResults[0] == 15) || (searchResults[1] == 15) || (searchResults[2] == 15));
    }
    TEST(ForestTest, RangeSearchTestIII)
    {
        Forest<size_t> trees(Point(0,0), Point(15,20), 2);

        trees.insert(Point(6.0f,11.0f),0);
        trees.insert(Point(6.0f, 7.0f),0);
        trees.insert(Point(4.0f, 4.0f),0);
        trees.insert(Point(9.0f, 2.0f),0);
        trees.insert(Point(4.0f, 7.0f),0);
        trees.insert(Point(10.0f,9.0f),0);
        trees.insert(Point(7.0f, 4.0f),0);
        trees.insert(Point(9.0f, 6.0f),0);

        std::vector<size_t> searchResults = trees.searchRange(Point(5.0f,5.0f), 3.0f);

        ASSERT_EQ(searchResults.size(), 4);
    }
    TEST(HydrologyTest, IDTest)
    {
        Hydrology hydrology(Point(0,0), Point(15,20), 1.0);

        Primitive node0 = *hydrology.addMouthNode  (
            Point(6.0f,11.0f), 0.0f, 0, 0
        );
        Primitive node1 = *hydrology.addMouthNode  (
            Point(6.0f, 7.0f), 0.0f, 0, 0
        );
        Primitive node2 = *hydrology.addMouthNode  (
            Point(4.0f, 4.0f), 0.0f, 0, 0
        );
        Primitive node3 = *hydrology.addMouthNode  (
            Point(9.0f, 2.0f), 0.0f, 0, 0
        );
        Primitive node4 = *hydrology.addRegularNode(
            Point(4.0f, 7.0f), 0.0f, 0, node0.getID()
        );
        Primitive node5 = *hydrology.addRegularNode(
            Point(10.0f,9.0f), 0.0f, 0, node1.getID()
        );
        Primitive node6 = *hydrology.addRegularNode(
            Point(7.0f, 4.0f), 0.0f, 0, node2.getID()
        );
        Primitive node7 = *hydrology.addRegularNode(
            Point(9.0f, 6.0f), 0.0f, 0, node3.getID()
        );

        ASSERT_EQ(node0.getID(), (size_t)0);
        ASSERT_EQ(node1.getID(), (size_t)1);
        ASSERT_EQ(node2.getID(), (size_t)2);
        ASSERT_EQ(node3.getID(), (size_t)3);
        ASSERT_EQ(node4.getID(), (size_t)4);
        ASSERT_EQ(node5.getID(), (size_t)5);
        ASSERT_EQ(node6.getID(), (size_t)6);
        ASSERT_EQ(node7.getID(), (size_t)7);
    }

    TEST(HydrologyTest, BallPointSearch)
    {
        Hydrology hydrology(Point(0,0), Point(15,20), 1.0);

        Primitive node0 = *hydrology.addMouthNode  (
            Point(6.0f,11.0f), 0.0f, 0, 0
        );
        Primitive node1 = *hydrology.addMouthNode  (
            Point(6.0f, 7.0f), 0.0f, 0, 0
        );
        Primitive node2 = *hydrology.addMouthNode  (
            Point(4.0f, 4.0f), 0.0f, 0, 0
        );
        Primitive node3 = *hydrology.addMouthNode  (
            Point(9.0f, 2.0f), 0.0f, 0, 0
        );
        /*Primitive node4 = */hydrology.addRegularNode(
            Point(4.0f, 7.0f), 0.0f, 0, node0.getID()
        );
        /*Primitive node5 = */hydrology.addRegularNode(
            Point(10.0f,9.0f), 0.0f, 0, node1.getID()
        );
        /*Primitive node6 = */hydrology.addRegularNode(
            Point(7.0f, 4.0f), 0.0f, 0, node2.getID()
        );
        /*Primitive node7 = */hydrology.addRegularNode(
            Point(9.0f, 6.0f), 0.0f, 0, node3.getID()
        );

        std::vector<Edge> edges = hydrology.queryArea(
            Point(5.0f,5.0f), 3.0f
        );

        ASSERT_EQ(edges.size(), 4);
    }

    TEST(HydrologyTest, BallPointSearchII)
    {
        Hydrology hydrology(Point(0,0), Point(15,20), 1.0);

        Primitive node0 = *hydrology.addMouthNode  (
            Point(3.0f, 7.0f), 0.0f, 0, 0
        );
        /*Primitive node1 = */hydrology.addRegularNode(
            Point(2.0f, 10.0f), 0.0f, 0, node0.getID()
        );
        Primitive node2 = *hydrology.addMouthNode  (
            Point(7.0f, 7.0f), 0.0f, 0, 0
        );
        /*Primitive node3 = */hydrology.addRegularNode(
            Point(8.0f, 10.0f), 0.0f, 0, node2.getID()
        );
        Primitive node4 = *hydrology.addMouthNode  (
            Point(3.0f, 3.0f), 0.0f, 0, 0
        );
        /*Primitive node5 = */hydrology.addRegularNode(
            Point(0.0f, 0.0f), 0.0f, 0, node4.getID()
        );
        Primitive node6 = *hydrology.addMouthNode  (
            Point(7.0f, 3.0f), 0.0f, 0, 0
        );
        /*Primitive node7 = */hydrology.addRegularNode(
            Point(10.0f, 0.0f), 0.0f, 0, node6.getID()
        );

        std::vector<Edge> edges = hydrology.queryArea(
            Point(5.0f,5.0f), 3.0f
        );

        ASSERT_EQ(edges.size(), 4);

        edges = hydrology.queryArea(
            Point(5.0f,5.0f), 6.0f
        );

        ASSERT_EQ(edges.size(), 8);
    }

    TEST(HydrologyFunctionsTests, SelectNodeTest)
    {
        HydrologyParameters testParams(Point(-1,-1), Point(1,1));
        testParams.zeta = 14.0f;
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 4.0, 1, 0
            )
        );
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 6.0, 2, 0
            )
        );
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 14.0, 3, 0
            )
        );
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 8.0, 3, 0
            )
        );
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 24.0, 1, 0
            )
        );
        testParams.candidates.push_back(
            testParams.hydrology.addMouthNode(
                Point(0,0), 23.0, 4, 0
            )
        );

        Primitive selected = selectNode(testParams);

        ASSERT_EQ(selected.getID(), (size_t)3);
    }
    TEST(HydrologyTest, NodeAdditionTest)
    {
        const float resolution = 13.5;
        HydrologyParameters params(Point(1500*resolution,1300*resolution), Point(1550*resolution,1400*resolution));
        params.resolution = resolution;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        ASSERT_EQ(params.hydrology.numNodes(), 4);
    }
    TEST(HydrologyFunctionsTests, IsAcceptablePositionAcceptableTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        // cv::Mat raw_dist( src.size(), CV_32F );
        // for( int i = 0; i < src.rows; i++ )
        // {
        //     for( int j = 0; j < src.cols; j++ )
        //     {
        //         raw_dist.at<float>(i,j) = (float)cv::pointPolygonTest( contours[0], cv::Point2f((float)j, (float)i), true );
        //     }
        // }
        // double minVal, maxVal;
        // cv::Point maxDistPt; // inscribed circle center
        // cv::minMaxLoc(raw_dist, &minVal, &maxVal, NULL, &maxDistPt);
        // minVal = abs(minVal);
        // maxVal = abs(maxVal);
        // cv::Mat drawing = cv::Mat::zeros( src.size(), CV_8UC3 );
        // for( int i = 0; i < src.rows; i++ )
        // {
        //     for( int j = 0; j < src.cols; j++ )
        //     {
        //         if( raw_dist.at<float>(i,j) < 0 )
        //         {
        //             drawing.at<cv::Vec3b>(i,j)[0] = (uchar)(255 - abs(raw_dist.at<float>(i,j)) * 255 / minVal);
        //         }
        //         else if( raw_dist.at<float>(i,j) > 0 )
        //         {
        //             drawing.at<cv::Vec3b>(i,j)[2] = (uchar)(255 - raw_dist.at<float>(i,j) * 255 / maxVal);
        //         }
        //         else
        //         {
        //             drawing.at<cv::Vec3b>(i,j)[0] = 255;
        //             drawing.at<cv::Vec3b>(i,j)[1] = 255;
        //             drawing.at<cv::Vec3b>(i,j)[2] = 255;
        //         }
        //     }
        // }
        // cv::circle(drawing, maxDistPt, (int)maxVal, cv::Scalar(255,255,255));
        // cv::imshow( "Source", src );
        // cv::imshow( "Distance and inscribed circle", drawing );
        // cv::waitKey(30000);

        HydrologyParameters params(Point(1500*2.0,1300*2.0),Point(1550*2.0,1400*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 1.1;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child2 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        EXPECT_TRUE(isAcceptablePosition(Point(1540*params.resolution,1415*params.resolution), 2 * params.edgeLength, 0 ,params));
    }
    TEST(HydrologyFunctionsTests, IsAcceptablePositionNotOnLandTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(1500*2.0,1300*2.0),Point(1550*2.0,1400*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 1.1;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1560*params.resolution,1330*params.resolution), 2 * params.edgeLength, 0,params)); //not on land
    }
    TEST(HydrologyFunctionsTests, IsAcceptablePositionTooCloseToNodeTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(1500*2.0,1300*2.0),Point(1550*2.0,1400*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 1.1;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1540*params.resolution,1410*params.resolution),2 * params.edgeLength,0,params)); // too close to a node
    }
    TEST(HydrologyFunctionsTests, IsAcceptablePositionTooCloseToEdgeTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(1500*2.0,1300*2.0),Point(1550*2.0,1400*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 1.1;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1560*params.resolution,1375*params.resolution),2 * params.edgeLength,0,params)); //too close to an edge
    }
    TEST(HydrologyFunctionsTests, IsAcceptablePositionTooCloseToSeeeTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(1500*2.0,1300*2.0),Point(1550*2.0,1400*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 1.1;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        EXPECT_FALSE(isAcceptablePosition(Point(250*params.resolution,134*params.resolution),2 * params.edgeLength,0,params)); //too close to the seeee
    }
    // TEST(HydrologyFunctionsTests, IsAcceptablePositionRealDataTestI)
    // {
    //     FILE *input = fopen("./binaryFile", "rb");

    //     if (input == NULL)
    //     {
    //         printf("Unable to open file\n");
    //         exit(1);
    //     }

    //     HydrologyParameters params = readParamsFromStream(input);

    //     Primitive testNode(10, 6, Point(105000.0, 70100.0), 0.0, 10);

    //     ASSERT_FALSE(isAcceptablePosition(testNode.loc, 6, params));
    // }
    // TEST(HydrologyFunctionsTests, IsAcceptablePositionRealDataTestII)
    // {
    //     FILE *input = fopen("./binaryFile", "rb");

    //     if (input == NULL)
    //     {
    //         printf("Unable to open file\n");
    //         exit(1);
    //     }

    //     HydrologyParameters params = readParamsFromStream(input);

    //     Point testLoc(
    //         params.hydrology.indexedNodes[6].loc.x-params.edgeLength,
    //         params.hydrology.indexedNodes[6].loc.y
    //     );

    //     ASSERT_TRUE(isAcceptablePosition(testLoc, 6, params));
    // }
    TEST(HydrologyFunctionsTests, CoastNormalTest) {
        /*
           NOTE: This test depends on cv::findContours always generating
           the contour in roughly the same order each time. Thus, this
           test could be a little flaky
        */

        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(500*2.0,1000*2.0),Point(3500*2.0,3000*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit

        const size_t contourIndex = 175;
        Primitive mouth = *params.hydrology.addMouthNode(
            Point(params.contour[contourIndex].x,params.contour[contourIndex].y),
            0.0f, 0, contourIndex
        );

        float angle = coastNormal(mouth,params);
        ASSERT_LT(abs(angle - M_PI_4f32), 0.1571);
    }
    TEST(HydrologyFunctionsTests, PickNewNodeLocTest)
    {
        /* Create a shoreline */
        const int r = 1000;
        cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
        std::vector<cv::Point2f> vert(6);
        vert[0] = cv::Point( 1500, 1340 );
        vert[1] = cv::Point( 1000, 2000 );
        vert[2] = cv::Point( 1500, 2860 );
        vert[3] = cv::Point( 2500, 2860 );
        vert[4] = cv::Point( 3000, 2000 );
        vert[5] = cv::Point( 2500, 1340 );
        for( int i = 0; i < 6; i++ )
        {
            cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
        }
        std::vector<std::vector<cv::Point> > contours;
        cv::findContours( 
            src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
        );

        HydrologyParameters params(Point(500*2.0,1000*2.0),Point(3500*2.0,3000*2.0));
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 0.95;
        params.sigma = 0.95;
        params.riverAngleDev = 0.05;
        params.maxTries = 15;

        Primitive mouth = *params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        Primitive child0 = *params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        Primitive child1 = *params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.getID()
        );
        /*Primitive child2 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.getID()
        );

        LockedPoint newLoc = pickNewNodeLoc(child0, params);
        ASSERT_TRUE(isAcceptablePosition(newLoc.getLoc(),2 * params.edgeLength, child0.getID(), params));
    }
    TEST(PrimitiveTests, ToBinaryIDTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(0.0,0.0), 0.0, 0, 0
        );
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(1);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        uint64_t nodeID;
        memcpy(&nodeID, buffer, sizeof(uint64_t));
        ASSERT_EQ((uint64_t)1, be64toh(nodeID));

        delete buffer;
    }
    TEST(PrimitiveTests, ToBinaryParentTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(0.0,0.0), 0.0, 0, 0
        );
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(1);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        uint64_t parent;
        memcpy(&parent, buffer + sizeof(uint64_t), sizeof(uint64_t));
        ASSERT_EQ((uint64_t)1, be64toh(parent));

        delete buffer;
    }
    TEST(PrimitiveTests, ToChildrenTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(0);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        uint8_t numChildren;
        memcpy(&numChildren, buffer + sizeof(uint64_t) * 3, sizeof(uint8_t));
        ASSERT_EQ(2, numChildren);

        uint64_t childID;
        for (size_t child = 0; child < numChildren; child++)
        {
            memcpy
            (
                &childID,
                buffer + sizeof(uint64_t)*3 + sizeof(uint8_t) + sizeof(uint64_t)*child,
                sizeof(uint64_t)
            );
            ASSERT_EQ(child+1, be64toh(childID));
        }

        delete buffer;
    }
    TEST(PrimitiveTests, ToBinaryLocXTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(0);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        float locX;
        memcpy(
            &locX,
            buffer
                + sizeof(uint64_t)*3
                + sizeof(uint8_t)
                + sizeof(uint64_t)*2,
            sizeof(float)
        );
        locX = float_swap(locX);
        ASSERT_LT(abs(3.14-locX), 0.001);

        delete buffer;
    }
    TEST(PrimitiveTests, ToBinaryLocYTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(0);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        float locY;
        memcpy(
            &locY,
            buffer
                + sizeof(uint64_t)*3
                + sizeof(uint8_t)
                + sizeof(uint64_t)*2
                + sizeof(float),
            sizeof(float)
        );
        locY = float_swap(locY);
        ASSERT_LT(abs(5.2-locY), 0.001);

        delete buffer;
    }
    TEST(PrimitiveTests, ToBinaryElevationTest)
    {
        Hydrology hydrology(Point(-1,-1),Point(4,6), 1.0);
        hydrology.addMouthNode(
            Point(3.14,5.2), 12.1,5,10
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );
        hydrology.addRegularNode(
            Point(0,0), 0.0, 0, 0
        );

        Primitive node = hydrology.getNode(0);

        uint8_t *buffer = new uint8_t[node.binarySize()];

        node.toBinary(buffer);

        float elevation;
        memcpy(
            &elevation,
            buffer
                + sizeof(uint64_t)*3
                + sizeof(uint8_t)
                + sizeof(uint64_t)*2
                + sizeof(float)*2,
            sizeof(float)
        );
        elevation = float_swap(elevation);
        ASSERT_LT(abs(12.1-elevation), 0.001);

        delete buffer;
    }
    // TEST(RasterTests, ReadRasterTest)
    // {
    //     FILE *input = fopen("./binaryFile", "rb");

    //     if (input == NULL)
    //     {
    //         printf("Unable to open file\n");
    //         exit(1);
    //     }

    //     HydrologyParameters params = readParamsFromStream(input);

    //     ASSERT_EQ(params.riverSlope.get(1095*100,100*100),78);
    //     ASSERT_EQ(params.riverSlope.get(682*100,1206*100),78);
    //     ASSERT_EQ(params.riverSlope.get(984*100,15*100),34);
    //     ASSERT_EQ(params.riverSlope.get(875*100,924*100),15);
    // }
    TEST(TerrainTests, DistanceTest)
    {
        Point p0(0,0);
        Point p1(1,0);

        float dist = distance(p0, p1);

        ASSERT_TRUE(abs(dist-1) < 0.01);
    }
    TEST(TerrainTests, DistanceTestI)
    {
        Q q0(Point(97360.9219,30977.2168),1239.06372,1084,std::vector<size_t>());
        Q q1(Point(97604.7656,31752.8867),1200.42749,1087,std::vector<size_t>());
        T t(Point(96373.5156,31288.375),800);

        float dist = distance(t.getLoc(),q0.getPosition());

        ASSERT_GT(dist, 1035);
        ASSERT_LT(dist, 1036);
    }
    TEST(TerrainTests, PointSegmentTest)
    {
        Point p0(0,0);
        Point p1(1,1);

        Point tst(0.75,0.25);

        endpointAndDistance ead = point_segment_distance(tst,p0,p1);

        ASSERT_TRUE(abs(ead.dist-0.353553) < 0.01);
    }
    TEST(TerrainTests, PointSegmentTestII)
    {
        Point p0(73527,32541);
        Point p1(73843,34327);
        Point p2(73833,34339);

        Point tst(73578,33562);

        endpointAndDistance eadI = point_segment_distance(tst,p0,p1);
        endpointAndDistance eadII = point_segment_distance(tst,p1,p2);

        ASSERT_TRUE(eadI.dist < eadII.dist);
    }
    TEST(TerrainTest, LERPRidgeTest)
    {
        Q q0(Point(97360.9219,30977.2168),1239.06372,1084,std::vector<size_t>());
        Q q1(Point(97604.7656,31752.8867),1200.42749,1087,std::vector<size_t>());
        T t(Point(96373.5156,31288.375),800);
        float dist = 1035.27661;

        float lerpedElev = lerpRidge(&q0, &q1, t, dist);

        ASSERT_GT(lerpedElev, q1.getElevation() - 0.0001);
        ASSERT_LT(lerpedElev, q0.getElevation() + 0.0001);
    }
}