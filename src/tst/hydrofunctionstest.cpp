#include "gtest/gtest.h"

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "../hydrologyFunctions.hpp"

namespace
{
    TEST(HydrologyFunctionsTests, SelectNodeTest)
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

        ASSERT_EQ(selected.id, (size_t)3);
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child2 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        EXPECT_TRUE(isAcceptablePosition(Point(1540*params.resolution,1415*params.resolution), 0 ,params));
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1560*params.resolution,1330*params.resolution), 0,params)); //not on land
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1540*params.resolution,1410*params.resolution),0,params)); // too close to a node
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        EXPECT_FALSE(isAcceptablePosition(Point(1560*params.resolution,1375*params.resolution),0,params)); //too close to an edge
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child0 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        EXPECT_FALSE(isAcceptablePosition(Point(250*params.resolution,134*params.resolution),0,params)); //too close to the seeee
    }
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit

        const size_t contourIndex = 175;
        Primitive mouth = params.hydrology.addMouthNode(
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

        HydrologyParameters params;
        params.contour = contours[0];
        params.resolution = 2.0; //space units / map unit
        params.edgeLength = 40.0; //space units
        params.eta = 1.5;
        params.sigma = 1.1;
        params.riverAngleDev = 0.05;
        params.maxTries = 15;

        Primitive mouth = params.hydrology.addMouthNode(
            Point(1530*params.resolution,1340*params.resolution), 0.0, 0, 0
        );
        Primitive child0 = params.hydrology.addRegularNode(
            Point(1520*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        Primitive child1 = params.hydrology.addRegularNode(
            Point(1540*params.resolution,1360*params.resolution), 0.0, 0, mouth.id
        );
        /*Primitive child2 = */params.hydrology.addRegularNode(
            Point(1540*params.resolution,1390*params.resolution), 0.0, 0, child1.id
        );

        Point newLoc = pickNewNodeLoc(child0, params);
        ASSERT_TRUE(isAcceptablePosition(newLoc, child0.id, params));
    }
}