//
// Created by randy on 19-6-14.
//

#ifndef VSLAM_COMMON_INCLUDE_H
#define VSLAM_COMMON_INCLUDE_H

// define the commonly included file to avoid a long include list
// for Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;
// for Sophus
#include <sophus/se3.h>
using Sophus::SE3;
// for cv
#include <opencv2/core/core.hpp>
using cv::Mat;

// std
#include <vector>
#include <list>
#include <memory>
#include <string>
#include <iostream>
#include <set>
#include <unordered_map>
#include <map>

using namespace std;
namespace vslam {
    struct myKeypoint {
        cv::KeyPoint cvKp;
        unsigned long kpID;     // use to index the keypoint of current frame among 3D world map points

    };
}


#endif
