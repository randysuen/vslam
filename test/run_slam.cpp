#include <fstream>
#include <boost/timer.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "config.h"
#include "bundleadjustment.h"

using namespace vslam;


void createPoses (int nPoses, vector<BAPose::Ptr> &truePoses, vector<BAPose::Ptr> &noisePoses);
void createPoints (int nPoints, int width, int height, Camera::Ptr &camera, vector<BAPose::Ptr> &truePoses,
                    vector<BAPoint::Ptr> &truePoints, vector<BAPoint::Ptr> &noisePoints);

int main ( int argc, char** argv )
{
    if ( argc != 2 )
    {
        cout<<"usage: run_vo parameter_file"<<endl;
        return 1;
    }

    Config::setParameterFile(argv[1]);
    Camera::Ptr camera(new Camera);

    BAOptimizer optimizer;
    optimizer.setCamera(camera);
    optimizer.setType(BA_FULL);
    optimizer.setVerbose(true);

    const int nPoses = 6;
    const int nPoints = 1000;
    const int width = 640;
    const int height = 480;

    // create pose (R and t) groundtruth, add noise to groundtruth to create pose estimation
    vector<BAPose::Ptr> truePoses, noisePoses;

    cout << "create groundtruth and estimiation for poses..." << endl;
    createPoses(nPoses,truePoses, noisePoses);
/*    for (int i = 0; i < truePoses.size(); i++) {
        cout << "true pose " << truePoses[i]->id_ << " is " << endl << truePoses[i]->pose_.matrix() << endl;
    }
    */

    // create 3d point groundtruth, add noise to generate estimated point
    vector<BAPoint::Ptr> truePoints, noisePoints;
    cout << "create groundtruth and estimiation for points..." << endl;
    createPoints(nPoints, width, height, camera, truePoses, truePoints, noisePoints);
   // cout << "number of point is " << truePoints.size() << endl;
   // cout << "points created done..." << endl;
   /* for (int i = 0; i < truePoints.size(); i++) {
        cout << "true point " << truePoints[i]->id_ << " is " << endl << truePoints[i]->point_ << endl;
    }*/

    /*********************** full BA test *******************************/
    int k = 0;
    for (int i = 0; i < noisePoses.size(); i++) {
        BAPose::Ptr &pose = noisePoses[i];
        //BAPose::Ptr &pose = truePoses[i];
        optimizer.addPose(pose);

        // check if the jth point can be seen at the ith true pose
        for (int j = 0; j < truePoints.size(); j++) {
            BAPoint::Ptr &point3d = truePoints[j];
            BAPoint::Ptr &noisepoint3d = noisePoints[j];
            Vector2d point2d = camera->world2pixel(point3d->point_,truePoses[i]->pose_);
           // cout << "j = " << j << endl;
            if (point2d(0) < 0 || point2d(1) < 0 || point2d(0) > (width-1) || point2d(1) > (height-1)) {
                continue;
            }
            // if can be seen, generate an error term
            else {
                BAError error;
                error.setId(k);
                error.setPose(pose);
                //error.setPoint(point3d);
                error.setPoint(noisepoint3d);
                error.setMeasurement(point2d);
                optimizer.addError(std::make_shared<BAError>(error));
                k++;
               // cout << "at pose " << pose->id_ << ", point " << point3d->id_ << " can be observed at pixel " << endl
               // << camera->world2pixel(point3d->point_, pose->pose_) << endl;
            }
        }
    }

    for (int i = 0; i < truePoints.size(); i++) {
        optimizer.addPoint(noisePoints[i]);
    }

/*
    for (int i = 0; i < optimizer.poses_.size(); i++) {
        cout << "true pose " << truePoses[i]->id_ << " is " << endl << truePoses[i]->pose_.matrix() << endl;
        cout << "estimated pose " << optimizer.poses_[i]->id_ << " is " << endl << optimizer.poses_[i]->pose_.matrix() << endl;
    }
    for (int i = 0; i < optimizer.points_.size(); i++) {
        cout << "estimated point " << optimizer.points_[i]->id_ << " is " << endl << optimizer.points_[i]->point_ << endl;
        for (int j = 0; j < optimizer.errors_.size(); j++) {
            if (optimizer.errors_[j]->point3d_->id_ == optimizer.points_[i]->id_){
                cout << "it can be observed at error " << optimizer.errors_[j]->id_
                << " at pose " << optimizer.errors_[j]->pose_->id_ <<
                " at pixel " << endl << optimizer.errors_[j]->point2d_ << endl;
            }
        }
    }
 */

    cout << "run optmize..." << endl;
    // run optimize
    optimizer.optimize(20);

    return 0;
}


void createPoses (int nPoses, vector<BAPose::Ptr> &truePoses, vector<BAPose::Ptr> &noisePoses) {
    const double angle_range = 0.1;  // rad, for rotation
    const double x_range = 1.0;  // meter, for translation
    const double y_range = 1.0;
    const double z_range = 0.5;
    const double rot_step = 0.1;
    const double trans_step = 0.1;
    const double rot_noise = 0.02;
    const double trans_noise = 0.02;
    cv::RNG rng((uint64)cv::getTickCount());
    vector<BAPose::Ptr>().swap(truePoses);
    vector<BAPose::Ptr>().swap(noisePoses);
    for (int i = 0; i < nPoses; i++) {
        // rotation, using Euler angle to generate
        Matrix3d R;
        Matrix3d R_x, R_y, R_z;
        // translation
        Vector3d t;
        double x, y, z;
        // SE3 matrix
        Sophus::SE3 T;
        Sophus::SE3 Tp; // previous T
        double yaw, pitch, roll;
        // define the first pose
        if (i == 0) {
            yaw = rng.uniform(-angle_range, angle_range);
            pitch = rng.uniform(-angle_range, angle_range);
            roll = rng.uniform(-angle_range, angle_range);
            R_z << cos(yaw), -sin(yaw), 0,
                    sin(yaw), cos(yaw), 0,
                    0, 0, 1;
            R_y << cos(pitch), 0, sin(pitch),
                    0, 1, 0,
                    -sin(pitch), 0, cos(pitch);
            R_x << 1, 0, 0,
                    0, cos(roll), -sin(roll),
                    0, sin(roll), cos(roll);
            R = R_z * R_y * R_x;
            x = rng.uniform(-x_range, x_range);
            y = rng.uniform(-y_range, y_range);
            z = rng.uniform(-z_range, z_range);
            t << x, y, z;
        }
        // generate the other poses using previous pose + step
        else {
            yaw = rng.gaussian(rot_step);
            pitch = rng.gaussian(rot_step);
            roll = rng.gaussian(rot_step);
            R_z << cos(yaw), -sin(yaw), 0,
                    sin(yaw), cos(yaw), 0,
                    0, 0, 1;
            R_y << cos(pitch), 0, sin(pitch),
                    0, 1, 0,
                    -sin(pitch), 0, cos(pitch);
            R_x << 1, 0, 0,
                    0, cos(roll), -sin(roll),
                    0, sin(roll), cos(roll);
            R = R_z * R_y * R_x * Tp.rotation_matrix();
            x = rng.gaussian(trans_step);
            y = rng.gaussian(trans_step);
            z = rng.gaussian(trans_step);
            t << x, y, z;
            t += Tp.translation();
        }

        T = Sophus::SE3(R,t);
        Tp = T;
        BAPose truePose, noisePose;
        truePose.setId(i);
        truePose.setEstimate(T);
        truePoses.push_back(std::make_shared<BAPose>(truePose));
        noisePose.setId(i);
        noisePose.setEstimate(T);

        // add noise to true pose to estimated pose
        yaw = rng.gaussian(rot_noise);
        pitch = rng.gaussian(rot_noise);
        roll = rng.gaussian(rot_noise);
        R_z << cos(yaw), -sin(yaw), 0,
                sin(yaw), cos(yaw), 0,
                0, 0, 1;
        R_y << cos(pitch), 0, sin(pitch),
                0, 1, 0,
                -sin(pitch), 0, cos(pitch);
        R_x << 1, 0, 0,
                0, cos(roll), -sin(roll),
                0, sin(roll), cos(roll);
        R = R_z * R_y * R_x;
        x = rng.gaussian(trans_noise);
        y = rng.gaussian(trans_noise);
        z = rng.gaussian(trans_noise);
        t << x, y, z;
        T = SE3(R,t);
        noisePose.setEstimate(T*noisePose.pose_);
        noisePoses.push_back(std::make_shared<BAPose>(noisePose));
    }
}

void createPoints (int nPoints, int width, int height, Camera::Ptr &camera, vector<BAPose::Ptr> &truePoses,
                    vector<BAPoint::Ptr> &truePoints, vector<BAPoint::Ptr> &noisePoints) {
    const double point_noise = 20.0;
    cv::RNG rng((uint64)cv::getTickCount());
    vector<BAPoint::Ptr>().swap(truePoints);
    vector<BAPoint::Ptr>().swap(noisePoints);

    for (int i = 0; i < nPoints; i++) {
        // the pixel related to the point
        Vector2d point2d;
        point2d(0) = cvRound(rng.uniform(30.0f, (float)(width-30)));
        point2d(1) = cvRound(rng.uniform(30.0f, (float)(height-30)));
        // select one pose as baseline, where the point can be seen at this pose
        int j = cvRound(rng.uniform(0.0f, (float)(truePoses.size())-1));
        BAPose pose = *truePoses[j];

        // convert pixel to point
        Vector3d point3d = camera->pixel2world(point2d, pose.pose_);
        BAPoint truePoint, noisePoint;
        truePoint.setId(i);
        truePoint.setEstimate(point3d);
        truePoints.push_back(std::make_shared<BAPoint>(truePoint));
        noisePoint.setId(i);
        noisePoint.setEstimate(point3d);

        // add noise to true point to generate estimated point
        Vector2d delta2d;
        delta2d(0) = cvRound(rng.uniform(-point_noise, point_noise));
        delta2d(1) = cvRound(rng.uniform(-point_noise, point_noise));
        point2d += delta2d;
        point3d = camera->pixel2world(point2d, pose.pose_);
        noisePoint.setEstimate(point3d);
        noisePoints.push_back(std::make_shared<BAPoint>(noisePoint));
    }
}
