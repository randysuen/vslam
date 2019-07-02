//
// Created by randy on 19-6-18.
//

#ifndef VSLAM_BUNDLEADJUSTMENT_H
#define VSLAM_BUNDLEADJUSTMENT_H
#include "common_include.h"
#include "camera.h"
//#include "frame.h"
//#include "mappoint.h"

namespace vslam {
    class BAPose {
    public:
        typedef shared_ptr<BAPose> Ptr;
        SE3 pose_;
        int id_;
    public:
        void setId(int id);
        void setEstimate(SE3 pose);
    };

    class BAPoint{
    public:
        typedef shared_ptr<BAPoint> Ptr;
        Vector3d point_;
        int id_;
    public:
        void setId(int id);
        void setEstimate(Vector3d point);

    };

    class BAError{
    public:
        typedef shared_ptr<BAError> Ptr;
        BAPose::Ptr pose_;
        BAPoint::Ptr point3d_;
        Vector2d point2d_;
        int id_;
    public:
        void setId(int id);
        void setPose(BAPose::Ptr pose);
        void setPoint(BAPoint::Ptr point3d);
        void setMeasurement(Vector2d point2d);
    };

    enum BA_TYPE {BA_FULL, BA_POSE, BA_POINT};


    class BAOptimizer {
    public:
        vector<BAPose::Ptr> poses_;
        vector<BAPoint::Ptr> points_;
        vector<BAError::Ptr> errors_;


    public:
        BAOptimizer();

        ~BAOptimizer();

        void setCamera(Camera::Ptr camera);

        void setType(BA_TYPE);

        void setVerbose(bool verbose);

        void addPose(BAPose::Ptr);

        void addPoint(BAPoint::Ptr);

        void addError(BAError::Ptr);

        void optimize(int niters);


    protected:
        void initial();

        Matrix<double, 2, 9> Jacobian(BAError::Ptr term);

        double error2(BAError::Ptr term);

        void computeHg();

        void constructEquation();

        void marginalize();

        void solveEquation();

        void updateStates();

        void recoverStates();

        Matrix3d antisymmetic(const Vector3d a) {
            Matrix3d a_as;
            a_as << 0, -a[2], a[1],
                    a[2], 0, -a[0],
                    -a[1], a[0], 0;
            return a_as;
        }

        // variables for normal equation
        vector<Matrix<double, 6, 6>> Hcc_;   // each element stores a 6x6 matrix for camera-camera, top-left part of Hessian matrix
        vector<Matrix3d> Hpp_;    // each element stores a 3x3 matrix for point-point, bottum-right part of Hessian matrix
        vector<vector<Matrix<double, 6, 3>>> Hcp_;// each element stores a 6x3 matrix for camera-point, top-right part of Hessian matrix
        VectorXd g_;//  -Jc*e for camera   -Jp*e for point
        vector<Matrix3d> HppInv_;    // store inverse of Hpp
        MatrixXd S_;                            //  top-left of Hessian after Shur elimination
        VectorXd gs_;       //  camera part of g after Shur elimination
        VectorXd delta_x_;  // pose and point increasement
        double maxDiaH_;                    // max value of the diagonal of Hessian

        // camera intrisics
        Camera::Ptr camera_;

        // convergence conditions
        const double epsilon1_ = 1e-6;
        const double epsilon2_ = 1e-6;
        const double tau_ = 1e-8;
        double upsilon_ = 2.0;
        double mu_;
        double sumError2_;


        BA_TYPE ba_type_;

        //MatrixXd HH_;
       // VectorXd gg_;
       // VectorXd DXXX_;
       // VectorXd e_;

        //double maxDiaHH_;
        bool verbose_ = true;

    };
}

#endif //VSLAM_BUNDLEADJUSTMENT_H
