//
// Created by randy on 19-6-14.
//

#include "bundleadjustment.h"
#include <algorithm>
#include <iomanip>
#include "run_timer.h"
//#include <chrono>

namespace vslam {
    void BAPose::setId(int id) {
        id_ = id;
    }

    void BAPose::setEstimate(SE3 pose) {
        pose_ = pose;
    }

    void BAPoint::setId(int id) {
        id_ = id;
    }

    void BAPoint::setEstimate(Vector3d point) {
        point_ = point;
    }

    void BAError::setId(int id) {
        id_ = id;
    }

    void BAError::setPose(BAPose::Ptr pose) {
        pose_ = pose;
    }

    void BAError::setPoint(BAPoint::Ptr point3d) {
        point3d_ = point3d;
    }

    void BAError::setMeasurement(Vector2d point2d) {
        point2d_ = point2d;
    }

    BAOptimizer::BAOptimizer() : camera_(nullptr), ba_type_(BA_FULL), upsilon_(2.0), verbose_(true) {}

    BAOptimizer::~BAOptimizer(){
    }

    void BAOptimizer::setCamera(Camera::Ptr camera){
        camera_ = camera;
    }

    void BAOptimizer::setType(BA_TYPE ba_type) {
        if (ba_type != BA_FULL && ba_type != BA_POSE && ba_type != BA_POINT) {
            cout << "usage: set bundle adjustment type as: FULL, POSE or POINT..." << endl;
        }
        else {
            ba_type_ = ba_type;
        }
    }

    void BAOptimizer::setHuber(double huber) {
        huber_ = huber;
    }

    void BAOptimizer::setVerbose(bool verbose) {
        verbose_ = verbose;
    }

    void BAOptimizer::addPose(BAPose::Ptr pose) {
        poses_.push_back(pose);
    }

    void BAOptimizer::addPoint(BAPoint::Ptr point) {
        points_.push_back(point);
    }
    void BAOptimizer::addError(BAError::Ptr error) {
        errors_.push_back(error);
    }

    void BAOptimizer::optimize(int niters) {
        initial();

        int k = 0;

        computeHg();

        double lastSumError2 = sumError2_;


        bool found = g_.lpNorm<Eigen::Infinity>() <= epsilon1_;

        mu_ = tau_ * maxDiaH_;

        double total_time = 0;

        while (!found && k < niters) {

            Runtimer t;
            t.start();
            k++;

            constructEquation();

            marginalize();

            solveEquation();


            if (delta_x_.norm() <= epsilon2_) {
                found = true;
            }
            else {
                //double curSumError2 = 0;
                updateStates();
                computeHg();

                double rho = (lastSumError2 - sumError2_)/(delta_x_.transpose()*(mu_*delta_x_+g_));
                if (rho > 0) {      // step acceptable
                    lastSumError2 = sumError2_;
                    mu_ = mu_ * std::max<double>(1.0/3.0f, 1-std::pow((2*rho-1), 3));
                    upsilon_ = 2;
                    found = g_.lpNorm<Eigen::Infinity>() <= epsilon1_;
                }
                else {
                    recoverStates();
                    computeHg();
                    mu_ = mu_ * upsilon_;
                    upsilon_ = 2 * upsilon_;
                }
            }

            t.stop();
            total_time += t.duration();

            if(verbose_)
            {
                std::ios::fmtflags f(cout.flags());
                cout << std::fixed << "Iter: " << std::left <<std::setw(4) << k
                          << " Cost: "<< std::left <<std::setw(20)  << std::setprecision(10) << sumError2_
                          << " Step: " << std::left <<std::setw(14) << std::setprecision(10) << delta_x_.norm()
                          << " Time " << std::left <<std::setw(10) << std::setprecision(3) << t.duration()
                          << " Total_time " << std::left <<std::setw(10) << std::setprecision(3) << total_time << endl;
                cout.flags(f);
            }

        }


    }

    void BAOptimizer::initial() {
        // determine the size of each Hessian block
        if (ba_type_ == BA_FULL) {
            Hcc_.resize(poses_.size());
            Hpp_.resize(points_.size());
            Hcp_.resize(poses_.size());
            for (int i = 0; i < Hcp_.size(); i++) {
                Hcp_[i].resize(points_.size());
            }
            g_.resize(6 * poses_.size() + 3 * points_.size());
            S_.resize(6 * Hcc_.size(), 6 * Hcc_.size());
            gs_.resize(6 * Hcc_.size());
            HppInv_.resize(Hpp_.size());
            delta_x_.resize(6*poses_.size()+3*points_.size());
        }
        if (ba_type_ == BA_POSE) {
            Hcc_.resize(poses_.size());
            g_.resize(6*poses_.size());
            S_.resize(6 * Hcc_.size(), 6 * Hcc_.size());
            delta_x_.resize(6*poses_.size());
        }
        if (ba_type_ == BA_POINT) {
            Hpp_.resize(points_.size());
            g_.resize(3*points_.size());
            S_.resize(3 * Hpp_.size(), 3 * Hpp_.size());
            delta_x_.resize(3*points_.size());
        }
    };

    Matrix<double, 2, 9> BAOptimizer::Jacobian(BAError::Ptr term) {
        Matrix3d R = term->pose_->pose_.rotation_matrix();
        Vector3d t = term->pose_->pose_.translation();
        Vector3d pointCamera = R * term->point3d_->point_ + t;
        double x = pointCamera[0], y = pointCamera[1], z = pointCamera[2];
        float fx = camera_->fx_, fy = camera_->fy_;
        Matrix<double, 2, 3> delta_e_pc;
        delta_e_pc << fx/z, 0, -fx*x/(z*z),
                0, fy/z, -fy*y/(z*z);
        delta_e_pc = -delta_e_pc;

        Matrix<double, 3, 6> delta_pc_xi;
        delta_pc_xi << Matrix3d::Identity(), -antisymmetic(pointCamera);

        // Jacobian for camera pose
        Matrix<double, 2, 6> J_pose;
        J_pose = delta_e_pc * delta_pc_xi;


        // Jacobian for point
        Matrix<double, 2, 3> J_point;
        J_point = delta_e_pc * R;

        Matrix<double, 2, 9> J;
        J << J_pose, J_point;

        return J;
    }

    double BAOptimizer::error2(BAError::Ptr term) {
        Vector2d error = term->point2d_ - camera_->world2pixel(term->point3d_->point_, term->pose_->pose_);
        return error.transpose()*error;
    }

    void BAOptimizer::computeHg(){

        sumError2_ = 0;
        for (int i = 0; i < Hcc_.size(); i++) {
            Hcc_[i].setZero();
        }
        for (int i = 0; i < Hpp_.size(); i++) {
            Hpp_[i].setZero();
        }
        for (int i = 0; i < Hcp_.size(); i++) {
            for (int j = 0; j < Hcp_[i].size(); j++) {
                Hcp_[i][j].setZero();
            }
        }
        g_.setZero();
        maxDiaH_ = 0;


        // compute J, H, error, J*error, using block matrix operations
        for (int i = 0; i < errors_.size(); i++) {
            BAError::Ptr term = errors_[i];

            Matrix<double, 2, 9> J = Jacobian(term);
            Matrix<double, 2, 6> Jc = J.block(0, 0, 2, 6);  // Jacobian for camera pose, using sub-matrix
            Matrix<double, 2, 3> Jp = J.block(0, 6, 2, 3);  // Jacobian for points in world reference
            Vector2d e;     // error = measurement - estimated
            e = term->point2d_ - camera_->world2pixel(term->point3d_->point_, term->pose_->pose_);

            if (ba_type_ == BA_FULL) {
                Hcc_[term->pose_->id_] += Jc.transpose() * Jc;
                Hpp_[term->point3d_->id_] += Jp.transpose() * Jp;
                Hcp_[term->pose_->id_][term->point3d_->id_] += Jc.transpose() * Jp;
                g_.segment(6*(term->pose_->id_), 6) -= Jc.transpose() * e;
                g_.segment(6*poses_.size() + 3*(term->point3d_->id_), 3) -= Jp.transpose() * e;
            }

            if (ba_type_ == BA_POSE) {
                Hcc_[term->pose_->id_] += Jc.transpose() * Jc;
                g_.segment(6*(term->pose_->id_), 6) -= Jc.transpose() * e;
            }

            if (ba_type_ == BA_POINT) {
                Hpp_[term->point3d_->id_] += Jp.transpose() * Jp;
                g_.segment(3*(term->point3d_->id_), 3) -= Jp.transpose() * e;

            }
            if (e.norm() < huber_)
                sumError2_ += e.transpose()*e;
            else
                sumError2_ += 2*huber_*e.norm() - huber_*huber_;
        }

        if (ba_type_ == BA_FULL || ba_type_ == BA_POSE) {
            for (int i = 0; i < Hcc_.size(); i++) {
                double tmaxDiaH = Hcc_[i].diagonal().lpNorm<Eigen::Infinity>();
                maxDiaH_ = std::max(tmaxDiaH, maxDiaH_);
            }
        }
        if (ba_type_ == BA_FULL || ba_type_ == BA_POINT) {
            for (int i = 0; i < Hpp_.size(); i++) {
                double tmaxDiaH = Hpp_[i].diagonal().lpNorm<Eigen::Infinity>();
                maxDiaH_ = std::max(tmaxDiaH, maxDiaH_);
            }
        }
    }

    void BAOptimizer::constructEquation() {
        if (ba_type_ == BA_FULL || ba_type_ == BA_POSE) {
            for (int i = 0; i < Hcc_.size(); i++) {
                Hcc_[i] += mu_ * MatrixXd::Identity(Hcc_[i].rows(), Hcc_[i].cols());
            }
        }
        if (ba_type_ == BA_FULL || ba_type_ == BA_POINT) {
            for (int i = 0; i < Hpp_.size(); i++) {
                Hpp_[i] += mu_ * MatrixXd::Identity(Hpp_[i].rows(), Hpp_[i].cols());
            }
        }

    }

    void BAOptimizer::marginalize() {

        S_.setZero();
        gs_.setZero();

        if (ba_type_ == BA_FULL) {
            // assign Hcc to diagonal block of S
            for (int k = 0; k < Hcc_.size(); k++) {
                S_.block(6 * k, 6 * k, 6, 6) = Hcc_[k];
            }
            // assign pose part of g to gs
            gs_ = g_.segment(0, 6*poses_.size());
            // compute inverse of Hpp
            for (int k = 0; k < Hpp_.size(); k++)
                HppInv_[k] = Hpp_[k].inverse();
            // compute S = Hcc - Hcp*HppInv*HcpT
            // compute gcs = gc - Hcp*HppInv*gp
            // S is symmetric, so just compute the upper trianglar and diagonal block
            for (int k = 0; k < points_.size(); k++) {
                for (int i = 0; i < poses_.size(); i++) {
                    Matrix<double, 6, 3> HcpiHppInvk = Hcp_[i][k] * HppInv_[k];
                    for (int j = i; j < poses_.size(); j++) {
                        S_.block(6 * i, 6 * j, 6, 6) -= HcpiHppInvk * Hcp_[j][k].transpose();
                    }
                    gs_.segment(6 * i, 6) -= HcpiHppInvk * g_.segment(6 * poses_.size() + 3 * k, 3);
                }
            }
            // determine the lower trianglar block of S
            for (int i = 0; i < poses_.size(); i++) {
                for (int j = i+1; j < poses_.size(); j++) {
                    S_.block(6 * j, 6 * i, 6, 6) = S_.block(6 * i, 6 * j, 6, 6).transpose();
                }
            }

        }

        if (ba_type_ == BA_POSE) {
            for (int k = 0; k < Hcc_.size(); k++) {
                S_.block(6 * k, 6 * k, 6, 6) = Hcc_[k];
            }
        }

        if (ba_type_ == BA_POINT) {
            for (int k = 0; k < Hpp_.size(); k++) {
                S_.block(3*k, 3*k, 3, 3) = Hpp_[k];
            }
        }
    }

    void BAOptimizer::solveEquation(){
        delta_x_.setZero();
        // solve increacement for camera pose
        if (ba_type_ == BA_FULL) {
            delta_x_.segment(0, 6*poses_.size()) = S_.fullPivHouseholderQr().solve(gs_);

            // back substitution to solve increacement for point
            for (int k = 0; k < points_.size(); k++) {
                Vector3d Temp = g_.segment(6*poses_.size() + 3 * k, 3);
                for (int i = 0; i < poses_.size(); i++) {
                    Temp -= Hcp_[i][k].transpose() * delta_x_.segment(6 * i, 6);
                }
                delta_x_.segment(6*poses_.size() + 3*k, 3) = HppInv_[k] * Temp;
            }
        }
        if (ba_type_ == BA_POSE || ba_type_ == BA_POINT) {
            delta_x_ = S_.fullPivHouseholderQr().solve(g_);
        }

    }

    void BAOptimizer::updateStates(){
        if (ba_type_ == BA_FULL) {
            // for camera pose, using SE3 left multiply
            for (int i = 0; i < poses_.size(); i++) {
                BAPose::Ptr & pose_ = poses_[i];
                VectorXd delta_pose_ = delta_x_.segment(6 * pose_->id_, 6);
                Sophus::SE3 delta_T_ = SE3::exp(delta_pose_);
                pose_->pose_ = delta_T_ * pose_->pose_;
            }
            // for map point, using plus
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr & point_ = points_[i];
                Vector3d delta_point_ = delta_x_.segment(6*poses_.size() + 3 * point_->id_, 3);
                point_->point_ = delta_point_ + point_->point_;
            }
        }
        if (ba_type_ == BA_POSE) {
            for (int i = 0; i < poses_.size(); i++) {
                BAPose::Ptr &pose_ = poses_[i];
                VectorXd delta_pose_ = delta_x_.segment(6 * pose_->id_, 6);
                Sophus::SE3 delta_T_ = SE3::exp(delta_pose_);
                pose_->pose_ = delta_T_ * pose_->pose_;
            }
        }
        if (ba_type_ == BA_POINT) {
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr & point_ = points_[i];
                Vector3d delta_point_ = delta_x_.segment(3 * point_->id_, 3);
                point_->point_ = delta_point_ + point_->point_;
            }
        }
    }

    void BAOptimizer::recoverStates(){
        if (ba_type_ == BA_FULL) {
            // for camera pose, using SE3 left multiply
            for (int i = 0; i < poses_.size(); i++) {
                BAPose::Ptr & pose_ = poses_[i];
                VectorXd delta_pose_ = -delta_x_.segment(6 * pose_->id_, 6);
                Sophus::SE3 delta_T_ = SE3::exp(delta_pose_);
                pose_->pose_ = delta_T_ * pose_->pose_;
            }
            // for map point, using plus
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr & point_ = points_[i];
                Vector3d delta_point_ = -delta_x_.segment(6*poses_.size() + 3 * point_->id_, 3);
                point_->point_ = delta_point_ + point_->point_;
            }
        }
        if (ba_type_ == BA_POSE) {
            for (int i = 0; i < poses_.size(); i++) {
                BAPose::Ptr &pose_ = poses_[i];
                VectorXd delta_pose_ = -delta_x_.segment(6 * pose_->id_, 6);
                Sophus::SE3 delta_T_ = SE3::exp(delta_pose_);
                pose_->pose_ = delta_T_ * pose_->pose_;
            }
        }
        if (ba_type_ == BA_POINT) {
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr & point_ = points_[i];
                Vector3d delta_point_ = -delta_x_.segment(3 * point_->id_, 3);
                point_->point_ = delta_point_ + point_->point_;
            }
        }
    }

}


