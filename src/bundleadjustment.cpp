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


     /*   if (verbose_) {
            cout << "*************************before optimization******************************" << endl;
            for (int i = 0; i < poses_.size(); i++) {
                BAPose::Ptr pose = poses_[i];
                cout << "pose " << pose->id_ << " is " << endl << pose->pose_.matrix() << endl;
            }
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr point = points_[i];
                cout << "point " << point->id_ << " is " << endl << point->point_ << endl;
            }

            for (int i = 0; i < errors_.size(); i++) {
                BAError::Ptr term = errors_[i];
                cout << "in error " << term->id_ << " pose " << term->pose_->id_ << " is " << endl
                     << term->pose_->pose_.matrix() << endl;
                cout << "point 3d " << term->point3d_->id_ << " is " << endl << term->point3d_->point_ << endl;
                cout << "point 2d is " << endl << term->point2d_ << endl;
            }
        }*/



        //cout << "g is " << endl << g_ << endl;
        //cout << "g max is " << g_.lpNorm<Eigen::Infinity>() << endl;

        bool found = g_.lpNorm<Eigen::Infinity>() <= epsilon1_;

        mu_ = tau_ * maxDiaH_;

        double total_time = 0;

        while (!found && k < niters) {

            Runtimer t;
            t.start();

         //   cout << "*********************** iterator " << k << " ***************************** " << endl;
            k++;

        //    cout << "before update, g is " << endl << g_ << endl;
         /*   cout << "before update, g max is " << g_.lpNorm<Eigen::Infinity>() << endl;
            for (int i=0; i < poses_.size(); i++) {
                cout << "pose " << i << " is " << endl << poses_[i]->pose_.matrix() << endl;
            }*/
         //   cout << "mu = " << mu_ << endl << "upsilon = " << upsilon_ << endl;

            constructEquation();
       //     cout << "normal equation building done..." << endl;
       //     cout <<"full HH is " << endl << HH_ << endl;
       //     cout << "full ee is " << endl << e_ << endl;
       //     cout << "full gg is" << endl << gg_ << endl;



            marginalize();
       //     cout << "marginalize done..." << endl;
        //    cout << "S = " << endl << S_ << endl;
        //    cout << "gs = " << endl << gs_ << endl;

            solveEquation();
        //    cout << "solve equation done..." << endl;
        //    cout << "delta x is " << endl << delta_x_ << endl;
        //    cout << "||delta x||2 is " << delta_x_.norm() << endl;

            if (delta_x_.norm() <= epsilon2_) {
                found = true;
            }
            else {
                //double curSumError2 = 0;
                updateStates();
                computeHg();

              //  cout << "after update, g is " << endl << g_ << endl;
              //  cout << "after update, g max is " << g_.lpNorm<Eigen::Infinity>() << endl;
            /*    for (int i=0; i < poses_.size(); i++) {
                    cout << "pose " << i << " is " << endl << poses_[i]->pose_.matrix() << endl;
                }*/
              /*  for (int i = 0; i < errors_.size(); i++) {
                    BAError::Ptr term = errors_[i];
                    curSumError2 += error2(term);
                } */
             //   cout << "last sum error2 is " << lastSumError2 << endl;
          //      cout << "current sum error2 is " << sumError2_ << endl;
                double rho = (lastSumError2 - sumError2_)/(delta_x_.transpose()*(mu_*delta_x_+g_));
             //   cout << "rho = " << rho << endl;
                if (rho > 0) {      // step acceptable
                    cout << "states update done..." << endl;
                    //cout << "temp = " << 1-std::pow((2*rho-1), 3) << endl;
                    lastSumError2 = sumError2_;
                    mu_ = mu_ * std::max<double>(1.0/3.0f, 1-std::pow((2*rho-1), 3));
                    upsilon_ = 2;
                    found = g_.lpNorm<Eigen::Infinity>() <= epsilon1_;
            //        cout << "mu = " << mu_ << endl << "upsilon = " << upsilon_ << endl;
                }
                else {
                    recoverStates();
                    computeHg();
            //        cout << "after recover, g is " << endl << g_ << endl;
             //       cout << "after recover, g max is " << g_.lpNorm<Eigen::Infinity>() << endl;

                  /*  for (int i=0; i < poses_.size(); i++) {
                        cout << "pose " << i << " is " << endl << poses_[i]->pose_.matrix() << endl;
                    }
                    for (int i=0; i < 5; i++) {
                        cout << "in error " << i << " pose " << errors_[i]->pose_->id_ << " is " << endl << errors_[i]->pose_->pose_.matrix() << endl;
                    }*/
                    mu_ = mu_ * upsilon_;
                    upsilon_ = 2 * upsilon_;
             //       cout << "mu = " << mu_ << endl << "upsilon = " << upsilon_ << endl;
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

       // cout << "************************* after optimization******************************" << endl;
       /* for (int i = 0; i < poses_.size(); i++){
            BAPose::Ptr pose = poses_[i];
            cout << "pose " << pose->id_ << " is " << endl << pose->pose_.matrix() << endl;
        } */
       /* for (int i = 0; i < points_.size(); i++){
            BAPoint::Ptr point= points_[i];
            cout << "point " << point->id_ << " is " << endl << point->point_ << endl;
        }*/
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



        //MatrixXd J_;
        //J_.resize(2*errors_.size(), 6*poses_.size()+3*points_.size());
        //J_.setZero();


        //e_.resize(2*errors_.size());
       // e_.setZero();

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



          //  J_.block(term->point3d_->id_*2,  term->pose_->id_ * 6, 2, 6) = Jc;

           // J_.block(term->point3d_->id_*2, 6*poses_.size()+ term->point3d_->id_*3, 2, 3) = Jp;

          //  e_.segment(term->point3d_->id_*2, 2) = e;


            if (ba_type_ == BA_POSE) {
                Hcc_[term->pose_->id_] += Jc.transpose() * Jc;
                g_.segment(6*(term->pose_->id_), 6) -= Jc.transpose() * e;
            }

            if (ba_type_ == BA_POINT) {
                Hpp_[term->point3d_->id_] += Jp.transpose() * Jp;
                g_.segment(3*(term->point3d_->id_), 3) -= Jp.transpose() * e;

            }
            sumError2_ += e.transpose()*e;
        }




      //  cout << "full J is" << endl << J_ << endl;
//
    //    HH_ = J_.transpose() * J_;


   //     gg_ = -J_.transpose() * e_;


   //     maxDiaHH_ = 0;

    //    maxDiaHH_ = HH_.diagonal().lpNorm<Eigen::Infinity>();
















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



  //      HH_ += tau_ * maxDiaHH_ * MatrixXd::Identity(HH_.rows(), HH_.cols());

//        MatrixXd BB_ = HH_.block(0, 0, 6*poses_.size(), 6*poses_.size());
//
 //       MatrixXd CC_ = HH_.block(6*poses_.size(), 6*poses_.size(), 3*points_.size(), 3*points_.size());

  //      MatrixXd EE_ = HH_.block(0, 6*poses_.size(), 6*poses_.size(), 3*points_.size());

   //     MatrixXd HHS_ = BB_ - EE_ * CC_.inverse() * EE_.transpose();

    //    cout << "HHS is" << endl << HHS_ << endl;







    }

    void BAOptimizer::marginalize() {

        S_.setZero();
        gs_.setZero();

//        typedef Matrix<double, 6, 1> Vector6d;
 //       vector<Vector6d> ECbp(poses_.size(), Vector6d::Zero());



 /************** for test **************************/
 /*       MatrixXd H_;
        H_.resize(6*poses_.size()+3*points_.size(), 6*poses_.size()+3*points_.size());
        H_.setZero();
        for (int i = 0; i < poses_.size(); i++) {
            H_.block(6*i, 6*i, 6, 6) = Hcc_[i];
        }
        for (int i = 0; i < points_.size(); i++) {
            H_.block(6*poses_.size()+3*i, 6*poses_.size()+3*i, 3, 3) = Hpp_[i];
        }
        for (int i = 0; i < poses_.size(); i++) {
            for (int j = 0; j < points_.size(); j++) {
                H_.block(6*i, 6*poses_.size()+3*j, 6, 3) = Hcp_[i][j];
            }
        }
        for (int j = 0; j < points_.size(); j++) {
            for (int i = 0; i < poses_.size(); i++) {
                H_.block(6*poses_.size()+3*j, 6*i,3,6) = Hcp_[i][j].transpose();
            }
        }


        cout << "H is" << endl << H_ << endl;

        MatrixXd B_ = H_.block(0, 0, 6*poses_.size(), 6*poses_.size());

        cout << "B is " << endl << B_ << endl;

        MatrixXd C_ = H_.block(6*poses_.size(), 6*poses_.size(), 3*points_.size(), 3*points_.size());
        cout << "C is " << endl << C_ << endl;


        MatrixXd E_ = H_.block(0, 6*poses_.size(), 6*poses_.size(), 3*points_.size());
        cout << "E is " << endl << E_ << endl;


        MatrixXd HS_ = B_ - E_ * C_.inverse() * E_.transpose();

        cout << "HS is" << endl << HS_ << endl;



        VectorXd GS_ = g_.segment(0, 6*poses_.size()) - E_ * C_.inverse() * g_.segment(6*poses_.size(), 3*points_.size());

        cout << "GS_ is" << endl << GS_ << endl;

        VectorXd Xc_ = HS_.fullPivLu().solve(GS_);


        VectorXd Xp_ = C_.inverse() * (g_.segment(6*poses_.size(), 3*points_.size()) - E_.transpose() * Xc_);

        VectorXd DX_;
        DX_.resize(g_.size());

        DX_ << Xc_, Xp_;

        cout << "Delta X is" << endl << DX_ << endl;


        DXXX_ = H_.fullPivLu().solve(g_);

        cout << "Full Delta X is" << endl << DXXX_ << endl;
*/


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

                  //  ECbp[i] += HcpiHppInvk * g_.segment(6 * poses_.size() + 3*k, 3);


                    for (int j = i; j < poses_.size(); j++) {
                        S_.block(6 * i, 6 * j, 6, 6) -= HcpiHppInvk * Hcp_[j][k].transpose();
                    }
                    gs_.segment(6 * i, 6) -= HcpiHppInvk * g_.segment(6 * poses_.size() + 3 * k, 3);
                }
            }
/*
            for (int k = 0; k < poses_.size(); k++) {
                gs_.segment(6 * k, 6) = g_.segment(6*k, 6) - ECbp[k];
            } */











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
                    //Temp += Hcp_[i][k].transpose() * delta_x_.segment(6 * i, 6);
                }
                delta_x_.segment(6*poses_.size() + 3*k, 3) = HppInv_[k] * Temp;
                //delta_x_.segment(6*poses_.size() + 3*k, 3) = HppInv_[k] * (g_.segment(6*poses_.size()+3*k,3) - Temp);

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
                //VectorXd delta_pose_ = DXXX_.segment(6*pose_->id_, 6);//delta_x_.segment(6 * pose_->id_, 6);
                Sophus::SE3 delta_T_ = SE3::exp(delta_pose_);
              //  cout << "delta pose is " << endl << delta_pose_ << endl;
              //  cout << "delta T is " << endl << delta_T_.matrix() << endl;
              //  cout << "pose is " << endl << pose_->pose_.matrix() << endl;

                pose_->pose_ = delta_T_ * pose_->pose_;
            }
            // for map point, using plus
            for (int i = 0; i < points_.size(); i++) {
                BAPoint::Ptr & point_ = points_[i];
                Vector3d delta_point_ = delta_x_.segment(6*poses_.size() + 3 * point_->id_, 3);
                //Vector3d delta_point_ = DXXX_.segment(6*poses_.size()+3*point_->id_,3);//delta_x_.segment(6*poses_.size() + 3 * point_->id_, 3);
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


