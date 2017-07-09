#include "kalman_filter.h"
#include "tools.h"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd y;
    MatrixXd S;
    MatrixXd K;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    MatrixXd Ht = H_.transpose();


    y = z - H_ * x_;
    S = H_ * P_ * Ht + R_;
    K = P_ * Ht * S.inverse();
    x_ = x_ + K * y;
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    VectorXd y;
    MatrixXd S;
    MatrixXd K;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    MatrixXd Ht = H_.transpose();

    float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
    float phi = atan2(x_(1), x_(0));
    float rho_d = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
    VectorXd h_prime = VectorXd(3);
    h_prime << rho, phi, rho_d;

    y = z - h_prime;
    // angle normalization
    while (y(1) >  M_PI) y(1) -= 2. * M_PI;
    while (y(1) < -M_PI) y(1) += 2. * M_PI;

    S = H_ * P_ * Ht + R_;
    K = P_ * Ht * S.inverse();
    x_ = x_ + K * y;
    P_ = (I - K * H_) * P_;
}
