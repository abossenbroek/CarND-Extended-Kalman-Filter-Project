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

    VectorXd z_conv = VectorXd(4);
    float rho = z(MeasurementPackage::RHO);
    float phi = z(MeasurementPackage::PHI);
    float rho_dot = z(MeasurementPackage::RHO_DOT);
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    */
    float c_phi = cos(phi);
    float s_phi = sin(phi);

    float px = rho * c_phi;
    float py = rho * s_phi;
    float vx = rho_dot * c_phi;
    float vy = rho_dot * s_phi;

    z_conv << px, py, vx, vy;

    Tools t;

    VectorXd y;
    MatrixXd S;
    MatrixXd K;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    MatrixXd Ht = H_.transpose();

    y = z - t.CalculateJacobian(x_) * x_;
    S = H_ * P_ * Ht + R_;
    K = P_ * Ht * S.inverse();
    x_ = x_ + K * y;
    P_ = (I - K * H_) * P_;
}
