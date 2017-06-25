#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    // Dynamic allocate a vector of the size equal to estimations' first element.
    VectorXd rmse = VectorXd::Zero(estimations[0].size());

    if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
        cerr << "wrongly specified arguments" << endl;
        return rmse;
    }

    VectorXd res;

    // accumulate squared residuals
    for (int i = 0; i < estimations.size(); ++i) {
        res = (estimations[i] - ground_truth[i]);
        res = res.array() * res.array();
        rmse += res;
    }

    rmse = rmse / (float) estimations.size();

    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //TODO: YOUR CODE HERE
    float px_sq_py_sq = px * px + py * py;

    //check division by zero
    if (abs(px_sq_py_sq) < 1.0e-9) {
        cerr << "division by zero" << std::endl;
        return Hj;
    }
    float sq_px_py = sqrt(px_sq_py_sq);
    float thtw_px_py = pow(px_sq_py_sq, 3.0 / 2.0);

    //compute the Jacobian matrix
    Hj << px / sq_px_py, py / sq_px_py, 0, 0,
            -py / px_sq_py_sq, px / px_sq_py_sq, 0, 0,
            py * (vx * py - vy * px) / thtw_px_py, px * (vy * px - vx * py) / thtw_px_py,
            px / sq_px_py, py / sq_px_py;


    return Hj;
}
