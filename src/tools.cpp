#include "tools.h"
#include <iostream>

#define EPS 0.0001     // A very small number
#define EPS2 0.0000001 // A very very small number

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using std::cout;
using std::endl;
Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) 
  {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  if (fabs(px) < EPS and fabs(py) < EPS){
    px = EPS;
    py = EPS;
  }

  double px_py_sqd = (px*px + py*py);

  if (fabs(px_py_sqd) < EPS2) {
    px_py_sqd = EPS2;
  }
  
  double a11 = (px / sqrt(px_py_sqd));
  double a12 = (py / sqrt(px_py_sqd));
  double a13 = 0;
  double a14 = 0;
  
  double a21 = -(py / (px_py_sqd));
  double a22 =  (px / (px_py_sqd));
  double a23 = 0;
  double a24 = 0;
  
  double a31 = py * (vx*py - vy*px) / sqrt(px_py_sqd * px_py_sqd * px_py_sqd);
  double a32 = px * (vy*px - vx*py) / sqrt(px_py_sqd * px_py_sqd * px_py_sqd);
  double a33 = a11;
  double a34 = a12;

    // compute the Jacobian matrix
  Hj << a11, a12, a13, a14,
        a21, a22, a23, a24,
        a31, a32, a33, a34;

  return Hj;
}
