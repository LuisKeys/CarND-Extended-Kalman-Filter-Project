#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  /**
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  // If rho near 0, skip the update step to avoid dividing by zero.
  if( px < 0.0001 && py < 0.0001 )
    return;

  H_ = tools.CalculateJacobian( x_ );

  VectorXd transf_x = VectorXd(3);

  float ro = sqrt( px * px + py * py );
  transf_x << ro, atan2( py, px ), ( px * vx + py * vy ) / ro;

  // Update the state using EKF
  VectorXd y = z - transf_x;

  if( y[1] > PI )
    y[1] -= 2.0f * PI;

  if( y[1] < -PI )
    y[1] += 2.0f * PI;

  MatrixXd Hjt = H_.transpose();
  MatrixXd S = H_* P_* Hjt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Hjt * Si;

  // Compute new state
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}
