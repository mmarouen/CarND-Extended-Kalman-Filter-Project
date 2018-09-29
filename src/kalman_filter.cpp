#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = P_ - K * H_*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	  float px = x_[0];
	  float py = x_[1];
	  float vx = x_[2];
	  float vy = x_[3];

	  // If rho == 0, skip the update step to avoid dividing by zero.
	  // This is crude but should be fairly robust on our data set.
	  if( px == 0. && py == 0. )
	    return;

	  VectorXd h_x(3);
	  float rho = sqrt( px*px + py*py );
	  h_x << rho, atan2( py, px ), ( px*vx + py*vy )/rho;

	  VectorXd y = z - h_x;
	  float pi=3.141592;
	  if( y[1] > pi )
	    y[1] =y[1]- 2.0*pi;
	  if( y[1] < -pi )
	    y[1] = y[1]+ 2.0*pi;
	  MatrixXd S = H_*P_*H_.transpose() + R_;
	  MatrixXd Sinv = S.inverse();
	  MatrixXd K =  P_*H_.transpose()*Sinv;

	  // Compute new state
	  x_ = x_ + K*y;
	  P_ = P_ - K*H_*P_;

  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
