#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);
  P_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  P_ << 200.0, 50.0, 30.0, 0.0, 
		50.0, 100.0, 0.0, 15.0,
        30.0, 0.0, 20.0, 5.0,
        0.0, 15.0, 5.0, 10.0;


  // Initialize transition matrix
  F_ << 1, 0, 0, 0,
	 0, 1, 0, 0,
       0, 0, 1, 0,
       0, 0, 0, 1;

  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;

  noise_ax = 9.0;
  noise_ay = 9.0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    VectorXd x_(4);
    x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    	double rho = measurement_pack.raw_measurements_[0]; //radius
    	double theta = measurement_pack.raw_measurements_[1]; //angle
        x_ << rho*cos(theta), rho*sin(theta),0.0,0.0;
        ekf_.Init(x_,P_,F_, Hj_,R_radar_,Q_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_<< measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],0.0,0.0;
      ekf_.Init(x_,P_,F_, H_laser_,R_laser_,Q_);
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // Initialize measurement matrix for laser measurements

    // Initialize ekf_ with the first state vector,
    // estimated initial state covariance matrix,
    // and an empty matrix for Q

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Prediction
   ***************************************************************************/

  // State transition matrix update
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4*noise_ax/4, 0, dt_3*noise_ax/2, 0,
  0, dt_4*noise_ay/4, 0, dt_3*noise_ay/2,
  dt_3*noise_ax/2, 0, dt_2*noise_ax, 0,
  0, dt_3*noise_ay/2, 0, dt_2*noise_ay;
  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  //return;
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_=R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
	  //	return;
	  ekf_.H_ = H_laser_;
	  ekf_.R_=R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
