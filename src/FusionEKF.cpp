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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // measurement matrix Laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // measurement matrix Radar
  Hj_ = MatrixXd(3, 4);
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

  // state vector
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 1, 1;

  // state covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  // state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // process covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);;

  //Delta time
  ekf_.dt_ = 0;
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
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float px = ro * cos(theta);
      float py = ro * sin(theta);

      ekf_.x_ << px, py, 
                0, 0;
      ekf_.timestamp_ = measurement_pack.raw_measurements_[3];
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 
                0, 0;
      ekf_.timestamp_ = measurement_pack.raw_measurements_[2];
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.timestamp_ = measurement_pack.raw_measurements_[3];
  }
  else {
    ekf_.timestamp_ = measurement_pack.raw_measurements_[2];
  }

  ekf_.dt_ = ((float)(ekf_.timestamp_ - previous_timestamp_)) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = ekf_.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float noise_ax = 9;
  float noise_ay = 9;

  ekf_.F_(0, 2) = ekf_.dt_;
  ekf_.F_(1, 3) = ekf_.dt_;
  
  //2
  float dt_2 = ekf_.dt_ * ekf_.dt_;
  float dt_3 = dt_2 * ekf_.dt_;
  float dt_4 = dt_3 * ekf_.dt_;

  ekf_.Q_ = MatrixXd(4, 4);

  //Q row 1
  ekf_.Q_(0, 0) = dt_4 / 4 * noise_ax;
  ekf_.Q_(0, 1) = 0;
  ekf_.Q_(0, 2) = dt_3 / 2 * noise_ax;
  ekf_.Q_(0, 3) = 0;

  //Q row 2
  ekf_.Q_(1, 0) = 0;
  ekf_.Q_(1, 1) = dt_4 / 4 * noise_ay;
  ekf_.Q_(1, 2) = 0;
  ekf_.Q_(1, 3) = dt_3 / 2 * noise_ay;

  //Q row 3
  ekf_.Q_(2, 0) = dt_3 / 2 * noise_ax;
  ekf_.Q_(2, 1) = 0;
  ekf_.Q_(2, 2) = dt_2 * noise_ax;
  ekf_.Q_(2, 3) = 0;

  //Q row 4
  ekf_.Q_(3, 0) = 0;
  ekf_.Q_(3, 1) = dt_3 / 2 * noise_ay;
  ekf_.Q_(3, 2) = 0;
  ekf_.Q_(3, 3) = dt_2* noise_ay;  

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
    // Radar updates
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0], 
         measurement_pack.raw_measurements_[1];

    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "dt_ = " << ekf_.dt_ << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout <<  endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
