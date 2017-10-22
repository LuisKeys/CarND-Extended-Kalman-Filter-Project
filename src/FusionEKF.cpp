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

  ekf_.Init( ekf_.x_,
          ekf_.P_,
          ekf_.F_,
          H_laser_,
          R_laser_,
          ekf_.Q_ );
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
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 
                0, 0;
    }

    long x_size = ekf_.x_.size();
    ekf_.I_ = MatrixXd::Identity(x_size, x_size);

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  ekf_.dt_ = ((float)(measurement_pack.timestamp_ - previous_timestamp_)) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

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

  /*
  cout << measurement_pack.raw_measurements_[0] << " " <<
          measurement_pack.raw_measurements_[1] << endl;
  cout << endl;

  cout << "x_ = " << ekf_.x_ << endl;
  cout << endl;

  cout << "measurement_pack.timestamp_ = " << measurement_pack.timestamp_ << endl;
  cout << endl;
  */

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
    VectorXd z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],
         measurement_pack.raw_measurements_[1],
         measurement_pack.raw_measurements_[2];
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
  cout << "Test 1 L";
   ekf_.Update(measurement_pack.raw_measurements_);
  cout << "Test 2 L";
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
