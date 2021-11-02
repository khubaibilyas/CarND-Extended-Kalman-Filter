#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

#define EPS 0.0001

/**
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

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_  << 1.0, 1.0, 1.0, 1.0,
          1.0, 1.0, 1.0, 1.0,
  		  1.0, 1.0, 1.0, 1.0;

  MatrixXd P_in = Eigen::MatrixXd(4, 4);
  P_in    << 10.0,0.0,0.0,0.0,
             0.0,10.0,0.0,0.0,
             0.0,0.0,100.0,0.0,
             0.0,0.0,0.0,100.0;
    // Initial covariance Matrix
  ekf_.P_ = P_in;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  int noise_ax = 9;
  int noise_ay = 9;
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      float rho_measured = measurement_pack.raw_measurements_(0);
      float phi_measured = measurement_pack.raw_measurements_(1);
      float rhodot_measured = measurement_pack.raw_measurements_(2);
      
      if(phi_measured < -M_PI)
      {
        while(phi_measured < -M_PI)
        {
          phi_measured += 2 * M_PI;
        }
      }

      if(phi_measured > M_PI)
      {
        while(phi_measured > M_PI)
        {
           phi_measured -= 2 * M_PI;
        }
      }
      
      float px = rho_measured * cos(phi_measured);
      float py = rho_measured * sin(phi_measured);
      float vx = rhodot_measured * cos(phi_measured);
      float vy = rhodot_measured * sin(phi_measured);

      Eigen::VectorXd x_in;
      x_in = Eigen::VectorXd(4);
      x_in << px, py, vx, vy;
      
      Eigen::MatrixXd F_in;      
      F_in = Eigen::MatrixXd(4, 4);
  	  F_in << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0;
      
      Eigen::MatrixXd Q_in;
      Q_in = Eigen::MatrixXd(4, 4);
  	  Q_in << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0;

      ekf_.Init(x_in, ekf_.P_, F_in, Hj_, R_radar_, Q_in);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);
      float vx = 0;
      float vy = 0;

      Eigen::VectorXd x_in;
      x_in = Eigen::VectorXd(4);
      x_in << px, py, vx, vy;

      Eigen::MatrixXd F_in;
      F_in = Eigen::MatrixXd(4, 4);
  	  F_in << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0;
      
      Eigen::MatrixXd Q_in;
      Q_in = Eigen::MatrixXd(4, 4);
  	  Q_in << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0;

      ekf_.Init(x_in, ekf_.P_, F_in, H_laser_, R_laser_, Q_in);

    }
    
    // Deal with special case where values are < than a very small number
    if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
      ekf_.x_(0) = EPS;
      ekf_.x_(1) = EPS;
    }
    
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;

    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  

  float delta_T = measurement_pack.timestamp_ - previous_timestamp_;
  delta_T /= 1000000.0;
  
  Eigen::MatrixXd F_in;
  F_in = Eigen::MatrixXd(4, 4);
  F_in << 1.0, 0.0, delta_T, 0.0,
          0.0, 1.0, 0.0, delta_T,
          0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 1.0;

  float delta_T_pow2 = delta_T * delta_T;

  float delta_T_pow3 = delta_T_pow2 * delta_T;

  float delta_T_pow4 = delta_T_pow3 * delta_T;


  Eigen::MatrixXd Q_in;
  Q_in = Eigen::MatrixXd(4, 4);
  Q_in <<  (delta_T_pow4/4 * noise_ax), 0, (delta_T_pow3/2 * noise_ax), 0,
           0, (delta_T_pow4/4 * noise_ay), 0, (delta_T_pow3/2 * noise_ay),
           (delta_T_pow3/2 * noise_ax), 0, (delta_T_pow2 * noise_ax), 0,
           0, (delta_T_pow3/2 * noise_ay), 0, (delta_T_pow2 * noise_ay);

  ekf_.F_ = F_in;
  ekf_.Q_ = Q_in;
  
  cout << "Predicting..." << endl;

  ekf_.Predict();
  
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
      cout << "Updating Radar..." << endl;

      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      previous_timestamp_ = measurement_pack.timestamp_;

  } else {
    // TODO: Laser updates
      cout << "Updating Laser..." << endl;

      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
      previous_timestamp_ = measurement_pack.timestamp_;

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  
}
