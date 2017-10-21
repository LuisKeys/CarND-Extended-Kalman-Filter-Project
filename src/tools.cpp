#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	VectorXd sum(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() == 0 ||
		estimations.size() != ground_truth.size()) {
		cout << "Incorrect estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		VectorXd estimation = estimations[i];
		VectorXd ground_truth_item = ground_truth[i];

		VectorXd residual = estimation - ground_truth_item;

		VectorXd square = residual.array() * residual.array();
		sum += square;
	}

	//calculate mean
	VectorXd mean = sum / estimations.size();

	//calculate the squared root
	rmse = mean.array().sqrt();

	//return the result
	return rmse;	
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//Precalc powers
	float px_2 = px * px;
	float py_2 = py * py;
	float diag_2 = px_2 + py_2;
	float g = sqrt(diag_2);
	float g_3_2 = g * g * g;

	//check division by zero
	if (diag_2 < 0.0001)
		cout << "Error, division by 0";

	//compute the Jacobian matrix
	//Row 1	
	Hj(0, 0) = px / g;
	Hj(0, 1) = py / g;
	Hj(0, 2) = 0;
	Hj(0, 3) = 0;

	//Row 2
	Hj(1, 0) = -py / diag_2;
	Hj(1, 1) = px / diag_2;
	Hj(1, 2) = 0;
	Hj(1, 3) = 0;

	//Row 3	
	Hj(2, 0) = py * (vx * py - vy * px) / g_3_2;
	Hj(2, 1) = py * (vy * px - vx * py) / g_3_2;
	Hj(2, 2) = px / g;
	Hj(2, 3) = py / g;


	return Hj;	
}
