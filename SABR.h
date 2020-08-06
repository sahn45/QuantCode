/*
Notes & References:
1. SABR Model
https://fincad.com/resources/resource-library/article/manage-smile-risk-sabr-model-stochastic-volatility

2. For Levenberg Marquardt impelementation in Eigen library
https://medium.com/@sarvagya.vaish/levenberg-marquardt-optimization-part-2-5a71f7db27a0
https://stackoverflow.com/questions/18509228/how-to-use-the-eigen-unsupported-levenberg-marquardt-implementation
*/

#pragma once
#include"Tools.cpp"
#include<Eigen/Dense>
#include<unsupported/Eigen/NonLinearOptimization>
#include<unsupported/Eigen/NumericalDiff>
#include<iostream>
#include<tuple>

class sabr {
private:
	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> calibPack_;

public:
	sabr(const std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> &calibPack);
	~sabr();
	void sabr_calibrate();
};

