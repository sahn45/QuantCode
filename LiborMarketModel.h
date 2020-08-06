/*
Notes & References:
1. Libor Market Model
http://lesniewski.us/papers/presentations/RiskCourse_LMM.pdf
https://docs.fincad.com/support/developerfunc/mathref/LIBORMarketModel.htm

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
#include<time.h>


template <class T>
class LMM {
private:
	Eigen::ArrayXd sigCapFloor_;
	//Eigen::MatrixXd sigSwaption_;
	//Eigen::MatrixXd corr_;
	Eigen::ArrayXd term_;
	Eigen::VectorXd z_;

public:
	LMM();
	~LMM();
	void calibrate(const Eigen::ArrayXd &sigCapFloor, const Eigen::ArrayXd &term);
	void MC(const float &horizon, const int &gridCnt, const int &scens);

};
