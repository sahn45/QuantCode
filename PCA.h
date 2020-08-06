#pragma once
#include<iostream>
#include<tuple>
#include<Eigen/Dense>

using namespace std;

class PCA {
private:

public:
	PCA();
	~PCA();
	Eigen::EigenSolver<Eigen::MatrixXf> run(const Eigen::MatrixXf &data, const float *target, const char *flag);

};