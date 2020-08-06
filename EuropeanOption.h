
/*
Notes & References:
1. Libor Market Model
http://lesniewski.us/papers/presentations/RiskCourse_LMM.pdf
https://docs.fincad.com/support/developerfunc/mathref/LIBORMarketModel.htm

2. For Levenberg Marquardt impelementation in Eigen library
https://medium.com/@sarvagya.vaish/levenberg-marquardt-optimization-part-2-5a71f7db27a0
https://stackoverflow.com/questions/18509228/how-to-use-the-eigen-unsupported-levenberg-marquardt-implementation

3. SABR Model
https://fincad.com/resources/resource-library/article/manage-smile-risk-sabr-model-stochastic-volatility

*/

#pragma once

#include<iostream>
#include<Eigen/Dense>
#include<vector>
#include<boost/math/distributions/normal.hpp>
#include"Tools.cpp"
#include<vector>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random.hpp>


template <class T>
class EuropeanOption {
private:
	T S_, K_, r_, sig_, T_;

public:
	EuropeanOption(const T &S, const T &K, const T &r, const T &sig, const T &T) ;
	~EuropeanOption();
	void Pricer_BSE(const char *callput);
	void Pricer_BM(const char *callput, const int &gridcount);
	void Pricer_FD(const char *callput, const char *fdType, const int &timegridcnt, const int &stockgridcnt, const int &upperMultiple);

};
