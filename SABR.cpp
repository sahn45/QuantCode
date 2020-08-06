#include"SABR.h"

sabr::sabr(const std::tuple < Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double > &calibPack) {
	std::cout << "Running constructor..." << std::endl;
	calibPack_ = calibPack;
}

sabr::~sabr() {
	std::cout << "Running destructor..." << std::endl;
}

void sabr::sabr_calibrate() {//tuple elements[forwardRates, tenors, sigmas, beta_assumed]

	std::cout << "Testing the SABR calibrator..." << std::endl;
	Eigen::VectorXd zInit(4); zInit << 0.1, 0.1, 0.1, 0.1;
	std::cout << "zInit: " << zInit.transpose() << std::endl;
	//Eigen::VectorXd zSoln(2); zSoln << 1.0, 3.0;
	//std::cout << "zSoln: " << zSoln.transpose() << std::endl;

	SABRFunctor functor(calibPack_);
	Eigen::NumericalDiff<SABRFunctor> numDiff(functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<SABRFunctor>, double> lm(numDiff);
	lm.parameters.maxfev = 1000;
	lm.parameters.xtol = 1.0e-10;
	std::cout << "max fun eval: " << lm.parameters.maxfev << std::endl;
	std::cout << "x tol: " << lm.parameters.xtol << std::endl;

	Eigen::VectorXd z = zInit;
	int ret = lm.minimize(z);
	std::cout << "iter count: " << lm.iter << std::endl;
	std::cout << "return status: " << ret << std::endl;
	std::cout << "zSolver: " << z.transpose() << std::endl;
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

};