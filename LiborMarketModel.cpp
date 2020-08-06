#include"LiborMarketModel.h"

template <class T>
LMM<T>::LMM() {
	std::cout << "Running constructor..." << std::endl;

	//sigCapFloor_ = sigCapFloor;
	//sigSwaption_ = sigSwaption;
	//corr_ = corr;
	//term_ = term;

}

template <class T>
LMM<T>::~LMM() {
	std::cout << "Running destructor..." << std::endl;
}

template <class T>
void LMM<T>::calibrate(const Eigen::ArrayXd& sigCapFloor, const Eigen::ArrayXd& term) {

	term_ = term;
	sigCapFloor_ = sigCapFloor;

	std::cout << "Testing the LMM calibration function..." << std::endl;
	Eigen::VectorXd zInit(4); zInit << 0, 0.2, 2, 0.2;
	std::cout << "zInit: " << zInit.transpose() << std::endl;
	//Eigen::VectorXd zSoln(2); zSoln << 1.0, 3.0;
	//std::cout << "zSoln: " << zSoln.transpose() << std::endl;

	LMMFunctor functor(sigCapFloor, term);
	Eigen::NumericalDiff<LMMFunctor> numDiff(functor);
	Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LMMFunctor>, double> lm(numDiff);
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
	std::cout << "Market Sigmas: " << sigCapFloor.transpose() << std::endl;
	Eigen::ArrayXd sigOut(sigCapFloor.size());
	for (int i = 0; i < sigCapFloor.size(); i++) {
		sigOut(i) = (z(0) + z(1) * term(i)) * exp(-z(2) * term(i)) + z(3);
	}
	std::cout << "Calibrated Sigmas: " << sigOut.transpose() << std::endl;
	z_ = z;
}

template<class T>
void LMM<T>::MC(const float& horizon, const int& gridCnt, const int& scens) {
	std::cout << "Running LMM Monte Carlo ..." << std::endl;

	float dt = horizon / gridCnt;
	float tempSig1, tempSig2, tempCorr, beta1, beta2;

	Eigen::MatrixXd cov(gridCnt, gridCnt);
	//Create covariance matrix
	//For correlation, will borrow function from https://docs.fincad.com/support/developerfunc/mathref/LIBORMarketModel.htm
	//rho(ij) = beta1 + (1-beta1) * exp(-beta2 * tau)
	//assume for now, beta1 = 0; beta2 = 0.1
	beta1 = 0, beta2 = 0.1;

	for (int ia = 0; ia < gridCnt;ia++) {
		for (int ib = 0;ib < gridCnt; ib++) {
			tempSig1 = (z_(0) + z_(1) * dt * (ia + 1)) * exp(-z_(2) * dt * (ia + 1)) + z_(3);
			tempSig2 = (z_(0) + z_(1) * dt * (ib + 1)) * exp(-z_(2) * dt * (ib + 1)) + z_(3);
			tempCorr = beta1 + (1 - beta1) * exp(-beta2 * abs(dt * (ia + 1) - dt * (ib + 1)));
			cov(ia, ib) = tempCorr * tempSig1 * tempSig2;
		}
	}
	std::cout << "cov: " << std::endl << cov << std::endl;

	Eigen::LLT<Eigen::MatrixXd> chol(cov);// = cov.llt();
	//Eigen::VectorXd c = chol.vectorD();
	//std::cout << "c: " << std::endl << c << std::endl;
	//int p = cov.cols();
	//Eigen::MatrixXd b = Eigen::MatrixXd::Identity(p, p);
	//b = chol.solve(b);
	//std::cout << "b: " << std::endl << b << std::endl;
	Eigen::MatrixXd L = chol.matrixL();
	std::cout << "L: " << std::endl << L << std::endl;


	Eigen::VectorXd zVec(gridCnt);
	boost::mt19937 rng(time(0));
	boost::normal_distribution<> nd(0.0, 1.0);
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > randNormal(rng, nd);
	double d = randNormal();
	for (int ia = 0;ia < scens;ia++) {
		for (int ib = 0; ib < gridCnt; ib++) {
			zVec(ib) = randNormal();
		}
		std::cout << "z vector: " << std::endl << zVec << std::endl;
		std::cout << "corr random: " << std::endl << L * zVec << std::endl;
	}
	
};

