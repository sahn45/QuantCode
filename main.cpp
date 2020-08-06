#include"EuropeanOption.h"
#include"EuropeanOption.cpp"
#include"LiborMarketModel.h"
#include"LiborMarketModel.cpp"
#include"SABR.h"
#include"PCA.h"
//#include"SABR.cpp"
#include<iostream>
#include<vector>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<boost/math/distributions/normal.hpp>
#include<boost/random.hpp>
#include<boost/random/normal_distribution.hpp>
#include<Eigen/Dense>
#include<unsupported/Eigen/NonLinearOptimization>
#include<unsupported/Eigen/NumericalDiff>
#include<string>
#include<time.h>
#include<tuple>

using boost::math::normal;

int main() {

	/*
	Eigen::MatrixXf data = Eigen::MatrixXf::Random(6, 6);
	float target = 0.95f;
	PCA p;
	p.run(data, &target, "pct");
	*/
	/*
	Eigen::VectorXd forward(6);
	Eigen::VectorXd tenor(6);
	Eigen::VectorXd capfloorSig(6);
	double beta = 1;
	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> calibPack;

	forward << 0.032, 0.034, 0.036, 0.037, 0.039, 0.0401;
	tenor << 1,2,3,4,5,6;
	capfloorSig<< 0.22, 0.24, 0.23, 0.25, 0.26, 0.265 ;
	calibPack = std::make_tuple(forward, tenor, capfloorSig, beta);

	sabr s(calibPack);
	s.sabr_calibrate();
	*/
	/*
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	Eigen::ArrayXd capfloorSig(6);
	capfloorSig << 0.22, 0.24, 0.23, 0.25, 0.26, 0.265;
	LMM_Calibrate(capfloorSig);

	Eigen::MatrixXd sig(5, 1);
	Eigen::MatrixXd corr(5, 5);
	sig << 0.22, 0.23, 0.23, 0.25, 0.26;
	std::cout << "sig: " << std::endl << sig << std::endl;
	corr << 1, 0.99, 0.99, 0.99, 0.99,
		0.99, 1, 0.99, 0.99, 0.99,
		0.99, 0.99, 1, 0.99, 0.99,
		0.99, 0.99, 0.99, 1, 0.99,
		0.99, 0.99, 0.99, 0.99, 1;
	std::cout << "corr: " << std::endl << corr << std::endl;
	LMM_MC(sig, corr);
	*/

	//testBoothFun();
	//testHimmelblauFun();

	
	Eigen::ArrayXd sigCapFloor(5);
	sigCapFloor << 0.225, 0.216, 0.206, 0.202, 0.193;
	Eigen::ArrayXd term(5);
	term << 1, 2, 3, 4, 5;

	LMM<float> L;
	L.calibrate(sigCapFloor, term);

	Eigen::MatrixXd corr(5, 5);
	corr << 1, 0.98, 0.98, 0.98, 0.98,
		0.98, 1, 0.98, 0.98, 0.98,
		0.98, 0.98, 1, 0.98, 0.98,
		0.98, 0.98, 0.98, 1, 0.98,
		0.98, 0.98, 0.98, 0.98, 1;
	L.MC(1, 10, 100);
	

	
	/*EuropeanOption<float> t1(100, 95, 0.02, 0.2, 1);

	t1.Pricer_BSE("call");
	t1.Pricer_BSE("put");
	t1.Pricer_BM("call", 100);
	t1.Pricer_BM("put", 100);
	t1.Pricer_FD("call", "explicit", 50, 50, 3);
	t1.Pricer_FD("put", "explicit", 50, 50, 3);
	t1.Pricer_FD("call", "implicit", 50, 50, 3);
	t1.Pricer_FD("put", "implicit", 50, 50, 3);
	t1.Pricer_FD("call", "crank-nicolson", 50, 50, 3);
	t1.Pricer_FD("put", "crank-nicolson", 50, 50, 3);*/
	

	return 0;
};

