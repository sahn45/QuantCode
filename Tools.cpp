#pragma once
#include<Eigen/Dense>
#include<boost/math/distributions/normal.hpp>
#include<boost/random.hpp>
#include<boost/random/normal_distribution.hpp>
#include<vector>

// Generic functor
// See http://eigen.tuxfamily.org/index.php?title=Functors
// C++ version of a function pointer that stores meta-data about the function
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{

	// Information that tells the caller the numeric type (eg. double) and size (input / output dim)
	typedef _Scalar Scalar;
	enum { // Required by numerical differentiation module
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};

	// Tell the caller the matrix sizes associated with the input, output, and jacobian
	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

	// Local copy of the number of inputs
	int m_inputs, m_values;

	// Two constructors:
	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	// Get methods for users to determine function input and output dimensions
	int inputs() const { return m_inputs; }
	int values() const { return m_values; }

};


// https://en.wikipedia.org/wiki/Test_functions_for_optimization
// LMM Function
// Implement f_sig(k(i),a,b,c,d) = ki[(a+b(ti-t)exp(-c(ti-t)))+d]
struct LMMFunctor : Functor<double>
{
	Eigen::ArrayXd sigCapFloor_, term_;

	// Simple constructor
	LMMFunctor(const Eigen::ArrayXd &sigCapFloor, const Eigen::ArrayXd &term) : Functor<double>(4, 5) {
		sigCapFloor_ = sigCapFloor;
		term_ = term;

	} //First parameter = number variables; second parameter = number of observations (i.e., count of fvec formulas used to calibrate)

	// Implementation of the objective function
	int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {
		//double x1 = z(0);   double x2 = z(1); double x3 = z(2); double x4 = z(3); double x5 = z(4);
		/*
		 * Evaluate the Booth function.
		 * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
		 * of squared terms. The algorithm takes this into account: do not do it yourself.
		 * In other words: objFun = sum(fvec(i)^2)
		 */

		for (int i = 0;i < sigCapFloor_.size();i++) {
			fvec(i) = ((z(0) + z(1) * term_(i))*exp(-z(2) * term_(i)) + z(3)) - sigCapFloor_(i);
		}

		return 0;
	}
};

struct SABRFunctor : Functor<double>
{
	Eigen::VectorXd forward_, tenor_, sigCapFloor_;
	double beta_;

	// Simple constructor
	SABRFunctor(const std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, double> &calibPack) : Functor<double>(3, 6) {
		forward_ = std::get<0>(calibPack);
		tenor_ = std::get<1>(calibPack);
		sigCapFloor_ = std::get<2>(calibPack);
		beta_ = std::get<3>(calibPack);
	} //First parameter = number variables; second parameter = number of observations (i.e., count of fvec formulas used to calibrate)

	// Implementation of the objective function
	int operator()(const Eigen::VectorXd &z, Eigen::VectorXd &fvec) const {

		for (int i = 0;i < forward_.size();i++) {
			fvec(i) = z(0) / pow(forward_(i), 1 - beta_) * (1.0 + (((1 - beta_)*(1 - beta_)*z(0)*z(0)) /
				(24 * pow(forward_(i), 2 - 2 * beta_)) + z(1)*beta_*z(0)*z(2) / (4 * pow(forward_(i), 1 - beta_)) + (2 - 3 * z(1)*z(1))*z(2)*z(2) / 24) * tenor_(i));
		}

		return 0;
	}
};

/*
template boost::normal_distribution<> nd(0.0, 1.0)
{
	return template boost::normal_distribution<>();
};
*/
/*
const double INF = 1.e100;
double interpolate(double x, std::vector<std::pair<double, double> > table) {
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > table.back().first) return INF;
	if (x < table[0].first) return -INF;
	std::vector<std::pair<double, double> >::iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
	// Corner case
	if (it == table.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first) / (it->first - it2->first);
};
*/