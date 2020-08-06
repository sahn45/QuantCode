#include "PCA.h"

PCA::PCA() {
	cout << "Running constuctor..." << endl;
}
PCA::~PCA() {
	cout << "Running destructor..." << endl;
}
Eigen::EigenSolver<Eigen::MatrixXf> PCA::run(const Eigen::MatrixXf &data, const float *target, const char *flag) {
	if (flag == "dim") {
		int dim = static_cast<int>(*target);
	}
	else if (flag == "pct") {

		Eigen::MatrixXf A = data;
		cout << "Random matrix: " << endl << A << endl;

		Eigen::EigenSolver<Eigen::MatrixXf> es(A);
		cout << "Eigenvalues are: " << endl << es.eigenvalues() << endl;
		cout << "Eigenvectors are : " << endl << es.eigenvectors() << endl;

		complex<float> lamda = es.eigenvalues()[0];
		cout << "Consider first eigenvalue, lamda = " << lamda << endl;
		Eigen::VectorXcf v = es.eigenvectors().col(0);
		cout << "If v is the corresponding eigenvector, then lamda * v = " << endl << lamda * v << endl;
		cout << "A * v = " << endl << A * v << endl;

		Eigen::MatrixXcf D = es.eigenvalues().asDiagonal();
		Eigen::MatrixXcf V = es.eigenvectors();
		cout << "Finally, V * D * V^(-1) = " << endl << V * D*V.inverse() << endl;

		return es;
	}

}
;