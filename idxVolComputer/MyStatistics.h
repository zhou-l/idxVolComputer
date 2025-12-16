#ifndef MY_STATISTICS_H_
#define MY_STATISTICS_H_
#include <vector>
#include <limits>
#include "../Eigen/Dense"
#include <iostream>


using std::cout;
using std::endl;
namespace MyStatistics
{
	void ComputeMean(const double* samples, const double* weights, const int num, const int dim, double** mean);
	void ComputeCovariance(const double* samples, const double* mean, const double* weights, const int num, int dim, double** cov);
	void ComputeCovariance(const std::vector<double>& samples, const double* mean, int dim, double* cov);
    double* randVec(int dim);
	// Using Eigen library
	void ComputeMean(const std::vector<Eigen::VectorXd>& samples, const double* weights, Eigen::VectorXd& mean);
	void ComputeCovariance(const std::vector<Eigen::VectorXd>& samples, Eigen::VectorXd mean, const double* weights, Eigen::MatrixXd& cov);
	void ComputeCovariance(const std::vector<Eigen::VectorXd>& samples, Eigen::VectorXd mean, Eigen::MatrixXd& cov);
	// How about template?
	template <class T>
	void ComputeMean(const std::vector<std::vector<T>>& samples, std::vector<T>& mean) // constant weights
	{
		if (samples.empty())
			return;
		int dim = int(samples[0].size());
		mean.resize(dim, T(0));
		for (size_t i = 0; i < samples.size(); i++)
		{
			for (int j = 0; j < dim; j++)
				mean[j] += samples[i][j];
		}
		T oneOverN = T(1.0 / double(samples.size()));
		for (int j = 0; j < dim; j++)
			mean[j] *= oneOverN;
	};

	// With VectorXd
	template <class T>
	void ComputeMean(const std::vector<std::vector<T>>& samples, Eigen::VectorXd& mean) // constant weights
	{
		if (samples.empty())
			return;
		int dim = int(samples[0].size());
		mean = Eigen::VectorXd(dim);
		mean.setConstant(0);
		for (size_t i = 0; i < samples.size(); i++)
		{
			for (int j = 0; j < dim; j++)
				mean(j) += double(samples[i][j]);
		}
		double oneOverN = 1.0 / double(samples.size());
		for (int j = 0; j < dim; j++)
			mean(j) *= oneOverN;
	};
	// With VectorXf
	template <class T>
	void ComputeMean(const std::vector<std::vector<T>>& samples, Eigen::VectorXf& mean) // constant weights
	{
		if (samples.empty())
			return;
		int dim = int(samples[0].size());
		mean = Eigen::VectorXf(dim);
		mean.setConstant(0);
		for (size_t i = 0; i < samples.size(); i++)
		{
			for (int j = 0; j < dim; j++)
				mean(j) += double(samples[i][j]);
		}
		double oneOverN = 1.0 / double(samples.size());
		for (int j = 0; j < dim; j++)
			mean(j) *= oneOverN;
	};
	// Weighted mean
	template <class T>
	void ComputeMean(const std::vector<std::vector<T>>& samples, const std::vector<T>& weights, std::vector<T>& mean) // constant weights
	{
		if (samples.empty())
			return;
		int dim = int(samples[0].size());
		mean.resize(dim, T(0));
		for (size_t i = 0; i < samples.size(); i++)
		{
			for (int j = 0; j < dim; j++)
				mean[j] += samples[i][j] * T(weights[i]);
		}
	};

	// covariance with MatrixXd
	template <class T>
	void ComputeCovariance(const std::vector<std::vector<T>>& samples, Eigen::VectorXd mean, Eigen::MatrixXd& cov)
	{
		if (samples.empty())
			return; // 
		int dim = mean.rows();
		cov = Eigen::MatrixXd(dim, dim);
		int num = int(samples.size());
		double oneOverNm1 = 1.0 / double(num -1);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				double s = 0.0;
				for (int k = 0; k < num; k++)
					s += (samples[k][j] - mean(j)) * (samples[k][i] - mean(i));
				s *= oneOverNm1;
				cov(i, j) = s;
				cov(j, i) = s;
			}
		}
	}; 

	// covariance with MatrixXf
	template <class T>
	void ComputeCovariance(const std::vector<std::vector<T>>& samples, Eigen::VectorXf mean, Eigen::MatrixXf& cov)
	{
		if (samples.empty())
			return; // 
		int dim = mean.rows();
		cov = Eigen::MatrixXf(dim, dim);
		int num = int(samples.size());
		double oneOverNm1 = 1.0 / double(num - 1);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				double s = 0.0;
				for (int k = 0; k < num; k++)
					s += (samples[k][j] - mean(j)) * (samples[k][i] - mean(i));
				s *= oneOverNm1;
				cov(i, j) = s;
				cov(j, i) = s;
			}
		}
	};

	
	//==========================================================================
	// Compute Correlations
	template <class T>
	double ComputeCorrelation(const std::vector<T>& X, const std::vector<T>& Y)
	{
		assert(!X.empty() && !Y.empty());
		assert(X.size() == Y.size());
		int N = int(X.size());
		// Use the simplified form
		double corr = 0.0;
		double sumXY = 0.0;
		double sumX = 0.0;
		double sumY = 0.0;
		double sumX2 = 0.0;
		double sumY2 = 0.0;
		for (int i = 0; i < N; i++)
		{
			double xd = double(X[i]);
			double yd = double(Y[i]);
			sumXY += xd * yd;
			sumX += xd;
			sumY += yd;
			sumX2 += xd * xd;
			sumY2 += yd * yd;
		}
		corr = (double(N) * sumXY - sumX * sumY) /
			(sqrt(double(N) * sumX2 - sumX * sumX) * sqrt(double(N) * sumY2 - sumY * sumY));
		return corr;
	}

	// Compute correlations from the covariance matrix
	template <class T>
	void ComputeCorrelationMat(const std::vector<std::vector<T>>& samples, const Eigen::VectorXf& mean, const Eigen::MatrixXf& cov, Eigen::MatrixXf& corr)
	{
		if (samples.empty())
			return; // 
		int dim = mean.rows();
		corr = cov;
		int num = int(samples.size());
		double delta = 1e-8;
		// 1. Compute standard deviation for each dimension
		Eigen::VectorXf stdev = mean;
		for (int i = 0; i < dim; i++)
		{
			stdev(i) = 0.0;
			for (int k = 0; k < num; k++)
			{
				stdev(i) += (samples[k][i] - mean(i)) * (samples[k][i] - mean(i));
			}
			stdev(i) = sqrt(1.0 / double(num) * stdev(i));
		}
		// 2. Compute correlations between every pair of dimensions using stdev and cov
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (stdev(i)*stdev(j) < delta) // If either stdev is very small, set the correlation to 0
					corr(i, j) = 0.0;
				else
					corr(i, j) = cov(i, j) / sqrt(stdev(i)*stdev(j));
			}
		}
	};
	//==========================================================================

	inline void EigenSolvCovMatrix(const std::vector<std::vector<float>>& valList, Eigen::VectorXf& eigVals, Eigen::MatrixXf& eigVecs)
	{
		Eigen::VectorXf mu;
		MyStatistics::ComputeMean<float>(valList, mu);
		//cout << "mu = " << endl<< mu << endl;
		Eigen::MatrixXf cov;
		MyStatistics::ComputeCovariance<float>(valList, mu, cov);
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(cov);
		// Get Eigenvals
		eigVals = eig.eigenvalues();
		// Get eigVecs
		eigVecs = eig.eigenvectors();
	}
	
	// Use Eigen::VectorXd to record multivariate value. 
	// Get eigen values, eigen vectors.
	// Also, return multivariate mean "mu" and correlation matrix as by-product
	inline void EigenSolvCovMatrix(const std::vector<Eigen::VectorXd>& valList, Eigen::VectorXd& mu, Eigen::VectorXd& eigVals, Eigen::MatrixXd& eigVecs)
	{
		assert(!valList.empty());
		// Compute the mean value
		for (std::vector<Eigen::VectorXd>::const_iterator IT = valList.begin(); IT != valList.end(); ++IT)
		{
			if (IT == valList.begin())
				mu.setZero();
			mu += *IT;
		}
		mu /= double(valList.size());
		// compute covariance 
		Eigen::MatrixXd cov;
		MyStatistics::ComputeCovariance(valList, mu, cov);
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
		// Get Eigenvals
		eigVals = eig.eigenvalues();
		// Get eigVecs
		eigVecs = eig.eigenvectors();
		
	}


	// Also outputs correlation information
	inline void EigenSolvCovMatrixGetCorr(const std::vector<std::vector<float>>& valList, Eigen::VectorXf& eigVals, Eigen::MatrixXf& eigVecs, Eigen::MatrixXf& corr)
	{
		Eigen::VectorXf mu;
		MyStatistics::ComputeMean<float>(valList, mu);
		//cout << "mu = " << endl<< mu << endl;
		Eigen::MatrixXf cov;
		MyStatistics::ComputeCovariance<float>(valList, mu, cov);
		// Compute correlations
		MyStatistics::ComputeCorrelationMat<float>(valList, mu, cov, corr);
		// Use eigen solver
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(cov);
		// Get Eigenvals
		eigVals = eig.eigenvalues();
		// Get eigVecs
		eigVecs = eig.eigenvectors();
	}



	inline void SVDCovMatrix(const std::vector<std::vector<float>>& valList, Eigen::MatrixXf& U, Eigen::VectorXf& S, Eigen::MatrixXf& V)
	{
		// NOTE: for Eigen SVD decomposition
		// http://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
		Eigen::VectorXf mu;
		MyStatistics::ComputeMean<float>(valList, mu);
		//cout << "mu = " << endl<< mu << endl;
		Eigen::MatrixXf cov;
		MyStatistics::ComputeCovariance<float>(valList, mu, cov);
		//cout << "Cov = " << endl<< cov << endl;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
		// Left singular vectors
		//cout << "U= " <<endl<< svd.matrixU() << endl;
		// Get the major eigen vector
		//cout << "Eigen vectors = " << endl;

		////////////////////
		// SVD gives us A = USV(^T), where
		// U: each row of U is an eigen vector
		// S: each diagonal element is the eigen value

		int maxEigId = -1;
		float maxEigVal = -std::numeric_limits<float>::max();

		// NOTE: no need to sort eigen values!!!
		// Eigen make sure singular values are always sorted in decreasing order!!!
		//For the SVD decomposition of a n - by - p matrix, letting m be the minimum of n and p, the returned vector has size m.Singular values are always sorted in decreasing order.

		// For test, take the first 2-component of the major eigen vector
		U = svd.matrixU();
		S = svd.singularValues();
		V = svd.matrixV();
	};

	//--------------------------------------------------------
	// Probablisitic PCA!

	// EM algorithm to derive W
	// input: S - covariance of input data Y
	//        k - number of PC components
	// output: W - ppca coefficients
	//         v - isotropic error
	void em_ppca(const Eigen::MatrixXf& S, int k, Eigen::MatrixXf& W, float& v);
	// Main function
	// Input:  Y - input p-dimensional data of N elements; N * p matrix
	//         k - number of PC components
	// Output: coeff - PCA coefficients; p * k matrix
	//         score - PCA scores;       N * k matrix // not used for now
	//         pcvar - PC variances;     eigenvalues
	//         mu    - estimated mean
	//         v     - isotropic noise; scalar
	void ppca(const Eigen::MatrixXf& Y, int k, Eigen::MatrixXf& coeff, Eigen::MatrixXf& score, Eigen::VectorXf& pcvar, Eigen::VectorXf& mu, float& v);	
	//-----------------------------------------------------------
}


// Eigen version Gaussian distribution
class GaussianXd
{
public:

	GaussianXd(Eigen::VectorXd mu, Eigen::MatrixXd cov);
	~GaussianXd();

	void Update(Eigen::VectorXd mu, Eigen::MatrixXd cov);
	double ProbFunc(Eigen::VectorXd x);
	void Fit(const std::vector<Eigen::VectorXd>& samples, double* weights, const int num);

	Eigen::VectorXd      m_mu;
	Eigen::MatrixXd      m_cov;
	int           m_dim;
};

//// A template version of the multivariat Gaussian inside the Eigen namespace
//namespace Eigen{
//	template <typename T>
//	class GaussianX
//	{
//		Matrix<T, Dynamic, Dynamic> m_cov;
//		Matrix<T, Dynamic, 1>        m_mu;
//		int                         m_dim;
//	public:
//		GaussianX(Matrix<T, Dynamic, 1> mu, Matrix<T, Dynamic, Dynamic> cov) :
//			m_mu(mu),
//			m_cov(cov)
//		{
//			m_dim = int(mu.SizeAtCompileTime());
//		}
//
//		T ProbFunc(Matrix<T, Dynamic, 1> x)
//		{
//			Matrix<T, Dynamic, Dynamic> invCov = m_cov.inverse();
//			double a = pow(2.0 * M_PI, m_dim / 2.0) * sqrt(m_cov.determinant());
//			Matrix<T, Dynamic, 1> dif(m_dim);
//			dif = x - m_mu;
//			double res = dif.transpose() * invCov * dif;
//			double N = exp(-0.5 * res);
//			N *= 1.0 / a;
//			return T(N);
//		}
//	};
//
//}

extern void drawSamplesFromGaussianXd(const GaussianXd& distr, int num_samples, Eigen::MatrixXd& samples); // Draw "num_samples" samples from the Gaussian distribution "distr"
#endif