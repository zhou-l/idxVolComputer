#include <cstdlib>
#include <memory>
#include <cmath>
#include "StdDefines.h"
#include "MersenneTwister.h"
#include "MyStatistics.h"
#include "eigenmvn.h"
using namespace std;
using namespace Eigen;
namespace MyStatistics
{
	void ComputeMean(const double* samples, const double* weights, const int num, const int dim, double** mean)
	{
		SAFE_DELETE_ARRAY(*mean);
		*mean = new double[dim];
		memset(*mean, 0, sizeof(double) * dim);
		for(int i = 0;  i < num; i++ )
		{
			double w = weights[i];
			for(int j = 0; j < dim; j++ )
				(*mean)[j] += samples[i * dim + j] * w;
		}
	}

	void ComputeCovariance(const double* samples, const double* mean, const double* weights, const int num, int dim, double** cov)
	{
		double sw = 1.0;
		for(int i = 0; i < num; i++ )
			sw -= weights[i] * weights[i];

		SAFE_DELETE_ARRAY(*cov);
		*cov = new double[dim * dim];
		memset(*cov, 0, sizeof(double) * dim * dim);
		for(int i = 0; i < dim; i++ )
		{
			for(int j = i; j < dim; j++ )
			{
				double s = 0.0;
				for(int k = 0; k < num; k++ )
					s += weights[k] * (samples[k * dim + j] - mean[j]) * (samples[k * dim + i] - mean[i]);
				s /= sw;
				(*cov)[i * dim + j] = s;
				(*cov)[j * dim + i] = s;
			}
		}
	}

	void ComputeCovariance(const vector<double>& samples, const double* mean, int dim, double* cov)
	{
        size_t totalNum = samples.size() / dim;
		for(int i = 0; i < dim; i++ )
		{
			for(int j = i; j < dim; j++ )
			{
				double s = 0.0;
				for(size_t k = 0; k < totalNum; k++ )
					s += (samples[k * dim + j] - mean[j]) * (samples[k * dim + i] - mean[i]);
				s /= ((double)totalNum - 1);
				cov[i * dim + j] = s;
				cov[j * dim + i] = s;
			}
		}
	}


	double* randVec(int dim)
	{
		MTRand rng;
		rng.seed((unsigned)0);
		double* rv = new double[dim];
		for(int i = 0; i < dim; i++ )
		    rv[i] = rng.rand();
		return rv;
	}

	void ComputeMean(const vector<VectorXd>& samples, const double* weights, VectorXd& mean)
	{
		mean.setConstant(0);
		size_t num = samples.size();
		for(size_t i = 0; i < num; i++ )
		{
			double w = weights[i];
			mean += samples[i] * w;
		}
	}

	void ComputeCovariance(const vector<VectorXd>& samples, VectorXd mean, const double* weights, MatrixXd& cov)
	{
		cov.setConstant(0);
		int dim = int(mean.size());
		double sw = 1.0;
		size_t num = samples.size();
		for(size_t i = 0; i < num; i++ )
			sw -= weights[i] * weights[i];

		for(int i = 0; i < dim; i++ )
		{
			for(int j = i; j < dim; j++ )
			{
				double s = 0.0;
				for(int k = 0; k < num; k++ )
					s += weights[k] * (samples[k](j) - mean(j)) * (samples[k](i) - mean(i));
				s /= sw;
				cov(i,j) = s;
				cov(j,i) = s;
			}
		}
	}

	void ComputeCovariance(const vector<VectorXd>& samples, VectorXd mean,  MatrixXd& cov)
	{
		cov = MatrixXd::Zero(mean.size(), mean.size());
		int dim = int(mean.size());
		double sw = 1.0;
		size_t num = samples.size();
		sw = 1.0 / double(num);

		for (int i = 0; i < dim; i++)
		{
			for (int j = i; j < dim; j++)
			{
				double s = 0.0;
				for (int k = 0; k < num; k++)
					s += (samples[k](j) - mean(j)) * (samples[k](i) - mean(i));
				s *= sw;
				cov(i, j) = s;
				cov(j, i) = s;
			}
		}
	}

	// ppca functions
	void em_ppca(const Eigen::MatrixXf& S, int k, Eigen::MatrixXf& W, float& v)
	{
		int d = S.cols();
		// Init
		W = Eigen::MatrixXf::Ones(d, k);
		v = 1.0f;
		float epsilon = 0.01f;
		// EM loop
		Eigen::MatrixXf M;
		Eigen::MatrixXf I = Eigen::MatrixXf::Identity(k, k);
		for (;;)
		{
			M = W.transpose() * W + v * v * I;
			Eigen::MatrixXf W_new = S * W * (v * v * I + M.inverse() * W.transpose() * S * W).inverse();
			float v_new = sqrt(1.0f / float(d) * (S - S * W * M.inverse() * W_new.transpose()).trace());
			Eigen::MatrixXf Wdiff = (W_new - W).cwiseAbs();
			if (abs(v_new - v) < epsilon && Wdiff.maxCoeff() < epsilon)
			{
				W = W_new;
				v = v_new;
				break;
			}
			W = W_new;
			v = v_new;
		}
	}

	void ppca(const Eigen::MatrixXf& Y, int k, Eigen::MatrixXf& coeff, Eigen::MatrixXf& score, Eigen::VectorXf& pcvar, Eigen::VectorXf& mu, float& v)
	{
		// dimension
		int d = Y.cols();
		// number of data
		int N = Y.rows();
		cout << "d = " << d << ", N = " << N << endl;
		// compute sample mean
		mu = Y.colwise().mean(); // get row vector

		cout << "Mu of Y=" << mu << endl<<"dim = " <<mu.cols()<<","<<mu.rows()<<endl;
		// compute covariance of the input data
		Eigen::MatrixXf S = Eigen::MatrixXf::Zero(d, d);
		for (int n = 0; n < N; n++)
		{
			//cout << "Y row dim = " << Y.row(n).cols() << "," << Y.row(n).rows() << endl;
			VectorXf t = Y.row(n).transpose() - mu;
			S = S + (t) * t.transpose();
		}
		S = 1.0f / float(N) * S;
		// EM estimation of W
		Eigen::MatrixXf W;
		em_ppca(S, k, W, v);
		cout << "W=" << W << endl;
		// Get X
		Eigen::MatrixXf WTW = W.transpose() * W;
		cout << "WTW=" << WTW << endl;
		Eigen::MatrixXf I = Eigen::MatrixXf::Identity(k, k);
		Eigen::MatrixXf M = WTW + v*v * I;
		Eigen::MatrixXf Ynorm = Y.rowwise() - mu.transpose();
		Eigen::MatrixXf X = M.inverse() * W.transpose() * Ynorm.transpose();
		cout << "X dim = " << X.rows()<<","<<X.cols()<<endl<<" val = "<< X << endl;
		Eigen::MatrixXf Y_hat = W * WTW.inverse() * (WTW + v * I) * X;

		// Orthogonalization
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(W, ComputeThinU | ComputeThinV);
		coeff = svd.matrixU();
		cout << "coeff=" << coeff << endl;
		score = Y_hat.transpose() * coeff; // the projected data in k-D space

		// Compute the eigen values
		Eigen::MatrixXf T = score.transpose() * score;
		Eigen::VectorXcf pcvar_complex = T.eigenvalues();
		cout << "pcvar complex = " << pcvar_complex << endl;

		pcvar = pcvar_complex.real();
		pcvar = pcvar / float(N - 1);
		cout << "pcvar = " << pcvar << endl;
	}

}

GaussianXd::GaussianXd(VectorXd mu, MatrixXd cov) :
m_mu(mu),
m_cov(cov)
{
	m_dim = m_mu.size();
}

GaussianXd::~GaussianXd()
{}

void GaussianXd::Update(VectorXd mu, MatrixXd cov)
{
	m_mu = mu;
	m_cov = cov;
}

void GaussianXd::Fit(const vector<VectorXd>& samples, double* weights, const int num)
{
	MyStatistics::ComputeMean(samples, weights, m_mu);
	MyStatistics::ComputeCovariance(samples, m_mu, weights, m_cov);
}

double GaussianXd::ProbFunc(VectorXd x)
{
	MatrixXd invCov = m_cov.inverse();
	double a = pow(2.0 * M_PI, m_dim / 2.0) * sqrt(m_cov.determinant());
	VectorXd dif(m_dim);
	dif = x - m_mu;
	double res = dif.transpose() * invCov * dif;
	double N = exp(-0.5 * res);
	N *= 1.0 / a;
	return N;
}

void drawSamplesFromGaussianXd(const GaussianXd& distr, int num_samples, MatrixXd& samples)
{
	// Use the eigenmvn routine
	EigenMultivariateNormal<double> normX_solver(distr.m_mu, distr.m_cov);
	samples = normX_solver.samples(num_samples);
}