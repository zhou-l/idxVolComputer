#include "XPCPdataAnalyzer.h"
#include <chrono>
#include <iomanip> // time
#include <iostream>
#include <fstream>
#include "ANN/ANN.h"
#include "MyStatistics.h"
#include "PCPInselberg.h"
#include <set>
#include "KDtree.h"
#include "pflatPt.h"
#include "UncertainXPCPSample.h"
#include "qsp.h"
#include "iqsort.h"
#include "VolumeData.h"
#include "PointMonteCarloSampler.h"
#include <QString>
#include "SysTools.h"
#include "progressBar.hpp"
#define OUTPUT_FILE 0 //1
#define IDX_PT_DEBUG 1

using namespace std;
using namespace Eigen;

Eigen::MatrixXf corr_mat; // a handy correlation matrix

vector<double> scaleLocs;
vector<double> scaleFactors;
const string defautSFfileName = "MirrorScaleRvL.txt";

// Polynomial fitted function of the scaling factors
const vector<double> polyScalingParams =
{
	//-0.738940869888870,	1.82144491769081, - 1.68803863213442,	1.53970674544351
	//0.443591625026322, - 1.08188255732298,	1.71564129882120, - 1.63277698493305,	1.54498792598118
	//0.217910102973902, -0.853227062338163, 1.77345046206235, -1.66980305421205, 1.54197728632631

	- 0.635419295094657,	1.27083859018931,	0.0598250374937986, - 0.695244332588451,	0.999315774099720
	//0.688400063140494 - 1.08298317867892	1.52898163124705 - 1.66056969842561	1.55468174727895
};
// u of the maximum vertical location after transformation
const double u_vertical_max = 0.778;

// Spline fit functions
const vector<vector<double>> spline_coefs =
{
	//{ -0.589786534624352, 2.64359160725805, -3.91086105722194, 2.90941000000000 },
	//{ -0.589786534624355, 1.78941733425189, -1.77078609362754, 1.57115529400000 },
	//{ -0.589786534624356, 0.444139947542413, -0.0725704813004136, 1 }
	// Coeffs 1
	//{0.508887654273929,	0.349519650248899, -0.731185390954345,	1},
	//{0.508887654273929,	0.731185390954345, -0.461009130653534,	0.847000000000000},
	//{-0.508887654273930,	1.11285113165979,	0,	0.785398173869293},
	//{-0.508887654273927,	0.731185390954345,	0.461009130653534,	0.847000000000000}
//// Coeffs 2
//{-28.5631865147271,	43.1215898860453, -22.4199983143409,	5},
//{-28.5631865147271,	21.6992000000000, -6.21480084282956,	1.64380000000000},
//{0.633532573635441,	0.276810113954683, -0.720798314340886,	1},
//{0.467339347820090,	0.751959544181265, -0.463605899806899,	0.847000000000000},
//{-0.467339347820088,	1.10246405504633, 0,	0.785398173869293},
//{-0.633532573635454,	0.751959544181266,	0.463605899806899,	0.847000000000000},
//{28.5631865147271,	0.276810113954680,	0.720798314340885,	1},
//{28.5631865147271,	21.6992000000000,	6.21480084282956,	1.64380000000000}

// Coeffs 3
//	{-28.6084129032258,	43.1555096774194, - 22.4256516129032,	5},
//{-28.6084129032258,	21.6992000000000 ,-6.21197419354839,	1.64380000000000},
//{0.859664516129035,	0.242890322580644, -0.726451612903226,	1},
//{8.88178419700125e-16,	0.887638709677419, - 0.443819354838710,	0.847000000000000},
//{-0.859664516129039,	0.887638709677420,	0.443819354838710,	0.847000000000000},
//{28.6084129032258,	0.242890322580651,	0.726451612903225,	1},
//{28.6084129032258,	21.6992000000000,	6.21197419354839,	1.64380000000000}

// Coeffs 4
//{14.3335225806452, -14.6509419354839,	2.14209032258065,	1.80000000000000},
//{14.3335225806451, -3.90079999999999, -2.49584516129032,	1.64380000000000},
//{-9.05001290322580,	6.84934193548387, -1.75870967741935,	1},
//{2.08166817117217e-16,	0.0618322580645168 - 0.0309161290322584,	0.847000000000000},
//{9.05001290322580,	0.0618322580645172,	0.0309161290322585,	0.847000000000000},
//{-14.3335225806452,	6.84934193548387,	1.75870967741935,	1},
//{-14.3335225806452, -3.90080000000000,	2.49584516129032,	1.64380000000000}

	// coeffs 5: -0.5, -0.25 0 0.25 0.75 1 1.25 1.5
//{-2.84283870967742,	2.13212903225807, -0.355354838709678,	1},
//{-2.84283870967742,	0,	0.177677419354839,	1},
//{4.42219354838710, -2.13212903225806, -0.355354838709678,	1},
//{0,	1.18451612903226, -0.592258064516129,	0.847000000000000},
//{-4.42219354838710,	1.18451612903226,	0.592258064516129,	0.847000000000000},
//{2.84283870967742, -2.13212903225807,	0.355354838709678,	1},
//{2.84283870967742,	0, -0.177677419354839,	1}

// coeffs 6: vertical scale = [1.2 1.1 1 0.847  0.847 1 1.1  1.2] for u=[-0.5, -0.25 0 0.25 0.75 1 1.25 1.5]
//	{-1.191225806,0.893419355,-0.548903226,1.2 },
//{-1.191225806,-1.78E-15,-0.325548387,1.1},
//{2.564129032,-0.893419355,-0.548903226,1},
//{-4.44e-16,1.029677419,-0.51483871,0.847},
//{-2.564129032,1.029677419,0.51483871,0.847},
//{1.191225806,-0.893419355,0.548903226,1},
//{1.191225806,-1.78e-15,0.325548387,1.1}

	// best known coeffs so far
//{-0.315870968,0.236903226,-0.651483871,1.306},
//{-0.315870968,-4.44e-16,-0.592258065,1.153},
//{1.579354839,-0.236903226,-0.651483871,1},
//{-2.22e-16,0.947612903,-0.473806452,0.847},
//{-1.579354839,0.947612903,0.473806452,0.847},
//{0.315870968,-0.236903226,0.651483871,1},
//{0.315870968,0.00E+00,0.592258065,1.153}

	// best coeffs so far
//{-0.614549707602338, 0.460912280701753, -0.688818713450292, 1.30600000000000},
//{-0.614549707602338, 0, -0.573590643274854, 1.15300000000000},
//{4.69099415204676, -0.460912280701750, -0.688818713450292, 1},
//{-0.282365172189726, 0.946385964912279, -0.640271345029240, 0.931200000000000},
//{-0.182643274853802, 0.819321637426901, -0.375415204678363, 0.855500000000000},
//{0.182643274853802, 0.682339181286549, 0, 0.810000000000000},
//{0.282365172189731, 0.819321637426900, 0.375415204678363, 0.855500000000000},
//{-4.69099415204678, 0.946385964912282, 0.640271345029239, 0.931200000000000},
//{0.614549707602336, -0.460912280701753, 0.688818713450292, 1},
//{0.614549707602341, 0, 0.573590643274854, 1.15300000000000}

//{-0.606128654970759, 0.454596491228069, -0.687766081871345, 1.30600000000000},
//{-0.606128654970759, 0, -0.574116959064328, 1.15300000000000},
//{4.52257309941520, -0.454596491228068, -0.687766081871345, 1},
//{0.134009096816115, 0.902175438596489, -0.643008187134503, 0.931200000000000},
//{-0.628959064327477, 0.962479532163740, -0.363309941520468, 0.855500000000000},
//{0.628959064327478, 0.490760233918132, -1.38777878078145e-17, 0.815000000000000},
//{-0.134009096816097, 0.962479532163737, 0.363309941520468, 0.855500000000000},
//{-4.52257309941521, 0.902175438596492, 0.643008187134503, 0.931200000000000},
//{0.606128654970760, -0.454596491228069, 0.687766081871345, 1},
//{0.606128654970757, 0, 0.574116959064328, 1.15300000000000}

{-0.611181286549707, 0.458385964912280, -0.688397660818713, 1.30600000000000},
{-0.611181286549707, 4.44089209850063e-16, -0.573801169590644, 1.15300000000000},
{4.62362573099413, -0.458385964912277, -0.688397660818713, 1},
{-0.115815464587385, 0.928701754385963, -0.641366081871345, 0.931200000000000},
{-0.361169590643276, 0.876584795321637, -0.370573099415205, 0.855500000000000},
{0.361169590643276, 0.605707602339181, -2.77555756156289e-17, 0.812000000000000},
{0.115815464587392, 0.876584795321637, 0.370573099415205, 0.855500000000000},
{-4.62362573099415, 0.928701754385967, 0.641366081871345, 0.931200000000000},
{0.611181286549705, -0.458385964912280, 0.688397660818714, 1},
{0.611181286549707, 4.44089209850063e-16, 0.573801169590643, 1.15300000000000}
};
const vector<double> spline_breaks =
{
	//-0.4828, 0, 0.7603, 1
	//0,	0.250000000000000,	0.500000000000000,	0.750000000000000,	1
	//-0.5, - 0.25,	0,	0.25,	0.5,	0.75,	1,	1.25,	1.5
	//-0.5, -0.25,	0,	0.25,	0.75,	1,	1.25,	1.5
	-0.5, -0.25,	0, 0.1,	0.25,	0.5, 0.75,	0.9, 1, 1.25,	1.5
};


XPCPdataAnalyzer::XPCPdataAnalyzer() :
//_kdTree_rawData(NULL),
//_kdTree_xpcp(NULL),
//_kdTree_1flat(NULL),
//_kdTree_2flat(NULL),
_pcp(g_pcp), // Use global PCP
_dimRawData(0),
_dimXPCPdata(0),
_numSamples(0),
_annKdTree_rawData(NULL)
{
}


XPCPdataAnalyzer::~XPCPdataAnalyzer()
{
	//SAFE_DELETE(_kdTree_rawData);
	//SAFE_DELETE(_kdTree_xpcp);
	//SAFE_DELETE(_kdTree_1flat);
	//SAFE_DELETE(_kdTree_2flat);
	_pcp = NULL;
	annClose();
}

void XPCPdataAnalyzer::buildKDTree(const vector<vector<float>>& pointCenterData, ANNkd_tree** pKdTree)
{
	SAFE_DELETE(*pKdTree);
	ANNpointArray dataPts;
	UINT64 nPts = pointCenterData.size();
	UINT64 dim = pointCenterData[0].size();
	// Setup ANN points
	dataPts = annAllocPts(nPts, dim);
	for (size_t i = 0; i < pointCenterData.size(); i++)
	{
		for (size_t j = 0; j < pointCenterData[i].size(); j++)
			dataPts[i][j] = double(pointCenterData[i][j]);
	}
	// build kd tree
	*pKdTree = new ANNkd_tree(dataPts,
		nPts,
		dim);

}

// The function for computing indexed points for spatial neighborhoods
bool XPCPdataAnalyzer::analyzeSpatialRawDataContIdxPts(const std::vector<VolumeData*>& volList, const std::vector<std::vector<float>>& samplePos, UINT64VECTOR3 idxVolDim, int numSamplesInNeighborhood, FLOATVECTOR3 neighborRadius, bool writeOutFile)
{
	if (volList.empty())
		return false;
	if (volList.size() < 3)
	{
		cout << "The number of volumes is less than 3 and 2-flats cannot be calculated! Failed to compute spatial domain indexed points." << endl;
		return false;
	}

	VolumeData* vol1 = volList[0];
	if (vol1 == NULL)
		return false;


	UINT64VECTOR3 volSize = vol1->getDim();
	FLOATVECTOR3 diskR = neighborRadius;

	UINT64 num_nb = numSamplesInNeighborhood;
	UINT64 dim_rawData = volList.size();

	_numSamples = UINT64(samplePos.size());
	_dimRawData = dim_rawData; // Raw data dimension
	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension

	//Preperations: Clear list
	// Create a temp vector for the mu value of each sample pos
	std::vector<Eigen::VectorXf> muVol(samplePos.size());

	g_xpcpData.resize(samplePos.size(), XPCPSample(dim_rawData));
	// build pcp if necessary
	if (_pcp == NULL)
	{
		vector<string> attribs(dim_rawData);
		for (size_t i = 0; i < attribs.size(); i++)
			attribs[i] = string("attr_") + number2String(int(i));
		_pcp = new PCPInselberg(attribs);
	}
	// Major eigenvectors
	g_majEigData.resize(samplePos.size());
	// NOTE: Record original eigenvector data as volumes 
	// Unnormalized eigenvector data!
	std::vector<Eigen::VectorXf> OrgMajEigData(samplePos.size());
	std::vector<Eigen::VectorXf> OrgSecEigData(samplePos.size());
	// use global variables to keep the p-flat records
	g_1flat_list.clear();

	ofstream ofDebug("neighborFitTest.txt");
	ofstream ofInterpEigDebug("eigVecInterp.csv");

	// Prepare variables to record min/max of p-flat strength in each subspace
	_1flatsMinMaxPerSubspace.clear();
	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
	{
		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	_2flatsMinMaxPerSubspace.clear();
	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
	{
		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	// Process each data point
	vector<md_val> pointCenterData(samplePos.size());
	string dummy;
	PointMonteCarloSampler* sampler = new PointMonteCarloSampler(volList, volSize, dummy);


	PR_METHOD method = PR_SUBSPACE;

	// Auxilary stuff
	size_t dumpPos = samplePos.size() / 2;
	progressbar bar(100);
	size_t stepsOnePercent = samplePos.size() / 100;

	// Storing eigenvalues
	vector<Eigen::VectorXf> eigValData;
	
	for (vector<vector<float>>::const_iterator IT = samplePos.begin(); IT != samplePos.end(); ++IT)
	{
		UINT64 i = IT - samplePos.begin();
		FLOATVECTOR3 normPos = FLOATVECTOR3((*IT)[0], (*IT)[1], (*IT)[2]);
		// sample the center 
		pointCenterData[i] = sampler->sampleVoxelWithInterpNormCoords(normPos);
		FLOATVECTOR3 pos = FLOATVECTOR3((*IT)[0] * float(volSize.x-1), 
			(*IT)[1] * float(volSize.y - 1), (*IT)[2] * float(volSize.z - 1)) ; // convert the position back to the volume's dimension
		// sample the neighborhood: Get all points for computation

		vector<md_val> neighbors(num_nb);
		UINT64 cnt = 0;
		while(cnt < num_nb)
		{
			FLOATVECTOR3 neighborLoc = FLOATVECTOR3(g_randGen.rand(), g_randGen.rand(), g_randGen.rand());
			neighborLoc = 2.0f * neighborLoc - FLOATVECTOR3(1.0f, 1.0f, 1.0f);
			neighborLoc = pos + FLOATVECTOR3(diskR.x * neighborLoc.x, diskR.y * neighborLoc.y, diskR.z * neighborLoc.z);

			if (neighborLoc.x < 0.0f || neighborLoc.y < 0.0f || neighborLoc.z < 0.0f ||
				neighborLoc.x >= float(volSize.x) || neighborLoc.y >= float(volSize.y) || neighborLoc.z >= float(volSize.z))
				continue; 
			//sample
			neighborLoc = FLOATVECTOR3(neighborLoc.x / (volSize.x - 1), 
				neighborLoc.y / (volSize.y - 1), neighborLoc.z / (volSize.z-1));
			neighborLoc.x = MAX(0.0f, MIN(1.0f, neighborLoc.x));
			neighborLoc.y = MAX(0.0f, MIN(1.0f, neighborLoc.y));
			neighborLoc.z = MAX(0.0f, MIN(1.0f, neighborLoc.z));
			md_val neighborVal = sampler->sampleVoxelWithInterpNormCoords(neighborLoc);
			neighbors[cnt++]= neighborVal;
		}


#if OUTPUT_FILE>0
			if (i == dumpPos) // check the central position
			{
				ofstream ofNeighbor0("neighborTest0.txt");
				ofNeighbor0 << "dumpPos pos= " << pos << endl;
				for (int i = 0; i < neighbors.size(); i++)
				{
					for (size_t d = 0; d < dim_rawData - 1; d++)
						ofNeighbor0 << neighbors[i][d] << ",";
					ofNeighbor0 << neighbors[i][dim_rawData - 1] << endl;
				}
				ofNeighbor0.close();
			}
#endif
			// Compute mean vector
			Eigen::VectorXf mu;
			MyStatistics::ComputeMean<float>(neighbors, mu);
#if OUTPUT_FILE>0
			if (i == dumpPos)
			{
				//cout << "Mean = " << endl << mu << endl;
				ofDebug << "Mean = " << endl << mu << endl;
			}
#endif
			// The correlation matrix recording correlations between every pair of attributes
			Eigen::MatrixXf corr;
			corr_mat = Eigen::MatrixXf::Zero(dim_rawData, dim_rawData);

			XPCPSample xpcp_tuple(dim_rawData);
			Eigen::VectorXf majEigV;
			Eigen::VectorXf secMajEigV;


			Eigen::VectorXf eigVal;
			Eigen::MatrixXf eigVec;
			//MyStatistics::EigenSolvCovMatrix(neighbors, eigVal, eigVec);
			MyStatistics::EigenSolvCovMatrixGetCorr(neighbors, eigVal, eigVec, corr);
#if OUTPUT_FILE > 0
			if (i == dumpPos)
			{
				//cout << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
				ofDebug << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
				// output correlations
		/*		cout << "Correlations: " << endl;
				for (int dy = 0; dy < corr.rows() - 1; dy++)
				{
					int dx = dy + 1;
					cout << "(" << dy << "," << dx << ") = ";
					cout << corr(dy, dx) << endl;
				}*/
			}
#endif

			majEigV = eigVec.col(eigVec.cols() - 1);
			// TODO: do we save the normalized eigenvectors or the originals???
			OrgMajEigData[i] = majEigV;
			// Do nomralization
			majEigV.normalize();

			secMajEigV = eigVec.col(eigVec.cols() - 2);
			OrgSecEigData[i] = secMajEigV;
			// TODO: do we save the normalized eigenvectors or the originals???
			secMajEigV.normalize();
#if OUTPUT_FILE > 0
			if (i == dumpPos)
			{
				//cout << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
				ofDebug << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
			}
#endif
			corr_mat = corr;
			// Save eigenvalues
			eigValData.push_back(eigVal);

			// compute 1-flat and 2-flat and record min/max of the strength of p-flats
			if (method == PR_SUBSPACE)
				calc_p_flats_subSpaceMethod(mu, majEigV, secMajEigV, xpcp_tuple);
			else if (method == PR_FROM_LOWER_DIM)
				calc_p_flats_fromLowerDimMethod(mu, majEigV, secMajEigV, xpcp_tuple);
			else
				calc_p_flats_subSpaceMethod(mu, majEigV, secMajEigV, xpcp_tuple);
	

			g_majEigData[i] = majEigV; // Set normalized Eigen vector data

			g_xpcpData[i] = xpcp_tuple;

			muVol[i] = mu;
			// Add a progress bar
			

			if (i % stepsOnePercent == 0)
			{
				// Update the progress bar
				bar.update();
				//cout << "progress: " << float(i) / float(samplePos.size()) * 100.0f <<"%"<<endl;
				
				//if (i == 0)
				//{
				//	for (size_t ii = 0; ii < majEigV.size(); ii++)
				//	{
				//		ofInterpEigDebug << "EigVec1_" << number2String(ii) << ",";
				//	}

				//	for (size_t ii = 0; ii < majEigV.size(); ii++)
				//	{
				//		ofInterpEigDebug << "EigVec2_" << number2String(ii) << ",";
				//	}
				//	for (size_t ii = 0; ii < majEigV.size(); ii++)
				//	{
				//		ofInterpEigDebug << "InterpEigVec_" << number2String(ii) << ",";
				//	}
				//	for (size_t ii = 0; ii < majEigV.size(); ii++)
				//	{
				//		ofInterpEigDebug << "MidVec_" << number2String(ii) << ",";
				//	}
				//	ofInterpEigDebug << "EigVal1, InterpEigVal, EigVal2, MidEigVal, InterpMidEigDotProd" << endl;
				//}
				//else 
				//{
				//	// test eigenvector interpolation
				//	Eigen::VectorXf vt;
				//	float t = 0.5f;
				//	Eigen::VectorXf v1 = majEigV;
				//	Eigen::VectorXf v2 = g_majEigData[i - 2];
				//	eigenVecInterp1D(vt, v1, v2, t);
				//	vector<float> eigVals;
				//	eigVals.push_back(eigVal[0]);
				//	eigVals.push_back(eigValData[i - 2][0]);
				//	vector<float> ts;
				//	ts.push_back(1 - t);
				//	ts.push_back(t);

				//	float lamt = 0.0f;
				//	eigenValInterpND(lamt, eigVals, ts);
				//	float diff = g_majEigData[i - 1].dot(vt);

				//	for (size_t ii = 0; ii < v1.size(); ii++)
				//		ofInterpEigDebug << v1(ii) << ",";
				//	for (size_t ii = 0; ii < v2.size(); ii++)
				//		ofInterpEigDebug << v2(ii) << ",";
				//	for (size_t ii = 0; ii < vt.size(); ii++)
				//		ofInterpEigDebug << vt(ii) << ",";
				//	for (size_t ii = 0; ii < g_majEigData[i - 1].size(); ii++)
				//		ofInterpEigDebug << g_majEigData[i - 1](ii) << ",";
				//	ofInterpEigDebug << eigVals[0] << "," << lamt << "," << eigVals[1] << "," << eigValData[i - 1][0] << diff<< endl;

				//	// Test interpvec and the midvec angle difference
	
				//	cout << "Angle diff (dot product)=" << diff << endl;
				//}
			}

			///////////////////////////////////////////////////////////////////////
			// !1 TODO: If we turn off normalization (!2), uncomment codes in the block
			//// Convert p-flat data and store to global variables
			//vector<pFlatPt*> oneFlats;
			//xpcp_tuple.to_pFlatRec(oneFlats, 1, UINT64(i));
			//// Add to global list
			//g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
			//vector<pFlatPt*> twoFlats;
			//xpcp_tuple.to_pFlatRec(twoFlats, 2, UINT64(i));
			//// Add to global list
			//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
			////////////////////////////////////////////////////////////////////////
		}
		SAFE_DELETE(sampler);
		ofDebug.close(); 
		ofInterpEigDebug.close();
		cout << endl << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
		
		/////////////////////////////////////////////////////
		// !2 Sep 22, 2016: Added normalization func for p-flats 
		// Do normalization after we compute all p-flats!
		int normMethod = 0;
		nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
		// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization

#ifdef IDX_PT_DEBUG 
		string debugFileName = g_attrib_names[0];
		//string idxDebugFileName = "idxPtDebug_compute";
		stringstream ssIdxDebugFileName; 
		ssIdxDebugFileName << "idxPtDebug_compute"<<debugFileName<<".txt";
		
		stringstream ss2flatIdxDebugFileName;
		ss2flatIdxDebugFileName << "twoFlats_idxPtDebug_compute" << debugFileName << ".txt";
		// show the concated file names
		cout << ssIdxDebugFileName.str() << ","<<ss2flatIdxDebugFileName.str()<< endl;
	    //****************************************
		// !!!NOTE:
		// Debug files are used for continuous indexed points in Matlab
		ofstream ofIdxDebug(ssIdxDebugFileName.str()); // ("idxPtDebug_compute.txt");
		ofstream of2FlatIdxDebug(ss2flatIdxDebugFileName.str());//("twoFlats_idxPtDebug_compute.txt");
		ofstream ofMajEigDebug("majEigVDebug_compute.txt");
		//ofIdxDebug << " Continuous idx points." << endl;
#endif
		for (size_t i = 0; i < g_xpcpData.size(); i++)
		{
			vector<pFlatPt*> oneFlats;
			g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
			// Add to global list
			g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
			vector<pFlatPt*> twoFlats;
			g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i), g_use_repeat_pcp);
			// Add to global list
			g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
#ifdef IDX_PT_DEBUG 
			// 1-flats
			ofIdxDebug << i << ": ";
			for (size_t jj = 0; jj < oneFlats.size(); jj++)
				ofIdxDebug << oneFlats[jj] << " ";
			ofIdxDebug << endl;
			// 2-flats
			of2FlatIdxDebug << i << ": ";
			for (size_t jj = 0; jj < twoFlats.size(); jj++)
				of2FlatIdxDebug << twoFlats[jj] << " ";
			of2FlatIdxDebug << endl;
			// major eigenvectors
			ofMajEigDebug << i << ":";
			for (int jj = 0; jj < muVol[i].size(); jj++) {
				ofMajEigDebug << muVol[i][jj];
				if (jj == muVol[i].size() - 1)
					ofMajEigDebug << ";";
				else
					ofMajEigDebug << ",";
			}
			for (int jj = 0; jj < g_majEigData[i].size(); jj++)
				ofMajEigDebug << g_majEigData[i][jj] << ",";
			ofMajEigDebug<<endl;

#endif
		}
#ifdef IDX_PT_DEBUG
		ofIdxDebug.close();
		of2FlatIdxDebug.close();
		ofMajEigDebug.close();
#endif
		////////////////////////
		// Set idx points volumes
     	// Get volume data list
		vector<VolumeData*> idxPtVols;
		vector<VolumeData*> idxVolsForWriteOut;
		vector<VolumeData*> idxPt2FlatsVols;
		vector<VolumeData*> idx2FlatsVolsForWriteOut;
		// Eigenvectors
		vector<VolumeData*> majEigVols;
		vector<VolumeData*> secEigVols;

		for (size_t i = 0; i < volList.size(); i++)
		{
			VolumeData* majEigVol = new VolumeData(idxVolDim, 4, NULL, true, true, 1);
			VolumeData* secEigVol = new VolumeData(idxVolDim, 4, NULL, true, true, 1);
			majEigVols.push_back(majEigVol);
			secEigVols.push_back(secEigVol);
			// Only the eigenvector volumes have dimensions of volList.size()-1
			if (i >= volList.size() - 1)
				break;
			// Set VolumeData
			// channels of the index volume: (u,v,strength_flat1_1, strength_flat1_2)
			VolumeData* idxVol = new VolumeData(idxVolDim, 4, NULL, true, true, 4);
			idxPtVols.push_back(idxVol);
			if (i < volList.size() - 2) {
				VolumeData* idx2FlatVol = new VolumeData(idxVolDim, 4, NULL, true, true, 4);
				idxPt2FlatsVols.push_back(idx2FlatVol);
			}
			// Set volumedata for NRRD write out
			VolumeData* idxPerCoordVol1 = new VolumeData(idxVolDim, 4, NULL, true, true);
			VolumeData* idxPerCoordVol2 = new VolumeData(idxVolDim, 4, NULL, true, true);
			idxVolsForWriteOut.push_back(idxPerCoordVol1);
			idxVolsForWriteOut.push_back(idxPerCoordVol2);

			VolumeData* idx2FlatsPerCoordVol1 = new VolumeData(idxVolDim, 4, NULL, true, true);
			VolumeData* idx2FlatsPerCoordVol2 = new VolumeData(idxVolDim, 4, NULL, true, true);
			idx2FlatsVolsForWriteOut.push_back(idx2FlatsPerCoordVol1);
			idx2FlatsVolsForWriteOut.push_back(idx2FlatsPerCoordVol2);

		}
		
		vector<FLOATVECTOR4> oneFlatUVrange(volList.size());
		for (size_t m = 0; m < oneFlatUVrange.size(); m++)
		{
			// range of u
			oneFlatUVrange[m].x = FLT_MAX; 
			oneFlatUVrange[m].y = -FLT_MAX;
			// range of v
			oneFlatUVrange[m].z = FLT_MAX;
			oneFlatUVrange[m].w = -FLT_MAX;
		}
		// Set idx points volumes
		for (size_t ii = 0; ii < samplePos.size(); ii++)
		{
			FLOATVECTOR3 normPos = FLOATVECTOR3((samplePos[ii])[0], (samplePos[ii])[1], (samplePos[ii])[2]);
			//idxPtPosFile << normPos.x << "," << normPos.y << "," << normPos.z << endl;
			UINT64VECTOR3 pos = UINT64VECTOR3(
				normPos.x * float(idxVolDim.x - 1 ) + 0.5f,
				normPos.y * float(idxVolDim.y - 1 ) + 0.5f,
				normPos.z * float(idxVolDim.z - 1) + 0.5f);
			// the indexed points data
			vector<FLOATVECTOR2> oneFlat = g_xpcpData[ii].pcp_1flat;
			vector<FLOATVECTOR2> twoFlat = g_xpcpData[ii].pcp_2flat;

			vector<FLOATVECTOR2> strength1flat = g_xpcpData[ii].strength_1flat;
			vector<FLOATVECTOR2> strength2flat = g_xpcpData[ii].strength_2flat;

			float oneMPer = g_1flat_list[ii]->oneMinusPercentile();
			float twoFlatOneMPer = g_2flat_list[ii]->oneMinusPercentile();

			for (size_t m = 0; m < volList.size(); m++)
			{
				majEigVols[m]->setVoxel(pos, OrgMajEigData[ii][m]);
				secEigVols[m]->setVoxel(pos, OrgSecEigData[ii][m]);

				if (m >= volList.size() - 1)
					break;
				// u
				oneFlatUVrange[m].x = MIN(oneFlatUVrange[m].x, oneFlat[m].x);
				oneFlatUVrange[m].y = MAX(oneFlatUVrange[m].y, oneFlat[m].x);
				// v
				oneFlatUVrange[m].z = MIN(oneFlatUVrange[m].z, oneFlat[m].y);
				oneFlatUVrange[m].w = MAX(oneFlatUVrange[m].w, oneFlat[m].y);

				// Set 1-flats
				// NAN may occur!
				idxPtVols[m]->setVoxel(pos, isnan(oneFlat[m].x) ? 0 : oneFlat[m].x, 0);
				idxPtVols[m]->setVoxel(pos, isnan(oneFlat[m].y) ? 0 : oneFlat[m].y, 1);
				//idxPtVols[m]->setVoxel(pos, strength1flat[m].x, 2);
				// set percentile data
				idxPtVols[m]->setVoxel(pos, oneMPer, 2);
				// set strength data
				idxPtVols[m]->setVoxel(pos, strength1flat[m].y, 3);

				// Set 2-flats
				if (m < volList.size() - 2) {
					idxPt2FlatsVols[m]->setVoxel(pos, isnan(twoFlat[m].x) ? 0 : twoFlat[m].x, 0);
					idxPt2FlatsVols[m]->setVoxel(pos, isnan(twoFlat[m].y) ? 0 : twoFlat[m].y, 1);
					//idxPtVols[m]->setVoxel(pos, strength1flat[m].x, 2);
					// set percentile data
					idxPt2FlatsVols[m]->setVoxel(pos, twoFlatOneMPer, 2);
					// set strength data
					idxPt2FlatsVols[m]->setVoxel(pos, strength2flat[m].y, 3);
				}

				// Set write-out volumes
				// 1-flats
				idxVolsForWriteOut[2 * m]->setVoxel(pos, isnan(oneFlat[m].x) ? 0 : oneFlat[m].x, 0);
				idxVolsForWriteOut[2 * m + 1]->setVoxel(pos, isnan(oneFlat[m].y) ? 0 : oneFlat[m].y, 0);		


				// 2-flats
				if (m >= volList.size() - 2) // skip impossible indices for 2-flats
					continue;
				idx2FlatsVolsForWriteOut[2 * m]->setVoxel(pos, isnan(twoFlat[m].x) ? 0 : twoFlat[m].x, 0);
				idx2FlatsVolsForWriteOut[2 * m + 1]->setVoxel(pos, isnan(twoFlat[m].y) ? 0 : twoFlat[m].y, 0);
			}		 
		}

		for (size_t m = 0; m < oneFlatUVrange.size(); m++)
		{
			cout << "attr " << m << ": "<< oneFlatUVrange[m] << endl;
		}
		// Set global parameters for rendering
		g_params.setIdxPtVols(idxPtVols);
		g_params.setIdxPt2flatsVols(idxPt2FlatsVols);
		g_params.setMajEigVols(majEigVols);
		g_params.setSecEigVols(secEigVols);

		//////////////////////////////////
		// Write out 1-flats as a file??
		if (writeOutFile)
		{

			auto now = std::chrono::system_clock::now();
			auto UTC = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();

			cout << "Start to write out idx points into files" << endl;
			auto in_time_t = std::chrono::system_clock::to_time_t(now);
			std::tm* struct_time = std::localtime(&in_time_t);
			std::stringstream datetime;
			//datetime << std::put_time(std::localtime(&in_time_t), "%Y_%m_%d_%X");
			datetime << struct_time->tm_year+1900 << "_" << struct_time->tm_mon+1 << "_" <<
				struct_time->tm_mday << "_" << struct_time->tm_hour << struct_time->tm_min << struct_time->tm_sec;

			// Time stamp is disabled
			std::string datetimeFixed = "";
			// Uncomment ot enable time stamp0
			// datetimeFixed = datetime.str();
			cout << datetimeFixed << endl;
			//std::ofstream idxPtPosFile("./pos_idx_pt" + datetime.str() + ".txt");
			// Write out volume data
			char* finalVolName = new char[255];

			for (size_t i = 0; i < volList.size(); i++)
			{
				// Write out eigenvector volumes
				string id = number2String(i);
	
				string majEigVolFile = "./majEigVol_" + id + "-xyz-" + number2String(idxVolDim.x) + string("X") + number2String(idxVolDim.y) + string("X") + number2String(idxVolDim.z) + ".nhdr";
				string secEigVolFile = "./secEigVol_" + id + "-xyz-" + number2String(idxVolDim.x) + string("X") + number2String(idxVolDim.y) + string("X") + number2String(idxVolDim.z) + ".nhdr";

				if (!majEigVols[i]->writeToNrrdFile(majEigVolFile.c_str()))
				{
					cout << "Writing file " << majEigVolFile << " failed!" << endl;
				}
				if (!secEigVols[i]->writeToNrrdFile(secEigVolFile.c_str()))
				{
					cout << "Writing file " << secEigVolFile << " failed!" << endl;
				}
				if (i >= idxPtVols.size())
					break;
				// Writeout indexed points volumes
				float minVal = FLT_MAX;
				float maxVal = -FLT_MAX;

				memset(finalVolName, 0, 255 * sizeof(char));

				string idxVolFile = "./vol_sd_idx_" + id + "_" + datetimeFixed;
				string idx2FlatVolFile = "./vol_2flat_idx_" + id + "_" + datetimeFixed;

				if (!idxPtVols[i]->writeToFile(idxVolFile.c_str(), finalVolName))
				{
					cout << "Writing file " << finalVolName << " failed!" << endl;
				}
				// write individual coordinates
				string idxPerVolFile = idxVolFile + "X";
				string idxPerVolFile2 = idxVolFile + "Y";

				string idx2flatPerVolFile = idx2FlatVolFile + "X";
				string idx2flatPerVolFile2 = idx2FlatVolFile + "Y";

				string finalNrrdName;
				// 1-flats
				if (!idxVolsForWriteOut[2 * i]->writeToFile(idxPerVolFile.c_str(), finalVolName))
				{
					cout << "Writing file " << idxPerVolFile << " failed!" << endl;
				}
				else {
					finalNrrdName = SysTools::ChangeExt(finalVolName, "nhdr");
					QString makeNrrdCmd = QString("unu make -s %1 %2 %3 -t float -i %4 -o %5")
						.arg(idxVolDim.x).arg(idxVolDim.y).arg(idxVolDim.z).arg(finalVolName).arg(finalNrrdName.c_str());
					system(makeNrrdCmd.toStdString().c_str());
				}

				if (!idxVolsForWriteOut[2 * i + 1]->writeToFile(idxPerVolFile2.c_str(), finalVolName))
				{
					cout << "Writing file " << idxPerVolFile2 << " failed!" << endl;
				}
				else {
					finalNrrdName = SysTools::ChangeExt(finalVolName, "nhdr");
					QString makeNrrdCmd = QString("unu make -s %1 %2 %3 -t float  -i %4 -o %5")
						.arg(idxVolDim.x).arg(idxVolDim.y).arg(idxVolDim.z).arg(finalVolName).arg(finalNrrdName.c_str());
					system(makeNrrdCmd.toStdString().c_str());
				}
				// 2-flats
				if (i >= idxPt2FlatsVols.size())
					continue;

				if (!idx2FlatsVolsForWriteOut[2 * i]->writeToFile(idx2flatPerVolFile.c_str(), finalVolName))
				{
					cout << "Writing file " << idx2flatPerVolFile << " failed!" << endl;
				}
				else {
					finalNrrdName = SysTools::ChangeExt(finalVolName, "nhdr");
					QString makeNrrdCmd = QString("unu make -s %1 %2 %3 -t float -i %4 -o %5")
						.arg(idxVolDim.x).arg(idxVolDim.y).arg(idxVolDim.z).arg(finalVolName).arg(finalNrrdName.c_str());
					system(makeNrrdCmd.toStdString().c_str());
				}

				if (!idx2FlatsVolsForWriteOut[2 * i + 1]->writeToFile(idx2flatPerVolFile2.c_str(), finalVolName))
				{
					cout << "Writing file " << idx2flatPerVolFile2 << " failed!" << endl;
				}
				else {
					finalNrrdName = SysTools::ChangeExt(finalVolName, "nhdr");
					QString makeNrrdCmd = QString("unu make -s %1 %2 %3 -t float  -i %4 -o %5")
						.arg(idxVolDim.x).arg(idxVolDim.y).arg(idxVolDim.z).arg(finalVolName).arg(finalNrrdName.c_str());
					system(makeNrrdCmd.toStdString().c_str());
				}

				//// Try nrrd write functions
				//// 1-flats
				////finalNrrdName = idxPerVolFile+ ".nhdr";
				//if (!idxVolsForWriteOut[2 * i]->writeToNrrdFile(string(idxPerVolFile+".nhdr").c_str()))
				//{
				//	cout << "Writing file " << finalNrrdName << " failed!" << endl;
				//}
				//
				////finalNrrdName = idxPerVolFile2 +".nhdr";
				//if (!idxVolsForWriteOut[2 * i+1]->writeToNrrdFile(string(idxPerVolFile2 + ".nhdr").c_str()))
				//{
				//	cout << "Writing file " << finalNrrdName << " failed!" << endl;
				//}
				//// 2-flats
				//// 

				////finalNrrdName = idx2flatPerVolFile + ".nhdr";
				//if (!idx2FlatsVolsForWriteOut[2 * i]->writeToNrrdFile(string(idx2flatPerVolFile + ".nhdr").c_str()))
				//{
				//	cout << "Writing file " << finalNrrdName << " failed!" << endl;
				//}

				////finalNrrdName = idx2flatPerVolFile2 + ".nhdr";
				//if (!idx2FlatsVolsForWriteOut[2 * i + 1]->writeToNrrdFile(string(idx2flatPerVolFile2 + ".nhdr").c_str()))
				//{
				//	cout << "Writing file " << finalNrrdName << " failed!" << endl;
				//}

			}
			//// write out sample position
			//std::ofstream idxPtPosFile("./pos_idx_pt.txt");
			//if (!idxPtPosFile.is_open())
			//	cout << "File open error! Cannot write out idx pt pos info!" << endl;
			//else
			//{
			//	idxPtPosFile << "x,y,z" << endl;
			//	for (size_t i = 0; i < samplePos.size(); i++)
			//	{
			//		FLOATVECTOR3 normPos = FLOATVECTOR3((samplePos[i])[0], (samplePos[i])[1], (samplePos[i])[2]);
			//		idxPtPosFile << normPos.x << "," << normPos.y << "," << normPos.z << endl;
			//	}
			//	idxPtPosFile.flush();
			//	idxPtPosFile.close();
			//}

			//std::ofstream idxPtValFile("./val_idx_pt" + datetimeFixed + ".txt");
			////std::ofstream idxPtValFile("./val_idx_pt.txt");
			//for (size_t i = 0; i < g_xpcpData.size(); i++)
			//{
			//	vector<float> pd;
			//	g_xpcpData[i].toStdVector(pd);
			//	for(size_t j = 0; j < pd.size(); j++)
			//		idxPtValFile << pd[j]<<",";
			//	idxPtValFile << endl;
			//}
			//idxPtValFile.flush();
			//idxPtValFile.close();

		}
		// Remove all temporary volumes
		for (size_t i = 0; i < idxVolsForWriteOut.size(); i++)
		{
			SAFE_DELETE(idxVolsForWriteOut[i]);
		}
		idxVolsForWriteOut.clear();

		for (size_t i = 0; i < idx2FlatsVolsForWriteOut.size(); i++)
		{
			SAFE_DELETE(idx2FlatsVolsForWriteOut[i]);
		}
		idx2FlatsVolsForWriteOut.clear();
		// //NOTE: At this point, we have set up the volume data, and we need to reduce points to be drawn for PCP
		// //If we use the grid samples, we need to keep only a small portion of them for PCP
		//UINT64 max_sample_for_pc = 100 * 100 * 20;
		//if (samplePos.size() > max_sample_for_pc)
		//{
		//	float sampleRateForPCP = 0.01f;
		//	size_t skipStepSize = size_t(1.0f/sampleRateForPCP) ;
		//	size_t reducedSize = size_t(sampleRateForPCP * g_xpcpData.size());
		//	vector<XPCPSample> reduceXPCPdata;		
		//	vector<pFlatPt*>   reduce1FlatList;	
		//	vector<pFlatPt*>   reduce2FlatList;
		//	vector<md_val>     reducePtCenterData;

		//	reduceXPCPdata.reserve(reducedSize);
		//	reduce1FlatList.reserve(reducedSize);
		//	reduce2FlatList.reserve(reducedSize);
		//	reducePtCenterData.reserve(reducedSize);

		//	size_t cnt = 0;
		//	for (size_t i = 0; i < g_xpcpData.size(); i+=skipStepSize)
		//	{
		//		reduceXPCPdata.push_back(g_xpcpData[i]);
		//		pFlatPt* oneflat = g_1flat_list[i];
		//		pFlatPt* twoflat = g_2flat_list[i];
		//		oneflat->setRawId(cnt);
		//		twoflat->setRawId(cnt);
		//		reduce1FlatList.push_back(oneflat);
		//		reduce2FlatList.push_back(twoflat);
		//		reducePtCenterData.push_back(pointCenterData[i]);
		//		cnt++;
		//	}

		//	for (size_t i = 0; i < g_xpcpData.size(); i++)
		//	{
		//		if (i % skipStepSize != 0) {
		//			SAFE_DELETE(g_1flat_list[i]);
		//			SAFE_DELETE(g_2flat_list[i]);
		//		}
		//	}
		//	g_xpcpData.clear();
		//	g_1flat_list.clear();
		//	g_2flat_list.clear();
		//	pointCenterData.clear();

		//	g_xpcpData = reduceXPCPdata;
		//	g_1flat_list = reduce1FlatList;
		//	g_2flat_list = reduce2FlatList;
		//	pointCenterData = reducePtCenterData;
		//}

		/////////////////////////////////////////////////////

// Delete raw kd tree?
		if (_annKdTree_rawData)
		{
			SAFE_DELETE(_annKdTree_rawData);
			_annKdTree_rawData = NULL;
		}
		//// Output major eigen vectors
		//ofstream ofMajEig("majEig.txt");
		//for (size_t i = 0; i < g_majEigData.size(); i++)
		//	ofMajEig << g_majEigData[i] << endl;
		//ofMajEig.close();
		////////////////////////////////////////////////
		// TODO: comment out the following until the separators when test is done!
		//// Write result of 1-flats
		//ofstream ofOneFlatTest("oneFlatTest.txt");
		//for (size_t i = 0; i < g_1flat_list.size(); i++)
		//	ofOneFlatTest << *(g_1flat_list[i]) << endl;
		//ofOneFlatTest.close();
		////////////////////////////////////////////////

		// Remove duplicate records

		removeDuplicateRecords(pointCenterData);

		// Clustering
		vector<int> labelData;
		if (g_params.doClustering())
			cluster(g_unique_sampl_raw, labelData, g_params.numClusters());
		/////////////////////////////////////////////
		// TODO: comment out when test is done
		/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
		for (size_t i = 0; i < g_1flat_list.size(); i++)
		ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
		ofUniqueOneFlatTest.close();*/
		////////////////////////////////////////////////

		// Build KD trees for user query
		g_0flat_list.resize(g_unique_sampl_raw.size());
		if (g_params.doClustering())
		{
			for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
				g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i], labelData[i]);
			}
		}
		else
		{
			for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
				g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
			}
		}

		if (g_kdTree_rawData) {
			g_kdTree_rawData->clear();
			g_kdTree_rawData = NULL;
		}
		g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
		g_kdTree_rawData->buildBalanceTree(g_0flat_list);

		// Build KD tree for 1Flat
		if (g_kdTree_1flat) {
			g_kdTree_1flat->clear();
			g_kdTree_1flat = NULL;
		}
		g_kdTree_1flat = new KD<pFlatPt*>(2);
		g_kdTree_1flat->buildBalanceTree(g_1flat_list);
		//////////////////////////////////////////////////
		// TODO: Uncomment to write to file!
		//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
		//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
		//ofOneFlatKdTree.close();
		//////////////////////////////////////////////////
		// Build KD tree for 2Flat
		if (g_kdTree_2flat) {
			g_kdTree_2flat->clear();
			g_kdTree_2flat = NULL;
		}
		g_kdTree_2flat = new KD<pFlatPt*>(2);
		g_kdTree_2flat->buildBalanceTree(g_2flat_list);
		//////////////////////////////////////////////////
		// TODO: Uncomment to write to file!
		//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
		//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
		//ofTwoFlatKdTree.close();
		/////////////////////////////////////////////////

		return true;
}

bool XPCPdataAnalyzer::processDataVol_load(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& neighborInfoVolList, const std::vector<std::vector<float>>& samplePos)
{
	return false;
}

bool XPCPdataAnalyzer::processDataVol_loadIdxVol(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& idxVolData, const std::vector<VolumeData*>& idx2FlatsVolData, const std::vector<std::vector<float>>& samplePos)
{
	//NOTE: we need to make sure that the strength of idx points are not all ZEROS to be drawn correctly!!!
	if (volList.empty())
		return false;

	VolumeData* vol1 = volList[0];
	if (vol1 == NULL)
		return false;

	UINT64VECTOR3 volSize = vol1->getDim();
	// use kdTree_rawData to find nearest neighbors to compute 1-flat &processVecRawData 2-flat

	UINT64 dim_rawData = volList.size();
	//Set dimensions of the data
	_numSamples = UINT64(samplePos.size());
	_dimRawData = dim_rawData; // Raw data dimension
	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
	// Do procssing
	double eps = 0.0; // error bound
	ANNpoint     queryPt; // query point
	ANNidxArray  nnIdx;   // near neighbor indices
	ANNdistArray dists;   // near neighbor distances   

	//Preperations: Clear list
	g_xpcpData.resize(samplePos.size(), XPCPSample(dim_rawData));
	// build pcp if necessary
	if (_pcp == NULL)
	{
		vector<string> attribs(dim_rawData);
		for (size_t i = 0; i < attribs.size(); i++)
			attribs[i] = string("attr_") + number2String(int(i));
		_pcp = new PCPInselberg(attribs);
	}

	//// Create a temporary structure to build the KD-tree
	//vector<vector<float>> vXpcpData(pointCenterData.size());
	// Allocate space for Major eigen data
	//g_majEigData.resize(samplePos.size());
	// use global variables to keep the p-flat records
	g_1flat_list.clear();
	g_2flat_list.clear();

	// Prepare variables to record min/max of p-flat strength in each subspace
	_1flatsMinMaxPerSubspace.clear();
	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
	{
		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	_2flatsMinMaxPerSubspace.clear();
	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
	{
		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	double r = 0.1;       // In a r-ball setting, set the radius respect to the full data range of [0,1].
	ANNdist sqRad = r * r; // Squared radius of the neighborhood. 
	// Compute 1-flat and 2-flat for each and every item in pointCenterData
	ofstream ofDebug("neighborFitTest.txt");
	vector<md_val> pointCenterData(samplePos.size());
	vector<md_val> oneflatData(samplePos.size());
	vector<md_val> twoflatData(samplePos.size());
	int numIdxPtChannels = 2;// idxVolData[0]->numChannels();
	string dummy;
	PointMonteCarloSampler* sampler = new PointMonteCarloSampler(volList, volSize, dummy);
	// TODO: 1. need to sample all voxels of the volume
	// 2. Use samplePos as downsampled 
	for (vector<vector<float>>::const_iterator IT = samplePos.begin(); IT != samplePos.end(); ++IT)
	{
		UINT64 i = IT - samplePos.begin();
		FLOATVECTOR3 normPos = FLOATVECTOR3((*IT)[0]/(volSize.x-1), (*IT)[1]/(volSize.y-1), (*IT)[2] / (volSize.z - 1));
		// sample the center 
		pointCenterData[i] = sampler->sampleVoxelWithInterpNormCoords(normPos);
		// Do not do interpolation for idx pt volumes
		oneflatData[i] = PointMonteCarloSampler::readVoxelNormalCoords(normPos, idxVolData);
		// two flats
		twoflatData[i] = PointMonteCarloSampler::readVoxelNormalCoords(normPos, idx2FlatsVolData);
//		// Compute mean vector
//		Eigen::VectorXf mu = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(pointCenterData[i].data(), pointCenterData[i].size());
//#if OUTPUT_FILE>0
//
//
//		if (i == 0)
//		{
//			cout << "Mean = " << endl << mu << endl;
//			ofDebug << "Mean = " << endl << mu << endl;
//		}
//#endif
//		// The correlation matrix recording correlations between every pair of attributes
//		Eigen::MatrixXf corr;
//		corr_mat = Eigen::MatrixXf::Zero(dim_rawData, dim_rawData);

		XPCPSample xpcp_tuple(dim_rawData);
	

		// compute 1-flat and 2-flat and record min/max of the strength of p-flats
		//calc_p_flats_fromLowFlats(mu, oneflatData[i][0], oneflatData[i][1], oneflatData[i][2], xpcp_tuple);

		// Setup 0- and 1-flats
		xpcp_tuple.x = pointCenterData[i];
		// 1-flat
		md_val flat1 = oneflatData[i];

		xpcp_tuple.pcp_1flat.resize(flat1.size()/numIdxPtChannels);
		for (size_t ii = 0; ii < xpcp_tuple.pcp_1flat.size(); ii++) {
			
			xpcp_tuple.pcp_1flat[ii] = FLOATVECTOR2(flat1[ii * numIdxPtChannels + 0], flat1[ii * numIdxPtChannels + 1]);
			/////////////////////////////////////////////////////////
		//NOTE: 
		// TODO:
		// We need to compute the Major eigenvectors anyway!
		//float eigVLen = Eigen::Vector2f(majEig[i], majEig[i + 1]).norm();
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
			xpcp_tuple.strength_1flat[ii].x = 1.0f; // eigVLen;
			xpcp_tuple.strength_1flat[ii].y = 1.0f;
		}

		// 2-flat
		md_val flat2 = twoflatData[i];
		xpcp_tuple.pcp_2flat.resize(flat2.size() / numIdxPtChannels);
		for (size_t ii = 0; ii < xpcp_tuple.pcp_2flat.size(); ii++) {
			xpcp_tuple.pcp_2flat[ii] = FLOATVECTOR2(flat2[ii * numIdxPtChannels + 0], flat2[ii * numIdxPtChannels + 1]);
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
			xpcp_tuple.strength_2flat[ii].x = 1.0f; // eigVLen;
			xpcp_tuple.strength_2flat[ii].y = 1.0f;
		}
		//Eigen::VectorXf majEigV;
		//majEigV = Eigen::VectorXf::Zero(dim_rawData); 
		//for(size_t k = 0; k < dim_rawData; k++)
		//	majEigV(k) = oneflatData[i][k];
		//// store 2-flat
		//if (has2flat) {
		//	xpcp_tuple.pcp_2flat.resize(flat2.size());
		//	for (size_t i = 0; i < flat2.size(); i++)
		//		xpcp_tuple.pcp_2flat[i] = FLOATVECTOR2(flat2[i].x, flat2[i].y);
		//}

		//g_majEigData[i] = majEigV; // Set Eigen vector data
		g_xpcpData[i] = xpcp_tuple;

		///////////////////////////////////////////////////////////////////////
		// !1 TODO: If we turn off normalization (!2), uncomment codes in the block
		//// Convert p-flat data and store to global variables
		//vector<pFlatPt*> oneFlats;
		//xpcp_tuple.to_pFlatRec(oneFlats, 1, UINT64(i));
		//// Add to global list
		//g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
		//vector<pFlatPt*> twoFlats;
		//xpcp_tuple.to_pFlatRec(twoFlats, 2, UINT64(i));
		//// Add to global list
		//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
		////////////////////////////////////////////////////////////////////////

	}
	ofDebug.close();
	cout << endl << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
	/////////////////////////////////////////////////////
	// !2 Sep 22, 2016: Added normalization func for p-flats 
	// Do normalization after we compute all p-flats!
	//int normMethod = 0;
	//nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
#ifdef IDX_PT_DEBUG 
	ofstream ofIdxDebug("idxPtDebug_load.txt");
	ofIdxDebug << " Continuous idx points." << endl;
#endif

	for (size_t i = 0; i < g_xpcpData.size(); i++)
	{
		vector<pFlatPt*> oneFlats;
		g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
		// Add to global list
		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
		vector<pFlatPt*> twoFlats;
		g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i), g_use_repeat_pcp);
		// Add to global list
		g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
#ifdef IDX_PT_DEBUG
		ofIdxDebug << i << ": ";
		for (size_t jj = 0; jj < oneFlats.size(); jj++)
			ofIdxDebug << oneFlats[jj] << " ";
		ofIdxDebug << endl;
#endif
	}

#ifdef IDX_PT_DEBUG
	ofIdxDebug.close();
#endif
	/////////////////////////////////////////////////////

	// Delete raw kd tree?
	if (_annKdTree_rawData)
	{
		SAFE_DELETE(_annKdTree_rawData);
		_annKdTree_rawData = NULL;
	}
	
	removeDuplicateRecords(pointCenterData);

	// Clustering
	vector<int> labelData;
	if (g_params.doClustering())
		cluster(g_unique_sampl_raw, labelData, g_params.numClusters());
	/////////////////////////////////////////////
	// TODO: comment out when test is done
	/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
	for (size_t i = 0; i < g_1flat_list.size(); i++)
	ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
	ofUniqueOneFlatTest.close();*/
	////////////////////////////////////////////////
	// Set indexed points volumes for Volume Rendering 
	// Need to reorganize the indexed volumes from single channel to multichannel volumes
	vector<VolumeData*> idxPtVols;
	vector<VolumeData*> idxPt2FlatsVols;

	// Set 1-flats
	for (size_t i = 0; i < idxVolData.size() ; i+=2)
	{
		float minY = FLT_MAX;
		float maxY = -FLT_MAX;
		VolumeData* idxVol = new VolumeData(idxVolData[i]->getDim(), 4, NULL, true, true, 4);
		for(size_t zz = 0; zz < idxVol->getDim().z; zz++)
			for(size_t yy = 0; yy < idxVol->getDim().y; yy++)
				for (size_t xx = 0; xx < idxVol->getDim().x; xx++)
				{
					UINT64VECTOR3 pos = UINT64VECTOR3(xx, yy, zz);
					float idxX = idxVolData[i]->getVoxel(pos);
					float idxY = idxVolData[i+1]->getVoxel(pos);
					minY = MIN(idxY, minY);
					maxY = MAX(idxY, maxY);
					idxVol->setVoxel(pos, idxX, 0);
					idxVol->setVoxel(pos, idxY, 1);

					// Channels 2, 3 are for indexed points strength
					idxVol->setVoxel(pos, idxVolData[i]->getVoxel(pos), 2);
					idxVol->setVoxel(pos, idxVolData[i + 1]->getVoxel(pos), 3);

				}
		idxPtVols.push_back(idxVol);
		cout << "min idxY = " << minY << ", max idxY = " << maxY << endl;
	}

	// Set 2-flats
	for (size_t i = 0; i < idx2FlatsVolData.size(); i += 2)
	{
		float minY = FLT_MAX;
		float maxY = -FLT_MAX;
		VolumeData* idxVol = new VolumeData(idx2FlatsVolData[i]->getDim(), 4, NULL, true, true, 4);
		for (size_t zz = 0; zz < idxVol->getDim().z; zz++)
			for (size_t yy = 0; yy < idxVol->getDim().y; yy++)
				for (size_t xx = 0; xx < idxVol->getDim().x; xx++)
				{
					UINT64VECTOR3 pos = UINT64VECTOR3(xx, yy, zz);
					float idxX = idx2FlatsVolData[i]->getVoxel(pos);
					float idxY = idx2FlatsVolData[i + 1]->getVoxel(pos);
					minY = MIN(idxY, minY);
					maxY = MAX(idxY, maxY);
					idxVol->setVoxel(pos, idxX, 0);
					idxVol->setVoxel(pos, idxY, 1);

					// Channels 2, 3 are for indexed points strength
					idxVol->setVoxel(pos, idx2FlatsVolData[i]->getVoxel(pos), 2);
					idxVol->setVoxel(pos, idx2FlatsVolData[i + 1]->getVoxel(pos), 3);

				}
		idxPt2FlatsVols.push_back(idxVol);
		cout << "min idxY = " << minY << ", max idxY = " << maxY << endl;
	}
	// Set indexed points volumes for OpenGL
	g_params.setIdxPtVols(idxPtVols); // 1-flats
	g_params.setIdxPt2flatsVols(idxPt2FlatsVols); // 2-flats
	// Build KD trees for user query
	g_0flat_list.resize(g_unique_sampl_raw.size());
	if (g_params.doClustering())
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i], labelData[i]);
		}
	}
	else
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
		}
	}

	///////////////////////////////
	// disable kd trees construction
  	if (g_kdTree_rawData) {
		g_kdTree_rawData->clear();
		g_kdTree_rawData = NULL;
	}
	g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
	g_kdTree_rawData->buildBalanceTree(g_0flat_list);

	// Build KD tree for 1Flat
	if (g_kdTree_1flat) {
		g_kdTree_1flat->clear();
		g_kdTree_1flat = NULL;
	}
	g_kdTree_1flat = new KD<pFlatPt*>(2);
	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
	//ofOneFlatKdTree.close();
	//////////////////////////////////////////////////
	// Build KD tree for 2Flat
	if (g_kdTree_2flat) {
		g_kdTree_2flat->clear();
		g_kdTree_2flat = NULL;
	}
	g_kdTree_2flat = new KD<pFlatPt*>(2);
	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
	//ofTwoFlatKdTree.close();
	/////////////////////////////////////////////////

	return true;
}

bool XPCPdataAnalyzer::processDataVol_loadEigVecVol(const std::vector<VolumeData*>& volList, const std::vector<VolumeData*>& majEigVols, const std::vector<VolumeData*>& secEigVols, const std::vector<std::vector<float>>& samplePos)
{
	//NOTE: we need to make sure that the strength of idx points are not all ZEROS to be drawn correctly!!!
	if (volList.empty())
		return false;

	if (majEigVols.empty())
	{
		cout << "No major eigenvector volumes provided! Force quit!" << endl;
		return false;
	}

	bool compute2Flats = true;
	if (secEigVols.empty()) {
		compute2Flats = false;
		cout << "No second major eigenvector volumes provided! 2-flats will not be comptued!" << endl;
	}
	VolumeData* vol1 = volList[0];
	if (vol1 == NULL)
		return false;
	if (majEigVols[0] == NULL)
		return false;
	// The size of the indexed volume 
	UINT64VECTOR3 idxVolDim = majEigVols[0]->getDim();
	// The size of the data volume 
	UINT64VECTOR3 volSize = volList[0]->getDim();
	// It is possible that volSize != idxVolDim

	// use kdTree_rawData to find nearest neighbors to compute 1-flat &processVecRawData 2-flat

	UINT64 dim_rawData = volList.size();


	//Set dimensions of the data
	_numSamples = UINT64(samplePos.size());
	_dimRawData = dim_rawData; // Raw data dimension
	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
	// Do procssing
	double eps = 0.0; // error bound
	ANNpoint     queryPt; // query point
	ANNidxArray  nnIdx;   // near neighbor indices
	ANNdistArray dists;   // near neighbor distances   

	//Preperations: Clear list
	g_xpcpData.resize(samplePos.size(), XPCPSample(dim_rawData));
	// build pcp if necessary
	if (_pcp == NULL)
	{
		vector<string> attribs(dim_rawData);
		for (size_t i = 0; i < attribs.size(); i++)
			attribs[i] = string("attr_") + number2String(int(i));
		_pcp = new PCPInselberg(attribs);
	}

	//// Create a temporary structure to build the KD-tree
	//vector<vector<float>> vXpcpData(pointCenterData.size());
	// Allocate space for Major eigen data
	g_majEigData.resize(samplePos.size());
	// use global variables to keep the p-flat records
	g_1flat_list.clear();
	g_2flat_list.clear();

	// Prepare variables to record min/max of p-flat strength in each subspace
	_1flatsMinMaxPerSubspace.clear();
	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
	{
		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	_2flatsMinMaxPerSubspace.clear();
	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
	{
		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	// Compute 1-flat and 2-flat for each and every item in pointCenterData
	vector<md_val> pointCenterData(samplePos.size());
	vector<md_val> majEigData(samplePos.size());
	vector<md_val> secEigData(samplePos.size());
	int numIdxPtChannels = 2;// idxVolData[0]->numChannels();
	string dummy;
	PointMonteCarloSampler* sampler = new PointMonteCarloSampler(volList, volSize, dummy);
	// TODO: 1. need to sample all voxels of the volume
	// 2. Use samplePos as downsampled 



	// vectors for indexed points computation
	Eigen::VectorXf mu(_dimRawData);
	Eigen::VectorXf majEigV(_dimRawData);
	Eigen::VectorXf secEigV(_dimRawData);
	if (!compute2Flats)
		secEigV.resize(0);
	
	for (vector<vector<float>>::const_iterator IT = samplePos.begin(); IT != samplePos.end(); ++IT)
	{
		UINT64 i = IT - samplePos.begin();
		UINT64VECTOR3 iPos = UINT64VECTOR3((*IT)[0], (*IT)[1], (*IT)[2]);
		FLOATVECTOR3 normPos = FLOATVECTOR3((*IT)[0] / (volSize.x - 1), (*IT)[1] / (volSize.y - 1), (*IT)[2] / (volSize.z - 1));
		// sample the center 
		pointCenterData[i] = sampler->sampleVoxelWithInterpNormCoords(normPos);
		// Do not do interpolation for eigenvector volumes
		majEigData[i] = PointMonteCarloSampler::readVoxelNormalCoords(normPos, majEigVols);
		if(compute2Flats)
			secEigData[i] = PointMonteCarloSampler::readVoxelNormalCoords(normPos, secEigVols) ;

		// set VectorXf with the convention of Eigen
		for (size_t k = 0; k < _dimRawData; k++) {
			majEigV(k) = majEigData[i][k];
			if(compute2Flats)
				secEigV(k) = secEigData[i][k];
			mu(k) = pointCenterData[i][k];
		}
		majEigV.normalize();
		if(compute2Flats)
			secEigV.normalize();
		// PCP information
		XPCPSample xpcp_tuple(_dimRawData);

		// compute 1-flat and 2-flat and record min/max of the strength of p-flats
		//if (g_params.PR == PR_SUBSPACE)
		//	calc_p_flats_subSpaceMethod(mu, majEigV, secEigV, xpcp_tuple);
		//else if (method == PR_FROM_LOWER_DIM)
		//	calc_p_flats_fromLowerDimMethod(mu, majEigV, secEigV, xpcp_tuple);
		//else
			calc_p_flats_subSpaceMethod(mu, majEigV, secEigV, xpcp_tuple);

		g_majEigData[i] = majEigV; // Set Eigen vector data

		g_xpcpData[i] = xpcp_tuple;


		// Set indexed volumes as multichannel volumes
		//for(size_t k = 0; k < idxPtVols.size(); k++)
		//{
		//	// Channels 0, 1 are for the indexed point position
		//	UINT64VECTOR3 pos = UINT64VECTOR3(xx, yy, zz);
		//	float idxX = idxVolData[i]->getVoxel(pos);
		//	float idxY = idxVolData[i + 1]->getVoxel(pos);
		//	minY = MIN(idxY, minY);
		//	maxY = MAX(idxY, maxY);
		//	idxVol->setVoxel(pos, idxX, 0);
		//	idxVol->setVoxel(pos, idxY, 1);

		//	// Channels 2, 3 are for indexed points strength
		//	idxVol->setVoxel(pos, idxVolData[i]->getVoxel(pos), 2);
		//	idxVol->setVoxel(pos, idxVolData[i + 1]->getVoxel(pos), 3);
		//}
	}

	cout << endl << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
	/////////////////////////////////////////////////////
	// !2 Sep 22, 2016: Added normalization func for p-flats 
	// Do normalization after we compute all p-flats!
	//int normMethod = 0;
	//nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
#ifdef IDX_PT_DEBUG 
	ofstream ofIdxDebug("idxPtDebug_load.txt");
	ofIdxDebug << " Continuous idx points before KDE." << endl;
#endif

	for (size_t i = 0; i < g_xpcpData.size(); i++)
	{
		vector<pFlatPt*> oneFlats;
		g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
		// Add to global list
		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
		vector<pFlatPt*> twoFlats;
		g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i), g_use_repeat_pcp);
		// Add to global list
		g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
#ifdef IDX_PT_DEBUG
		ofIdxDebug << i << ": ";
		for (size_t jj = 0; jj < oneFlats.size(); jj++)
			ofIdxDebug << oneFlats[jj] << " ";
		ofIdxDebug << endl;
#endif
	}

#ifdef IDX_PT_DEBUG
	ofIdxDebug.close();
#endif
	// Set indexed points volumes
// Need to reorganize the indexed volumes from single channel to multichannel volumes
	vector<VolumeData*> idxPtVols;

	for (size_t i = 0; i < volList.size(); i ++)
	{
		if (i >= volList.size() - 1)
			break;
		VolumeData* idxVol = new VolumeData(idxVolDim, 4, NULL, true, true, 4);
		idxPtVols.push_back(idxVol);
	}
	// Set idx points volumes
	for (size_t ii = 0; ii < samplePos.size(); ii++)
	{
		FLOATVECTOR3 normPos = FLOATVECTOR3((samplePos[ii])[0]/(volSize.x-1), 
			(samplePos[ii])[1] / (volSize.y - 1), 
			(samplePos[ii])[2] / (volSize.z - 1));
		//idxPtPosFile << normPos.x << "," << normPos.y << "," << normPos.z << endl;
		UINT64VECTOR3 pos = UINT64VECTOR3(
			normPos.x * float(idxVolDim.x - 1) + 0.5f,
			normPos.y * float(idxVolDim.y - 1) + 0.5f,
			normPos.z * float(idxVolDim.z - 1) + 0.5f);
		// the indexed points data
		vector<FLOATVECTOR2> oneFlat = g_xpcpData[ii].pcp_1flat;
		vector<FLOATVECTOR2> twoFlat = g_xpcpData[ii].pcp_2flat;

		vector<FLOATVECTOR2> strength1flat = g_xpcpData[ii].strength_1flat;
		vector<FLOATVECTOR2> strength2flat = g_xpcpData[ii].strength_2flat;

		float oneMPer = g_1flat_list[ii]->oneMinusPercentile();
		float twoFlatOneMPer = g_2flat_list[ii]->oneMinusPercentile();

		for (size_t m = 0; m < volList.size(); m++)
		{

			if (m >= volList.size() - 1)
				break;

			// Set 1-flats
			// NAN may occur!
			// Channels 0, 1 for the 2D location 
			idxPtVols[m]->setVoxel(pos, isnan(oneFlat[m].x) ? 0 : oneFlat[m].x, 0);
			idxPtVols[m]->setVoxel(pos, isnan(oneFlat[m].y) ? 0 : oneFlat[m].y, 1);
			//idxPtVols[m]->setVoxel(pos, strength1flat[m].x, 2);
			// Channel 2: set percentile data
			idxPtVols[m]->setVoxel(pos, oneMPer, 2);
			// Channel 3: set strength data
			idxPtVols[m]->setVoxel(pos, strength1flat[m].y, 3);
		}
	}
	// set the global variable
	g_params.setIdxPtVols(idxPtVols);
	/////////////////////////////////////////////////////

	// Delete raw kd tree?
	if (_annKdTree_rawData)
	{
		SAFE_DELETE(_annKdTree_rawData);
		_annKdTree_rawData = NULL;
	}

	removeDuplicateRecords(pointCenterData);

	// Clustering
	vector<int> labelData;
	if (g_params.doClustering())
		cluster(g_unique_sampl_raw, labelData, g_params.numClusters());
	/////////////////////////////////////////////
	// TODO: comment out when test is done
	/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
	for (size_t i = 0; i < g_1flat_list.size(); i++)
	ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
	ofUniqueOneFlatTest.close();*/
	////////////////////////////////////////////////
	
	// Build KD trees for user query
	g_0flat_list.resize(g_unique_sampl_raw.size());
	if (g_params.doClustering())
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i], labelData[i]);
		}
	}
	else
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
		}
	}

	///////////////////////////////
	// disable kd trees construction
	if (g_kdTree_rawData) {
		g_kdTree_rawData->clear();
		g_kdTree_rawData = NULL;
	}
	//g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
	//g_kdTree_rawData->buildBalanceTree(g_0flat_list);

	// Build KD tree for 1Flat
	if (g_kdTree_1flat) {
		g_kdTree_1flat->clear();
		g_kdTree_1flat = NULL;
	}
	//g_kdTree_1flat = new KD<pFlatPt*>(2);
	//g_kdTree_1flat->buildBalanceTree(g_1flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
	//ofOneFlatKdTree.close();
	//////////////////////////////////////////////////
	// Build KD tree for 2Flat
	if (g_kdTree_2flat) {
		g_kdTree_2flat->clear();
		g_kdTree_2flat = NULL;
	}
	//g_kdTree_2flat = new KD<pFlatPt*>(2);
	//g_kdTree_2flat->buildBalanceTree(g_2flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
	//ofTwoFlatKdTree.close();
	/////////////////////////////////////////////////

	return true;
}

void XPCPdataAnalyzer::buildKDTreeInSubSpaces(const vector<vector<float>>& pointCenterData, int subspaceId, int subspaceDim, ANNkd_tree** pKdTree)
{
	SAFE_DELETE(*pKdTree);
	if (subspaceId + subspaceDim > pointCenterData[0].size())
	{
		cout << "Subspace does not exist." << endl;
		return;
	}
	ANNpointArray dataPts;
	UINT64 nPts = pointCenterData.size();
	UINT64 dim = subspaceDim;
	// Setup ANN points
	dataPts = annAllocPts(nPts, dim);
	for (size_t i = 0; i < pointCenterData.size(); i++)
	{
		for (int j = subspaceId; j < subspaceId + subspaceDim; j++)
			dataPts[i][j - subspaceId] = double(pointCenterData[i][j]);
	}
	// build kd tree
	*pKdTree = new ANNkd_tree(dataPts,
		nPts,
		dim);

}

void XPCPdataAnalyzer::processRawData(const vector<vector<float>>& pointCenterData, int num_nearest_neighbors, PR_METHOD method)
{

	if (pointCenterData.empty())
	{
		cout << "Empty raw data!" << endl;
		return;
	}
	if (_annKdTree_rawData == NULL)
		buildKDTree(pointCenterData, &_annKdTree_rawData);
	// use kdTree_rawData to find nearest neighbors to compute 1-flat &processVecRawData 2-flat
	int k_0 = num_nearest_neighbors;
	int dim_rawData = pointCenterData[0].size();
	//Set dimensions of the data
	_numSamples = UINT64(pointCenterData.size());
	_dimRawData = dim_rawData; // Raw data dimension
	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
	// Do procssing
	double eps = 0.0; // error bound
	ANNpoint     queryPt; // query point
	ANNidxArray  nnIdx;   // near neighbor indices
	ANNdistArray dists;   // near neighbor distances   
	queryPt = annAllocPt(dim_rawData);
	nnIdx = new ANNidx[k_0];
	dists = new ANNdist[k_0];
	//Preperations: Clear list
	g_xpcpData.resize(pointCenterData.size(), XPCPSample(dim_rawData));
	// build pcp if necessary
	if (_pcp == NULL)
	{
		vector<string> attribs(dim_rawData);
		for (size_t i = 0; i < attribs.size(); i++)
			attribs[i] = string("attr_") + number2String(int(i));
		_pcp = new PCPInselberg(attribs);
	}

	//// Create a temporary structure to build the KD-tree
	//vector<vector<float>> vXpcpData(pointCenterData.size());
	// Allocate space for Major eigen data
	g_majEigData.resize(pointCenterData.size());
	// use global variables to keep the p-flat records
	g_1flat_list.clear();
	g_2flat_list.clear();

	// Prepare variables to record min/max of p-flat strength in each subspace
	_1flatsMinMaxPerSubspace.clear();
	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
	{
		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	_2flatsMinMaxPerSubspace.clear();
	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
	{
		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
	}
	double r = 0.1;       // In a r-ball setting, set the radius respect to the full data range of [0,1].
	ANNdist sqRad = r * r; // Squared radius of the neighborhood. 
	// Compute 1-flat and 2-flat for each and every item in pointCenterData
	ofstream ofDebug("neighborFitTest.txt");
	for (UINT64 i = 0; i < UINT64(pointCenterData.size()); i++)
	{
		for (int d = 0; d < dim_rawData; d++)
			queryPt[d] = pointCenterData[i][d];

		// Search for nearest neighbors
		int k = k_0;
		switch (g_params.NnQueryMethod())
		{
		case NNQ_KNN:
			_annKdTree_rawData->annkSearch(// search
				queryPt,	// query point
				k,			// number of near neighbors
				nnIdx,		// nearest neighbors (returned)
				dists,		// distance (returned)
				eps			// error bound
				);
			break;
		case NNQ_RBALL:
			// Search for neighbors inside a hyper-sphere
			k = _annKdTree_rawData->annkFRSearch(// search
				queryPt, // query point
				sqRad, //squared radius
				k_0, // number of near neighbors
				nnIdx, // nearest neighbors (returned)
				dists, // distance (returned)
				eps);// error bound
			// The actual sample from the search may be greater than our k0
			if (k > k_0)
			{
				// In that case, query again to get all neighbors
				SAFE_DELETE_ARRAY(nnIdx); // Need to destroy current arrays and reallocate with size k
				SAFE_DELETE_ARRAY(dists);
				nnIdx = new ANNidx[k];
				dists = new ANNdist[k];
				k = _annKdTree_rawData->annkFRSearch(// search
					queryPt, // query point
					sqRad, //squared radius
					k, // number of near neighbors
					nnIdx, // nearest neighbors (returned)
					dists, // distance (returned)
					eps);// error bound
			}
			//k = k_0; // In that case, we take only k0 samples from the query result
			break;
		default: // The default is KNN
			_annKdTree_rawData->annkSearch(// search
				queryPt,	// query point
				k,			// number of near neighbors
				nnIdx,		// nearest neighbors (returned)
				dists,		// distance (returned)
				eps			// error bound
				);
			break;
		}



		// Get all points for computation

		vector<vector<float>> neighbors(k);
		//---------------------------------------------------------------
		// From ANN manual page 9:
		//https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf
		// the neighbors are sorted by their distances to the query point: i.e., nnIdx[0] is the neareast one, in our case, the query point itself!
		for (int nn = 0; nn < k; nn++)
		{
			int id = nnIdx[nn];
			neighbors[nn] = pointCenterData[id];
		}
#if OUTPUT_FILE>0
		if (i == 0)
		{
			ofstream ofNeighbor0("neighborTest0.txt");
			for (int i = 0; i < k; i++)
			{
				for (size_t d = 0; d < dim_rawData - 1; d++)
					ofNeighbor0 << neighbors[i][d] << ",";
				ofNeighbor0 << neighbors[i][dim_rawData - 1] << endl;
			}
			ofNeighbor0.close();
		}
#endif
		// Compute mean vector
		Eigen::VectorXf mu;
		MyStatistics::ComputeMean<float>(neighbors, mu);
#if OUTPUT_FILE>0
		if (i == 0)
		{
			cout << "Mean = " << endl << mu << endl;
			ofDebug << "Mean = " << endl << mu << endl;
		}
#endif
		// The correlation matrix recording correlations between every pair of attributes
		Eigen::MatrixXf corr;
		corr_mat = Eigen::MatrixXf::Zero(dim_rawData, dim_rawData);

		//// Compute the covariance matrix and then the SVD of that cov matrix
		//Eigen::MatrixXf svdU, svdV;
		//Eigen::VectorXf svdS;
		//MyStatistics::SVDCovMatrix(neighbors, svdU, svdS, svdV);
		//Eigen::VectorXf majSingV = svdU.col(0) * svdS(0);
		//Eigen::VectorXf secMajSingV = svdU.col(1) * svdS(1);
		//if (i == 0)
		//{
		//	cout << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
		//	ofDebug << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
		//}
		//Eigen::VectorXf majEigV = majSingV;
		//Eigen::VectorXf secMajEigV = secMajSingV;

		XPCPSample xpcp_tuple(dim_rawData);
		Eigen::VectorXf majEigV;
		Eigen::VectorXf secMajEigV;
		// For the case of there is no neighbor
		if (k == 1)
		{
			// Set only the x-component of the tuple, set 0 for all other components in the tuple!
			xpcp_tuple.x = neighbors[0];
			majEigV = Eigen::VectorXf::Zero(dim_rawData); // Simple space holder
			secMajEigV = Eigen::VectorXf::Zero(dim_rawData);
		}
		else
		{
			Eigen::VectorXf eigVal;
			Eigen::MatrixXf eigVec;
			//MyStatistics::EigenSolvCovMatrix(neighbors, eigVal, eigVec);
			MyStatistics::EigenSolvCovMatrixGetCorr(neighbors, eigVal, eigVec, corr);
#if OUTPUT_FILE > 0
			if (i == 0)
			{
				cout << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
				ofDebug << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
				// output correlations
				cout << "Correlations: " << endl;
				for (int dy = 0; dy < corr.rows() - 1; dy++)
				{
					int dx = dy + 1;
					cout << "(" << dy << "," << dx << ") = ";
					cout << corr(dy, dx) << endl;
				}
			}
#endif

			majEigV = eigVec.col(eigVec.cols() - 1);
			majEigV.normalize();
			secMajEigV = eigVec.col(eigVec.cols() - 2);
			secMajEigV.normalize();
#if OUTPUT_FILE > 0
			if (i == 0)
			{
				cout << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
				ofDebug << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
			}
#endif
			corr_mat = corr;
		}


		// compute 1-flat and 2-flat and record min/max of the strength of p-flats
		if (method == PR_SUBSPACE)
			calc_p_flats_subSpaceMethod(mu, majEigV, secMajEigV, xpcp_tuple);
		else if (method == PR_FROM_LOWER_DIM)
			calc_p_flats_fromLowerDimMethod(mu, majEigV, secMajEigV, xpcp_tuple);
		else
			calc_p_flats_subSpaceMethod(mu, majEigV, secMajEigV, xpcp_tuple);

		g_majEigData[i] = majEigV; // Set Eigen vector data
		g_xpcpData[i] = xpcp_tuple;

		///////////////////////////////////////////////////////////////////////
		// !1 TODO: If we turn off normalization (!2), uncomment codes in the block
		//// Convert p-flat data and store to global variables
		//vector<pFlatPt*> oneFlats;
		//xpcp_tuple.to_pFlatRec(oneFlats, 1, UINT64(i));
		//// Add to global list
		//g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
		//vector<pFlatPt*> twoFlats;
		//xpcp_tuple.to_pFlatRec(twoFlats, 2, UINT64(i));
		//// Add to global list
		//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
		////////////////////////////////////////////////////////////////////////

	}
	ofDebug.close();
	cout << endl << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
	/////////////////////////////////////////////////////
	// !2 Sep 22, 2016: Added normalization func for p-flats 
	// Do normalization after we compute all p-flats!
	int normMethod = 0;
	nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
	for (size_t i = 0; i < g_xpcpData.size(); i++)
	{
		vector<pFlatPt*> oneFlats;
		g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
		// Add to global list
		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
		vector<pFlatPt*> twoFlats;
		g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i), g_use_repeat_pcp);
		// Add to global list
		g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
	}
	/////////////////////////////////////////////////////

	// Delete raw kd tree?
	if (_annKdTree_rawData)
	{
		SAFE_DELETE(_annKdTree_rawData);
		_annKdTree_rawData = NULL;
	}
	//// Output major eigen vectors
	//ofstream ofMajEig("majEig.txt");
	//for (size_t i = 0; i < g_majEigData.size(); i++)
	//	ofMajEig << g_majEigData[i] << endl;
	//ofMajEig.close();
	////////////////////////////////////////////////
	// TODO: comment out the following until the separators when test is done!
	//// Write result of 1-flats
	//ofstream ofOneFlatTest("oneFlatTest.txt");
	//for (size_t i = 0; i < g_1flat_list.size(); i++)
	//	ofOneFlatTest << *(g_1flat_list[i]) << endl;
	//ofOneFlatTest.close();
	////////////////////////////////////////////////

	// Remove duplicate records

	removeDuplicateRecords(pointCenterData);

	// Clustering
	vector<int> labelData;
	if (g_params.doClustering())
		cluster(g_unique_sampl_raw, labelData, g_params.numClusters());
	/////////////////////////////////////////////
	// TODO: comment out when test is done
	/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
	for (size_t i = 0; i < g_1flat_list.size(); i++)
	ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
	ofUniqueOneFlatTest.close();*/
	////////////////////////////////////////////////

	// Build KD trees for user query
	g_0flat_list.resize(g_unique_sampl_raw.size());
	if (g_params.doClustering())
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i], labelData[i]);
		}
	}
	else
	{
		for (size_t i = 0; i < g_unique_sampl_raw.size(); i++) {
			g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
		}
	}

	if (g_kdTree_rawData){
		g_kdTree_rawData->clear();
		g_kdTree_rawData = NULL;
	}
	g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
	g_kdTree_rawData->buildBalanceTree(g_0flat_list);

	// Build KD tree for 1Flat
	if (g_kdTree_1flat){
		g_kdTree_1flat->clear();
		g_kdTree_1flat = NULL;
	}
	g_kdTree_1flat = new KD<pFlatPt*>(2);
	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
	//ofOneFlatKdTree.close();
	//////////////////////////////////////////////////
	// Build KD tree for 2Flat
	if (g_kdTree_2flat){
		g_kdTree_2flat->clear();
		g_kdTree_2flat = NULL;
	}
	g_kdTree_2flat = new KD<pFlatPt*>(2);
	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
	//ofTwoFlatKdTree.close();
	/////////////////////////////////////////////////

	SAFE_DELETE_ARRAY(nnIdx);
	SAFE_DELETE_ARRAY(dists);
}

//void XPCPdataAnalyzer::processUncertainRawData(const vector<GaussianXd>& distrPointData, int num_nearest_neighbors)
//{
//	// Create a certain point-based raw data by taking the center of distrPointData
//	vector<vector<float>> pointCenterData;
//
//	uncertain_getPointCenters(distrPointData, pointCenterData);
//	// Copy to g_sampledRaw
//	g_sampled_raw = pointCenterData;
//
//	if (_annKdTree_rawData == NULL)
//		buildKDTree(pointCenterData, &_annKdTree_rawData);
//	// use kdTree_rawData to find nearest neighbors to compute 1-flat & 2-flat
//	int k_0 = num_nearest_neighbors;
//	int dim_rawData = pointCenterData[0].size();
//	//Set dimensions of the data
//	_numSamples = UINT64(pointCenterData.size());
//	_dimRawData = dim_rawData; // Raw data dimension
//	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
//	// Do procssing
//	double eps = 0.0; // error bound
//	ANNpoint     queryPt; // query point
//	ANNidxArray  nnIdx;   // near neighbor indices
//	ANNdistArray dists;   // near neighbor distances   
//	queryPt = annAllocPt(dim_rawData);
//	nnIdx = new ANNidx[k_0];
//	dists = new ANNdist[k_0];
//	//Preperations: Clear list
//	g_xpcpData.resize(pointCenterData.size(), UncertainXPCPSample(dim_rawData, g_params.Uncertain_samples_per_distr()));
//	// build pcp if necessary
//	if (_pcp == NULL)
//	{
//		vector<string> attribs(dim_rawData);
//		for (size_t i = 0; i < attribs.size(); i++)
//			attribs[i] = string("attr_") + number2String(int(i));
//		_pcp = new PCPInselberg(attribs);
//	}
//
//	//// Create a temporary structure to build the KD-tree
//	//vector<vector<float>> vXpcpData(pointCenterData.size());
//	// Allocate space for Major eigen data
//	g_majEigData.resize(pointCenterData.size());
//	// use global variables to keep the p-flat records
//	g_1flat_list.clear();
//	g_2flat_list.clear();
//
//	// Prepare variables to record min/max of p-flat strength in each subspace
//	_1flatsMinMaxPerSubspace.clear();
//	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
//	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
//	{
//		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//	_2flatsMinMaxPerSubspace.clear();
//	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
//	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
//	{
//		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//	double r = 0.1;       // In a r-ball setting, set the radius respect to the full data range of [0,1].
//	ANNdist sqRad = r * r; // Squared radius of the neighborhood. 
//	// Compute 1-flat and 2-flat for each and every item in pointCenterData
//	ofstream ofDebug("neighborFitTest.txt");
//	ofstream ofActualNeighbors("samples_for_eigvec.txt");
//	for (UINT64 i = 0; i < UINT64(pointCenterData.size()); i++)
//	{
//		for (int d = 0; d < dim_rawData; d++)
//			queryPt[d] = pointCenterData[i][d];
//
//		// Search for nearest neighbors
//		int k = k_0;
//		switch (g_params.NnQueryMethod())
//		{
//		case NNQ_KNN:
//			_annKdTree_rawData->annkSearch(// search
//				queryPt,	// query point
//				k,			// number of near neighbors
//				nnIdx,		// nearest neighbors (returned)
//				dists,		// distance (returned)
//				eps			// error bound
//				);
//			break;
//		case NNQ_RBALL:
//			// Search for neighbors inside a hyper-sphere
//			k = _annKdTree_rawData->annkFRSearch(// search
//				queryPt, // query point
//				sqRad, //squared radius
//				k_0, // number of near neighbors
//				nnIdx, // nearest neighbors (returned)
//				dists, // distance (returned)
//				eps);// error bound
//			// The actual sample from the search may be greater than our k0
//			if (k > k_0)
//			{
//				// In that case, query again to get all neighbors
//				SAFE_DELETE_ARRAY(nnIdx); // Need to destroy current arrays and reallocate with size k
//				SAFE_DELETE_ARRAY(dists);
//				nnIdx = new ANNidx[k];
//				dists = new ANNdist[k];
//				k = _annKdTree_rawData->annkFRSearch(// search
//					queryPt, // query point
//					sqRad, //squared radius
//					k, // number of near neighbors
//					nnIdx, // nearest neighbors (returned)
//					dists, // distance (returned)
//					eps);// error bound
//			}
//			//k = k_0; // In that case, we take only k0 samples from the query result
//			break;
//		default: // The default is KNN
//			_annKdTree_rawData->annkSearch(// search
//				queryPt,	// query point
//				k,			// number of near neighbors
//				nnIdx,		// nearest neighbors (returned)
//				dists,		// distance (returned)
//				eps			// error bound
//				);
//			break;
//		}
//
//
//
//		// Get all points for computation
//		// For uncertain datasets, we need several different sampled neighbors
//
//		vector<Eigen::MatrixXd> sampled_neighbors(k); // Each entry contains a matrix that records the m-dimensional samples of the neighbor
//		// Monte-carlo sampling of the neighborhood 
//		uncertain_sampleNeighborhood(distrPointData, nnIdx, k, sampled_neighbors, g_params.Uncertain_samples_per_distr());
//#if OUTPUT_FILE > 0
//		// Print out the result for the first neighbor
//		if (i == 0)
//		{
//			ofstream ofNeighbor0("uncertain_neighborTest0.txt");
//			for (int j = 0; j < sampled_neighbors.size(); j++)
//			{
//				ofNeighbor0 << "neighbor " << j << endl;
//				ofNeighbor0 << sampled_neighbors[j] << endl;
//				//for (int c = 0; c < sampled_neighbors[j].cols(); c++)
//				//{
//				//	ofNeighbor0 << sampled_neighbors[j].col(c).transpose() << endl;
//				//}
//				ofNeighbor0 << endl;
//			}
//			ofNeighbor0.close();
//		}
//#endif
//		// For each sampling, we need to fit a 1-flat!
//		Eigen::VectorXd majEigV;
//		Eigen::VectorXd secMajEigV;
//		UncertainXPCPSample xpcp_tuple(dim_rawData, g_params.Uncertain_samples_per_distr());
//
//
//
//		for (int s = 0; s < g_params.Uncertain_samples_per_distr(); s++)
//		{
//			vector<Eigen::VectorXd> neighbors_one_sample;
//			// Get one setting of the neighborhood sampling
//			//uncertain_get_neighbors_one_sample(sampled_neighbors, s, neighbors_one_sample);
//			uncertain_get_neighbors_from_all_samples(sampled_neighbors, neighbors_one_sample);
//#if OUTPUT_FILE > 0
//			if (i == 0)
//			{
//				ofActualNeighbors << "Setting " << s << endl;
//				for (size_t jj = 0; jj < neighbors_one_sample.size(); jj++)
//					ofActualNeighbors << neighbors_one_sample[jj].transpose() << endl;
//			}
//#endif
//			Eigen::VectorXd mu(dim_rawData);
//
//			// For the case that there is no neighbor
//			if (k == 1)
//			{
//				//// Set only the x-component of the tuple, set 0 for all other components in the tuple!
//				//xpcp_tuple.x = vector<float>(neighbors_one_sample[0].data(), neighbors_one_sample[0].data() + neighbors_one_sample[0].size());
//				majEigV = Eigen::VectorXd::Zero(dim_rawData); // Simple space holder
//				secMajEigV = Eigen::VectorXd::Zero(dim_rawData);
//			}
//			else
//			{
//				Eigen::VectorXd eigVal;
//				Eigen::MatrixXd eigVec;
//				//MyStatistics::EigenSolvCovMatrix(neighbors, eigVal, eigVec);
//				MyStatistics::EigenSolvCovMatrix(neighbors_one_sample, mu, eigVal, eigVec);
//				//if (i == 0)
//				//{
//				//	ofDebug << "mu = " << endl << mu << endl;
//				//	cout << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
//				//	ofDebug << "EigVals: " << endl << eigVal << endl << "EigVecs: " << endl << eigVec << endl;
//				//}
//
//
//				majEigV = eigVec.col(eigVec.cols() - 1);
//				majEigV.normalize();
//				secMajEigV = eigVec.col(eigVec.cols() - 2);
//				secMajEigV.normalize();
//#if OUTPUT_FILE > 0
//				if (i == 0)
//				{
//					ofDebug << "mu = " << endl << mu << endl;
//					cout << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
//					ofDebug << "Maj Eig Vec: " << endl << majEigV << endl << "Sec Maj Eig Vec: " << endl << secMajEigV << endl;
//				}
//#endif
//			}
//			if (s == 0)
//				g_majEigData[i] = majEigV.cast<float>(); // Set Eigen vector data for only one sample
//
//			// compute 1-flat and 2-flat and record min/max of the strength of p-flats
//			uncertain_calc_p_flats_subSpaceMethod(mu, majEigV, secMajEigV, xpcp_tuple, s);
//		}
//
//		// Set uncertain xpcp p-flats
//		g_xpcpData[i] = xpcp_tuple;
//
//	}
//	ofActualNeighbors.close();
//	ofDebug.close();
//	cout << endl << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
//	/////////////////////////////////////////////////////
//	// !2 Sep 22, 2016: Added normalization func for p-flats 
//	// Do normalization after we compute all p-flats!
//	int normMethod = 0;
//	nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
//	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
//	for (size_t i = 0; i < g_xpcpData.size(); i++)
//	{
//		vector<pFlatPt*> oneFlats;
//		UncertainXPCPSample* uss = reinterpret_cast<UncertainXPCPSample*>(&g_xpcpData[i]);
//		uss->UncertainXPCPSample::to_pFlatRec(oneFlats, 1, UINT64(i));
//		// Add to global list
//		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
//		vector<pFlatPt*> twoFlats;
//		uss->UncertainXPCPSample::to_pFlatRec(twoFlats, 2, UINT64(i), g_use_repeat_pcp);
//		// Add to global list
//		g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
//	}
//	/////////////////////////////////////////////////////
//
//	// Delete raw kd tree?
//	if (_annKdTree_rawData)
//	{
//		SAFE_DELETE(_annKdTree_rawData);
//		_annKdTree_rawData = NULL;
//	}
//	////////////////////////////////////////////////
//
//	// Build KD trees for user query
//	g_0flat_list.resize(g_unique_sampl_raw.size());
//	for (size_t i = 0; i < g_unique_sampl_raw.size(); i++){
//		g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
//	}
//	if (g_kdTree_rawData){
//		g_kdTree_rawData->clear();
//		g_kdTree_rawData = NULL;
//	}
//	g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
//	g_kdTree_rawData->buildBalanceTree(g_0flat_list);
//
//	// Build KD tree for 1Flat
//	if (g_kdTree_1flat){
//		g_kdTree_1flat->clear();
//		g_kdTree_1flat = NULL;
//	}
//	g_kdTree_1flat = new KD<pFlatPt*>(2);
//	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
//	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
//	//ofOneFlatKdTree.close();
//	//////////////////////////////////////////////////
//	// Build KD tree for 2Flat
//	if (g_kdTree_2flat){
//		g_kdTree_2flat->clear();
//		g_kdTree_2flat = NULL;
//	}
//	g_kdTree_2flat = new KD<pFlatPt*>(2);
//	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
//	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
//	//ofTwoFlatKdTree.close();
//	/////////////////////////////////////////////////
//
//	SAFE_DELETE_ARRAY(nnIdx);
//	SAFE_DELETE_ARRAY(dists);
//}
//
//void XPCPdataAnalyzer::processVecRawData(const vector<vector<float>>& pointCenterData, const std::vector<std::vector<float>>& samplePos)
//{
//	if (pointCenterData.empty())
//	{
//		cout << "Empty raw data!" << endl;
//		return;
//	}
//	// Don't use nearest neighbors! Use vector values for p-flats directly!
//	//Set dimensions of the data
//	int dim_rawData = pointCenterData[0].size();
//	_numSamples = UINT64(pointCenterData.size());
//	_dimRawData = dim_rawData; // Raw data dimension
//	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
//	//Preperations: Clear list
//	g_xpcpData.resize(pointCenterData.size(), XPCPSample(dim_rawData));
//	// build pcp if necessary
//	if (_pcp == NULL)
//	{
//		vector<string> attribs(dim_rawData);
//		for (size_t i = 0; i < attribs.size(); i++)
//			attribs[i] = string("attr_") + number2String(int(i));
//		_pcp = new PCPInselberg(attribs);
//	}
//
//	// Allocate space for Major eigen data
//	g_majEigData.resize(pointCenterData.size());
//	// use global variables to keep the p-flat records
//	g_1flat_list.clear();
//	g_2flat_list.clear();
//
//	// Prepare variables to record min/max of p-flat strength in each subspace
//	_1flatsMinMaxPerSubspace.clear();
//	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
//	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
//	{
//		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//	_2flatsMinMaxPerSubspace.clear();
//	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
//	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
//	{
//		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//	ANNdist sqRad = 0.1 * 0.1;
//
//	// Compute 1-flat and 2-flat for each and every item in pointCenterData
//	ofstream ofDebug("neighborFitTest.txt");
//	for (UINT64 i = 0; i < pointCenterData.size(); i++)
//	{
//		// Just need to take the pointCenterData[i]
//		vector<float> vval = pointCenterData[i];
//		vector<float> pos = samplePos[i];
//		Eigen::VectorXf vecVval(vval.size());
//		for (size_t j = 0; j < vval.size(); j++)
//			vecVval[j] = vval[j];
//		XPCPSample xpcp_tuple(dim_rawData);
//		calc_p_flats_fromVecData(vval, pos, xpcp_tuple);
//
//		g_xpcpData[i] = xpcp_tuple;
//		g_majEigData[i] = vecVval; // Set Eigen vector data
//		///////////////////////////////////////////////////////////////////////
//		// !1 TODO: If we turn off normalization (!2), uncomment codes in the block
//		//// Convert p-flat data and store to global variables
//		//vector<pFlatPt*> oneFlats;
//		//xpcp_tuple.to_pFlatRec(oneFlats, 1, UINT64(i));
//		//// Add to global list
//		//g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
//		//vector<pFlatPt*> twoFlats;
//		//xpcp_tuple.to_pFlatRec(twoFlats, 2, UINT64(i));
//		//// Add to global list
//		//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
//		////////////////////////////////////////////////////////////////////////
//
//	}
//	ofDebug.close();
//	cout << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
//	/////////////////////////////////////////////////////
//	// !2 Sep 22, 2016: Added normalization func for p-flats 
//	// Do normalization after we compute all p-flats!
//	int normMethod = 0;
//	nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
//	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
//	for (size_t i = 0; i < g_xpcpData.size(); i++)
//	{
//		vector<pFlatPt*> oneFlats;
//		g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
//		// Add to global list
//		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
//		vector<pFlatPt*> twoFlats;
//		g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i));
//		// Add to global list
//		//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
//	}
//	////////////////////////////////////////////////
//	// Remove duplicate records
//	removeDuplicateRecords(pointCenterData);
//	/////////////////////////////////////////////
//	// TODO: comment out when test is done
//	/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
//	for (size_t i = 0; i < g_1flat_list.size(); i++)
//	ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
//	ofUniqueOneFlatTest.close();*/
//	////////////////////////////////////////////////
//
//	// Build KD trees for user query
//	g_0flat_list.resize(g_unique_sampl_raw.size());
//	for (size_t i = 0; i < g_unique_sampl_raw.size(); i++){
//		g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
//	}
//	if (g_kdTree_rawData){
//		g_kdTree_rawData->clear();
//		g_kdTree_rawData = NULL;
//	}
//	g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
//	g_kdTree_rawData->buildBalanceTree(g_0flat_list);
//
//	// Build KD tree for 1Flat
//	if (g_kdTree_1flat){
//		g_kdTree_1flat->clear();
//		g_kdTree_1flat = NULL;
//	}
//	g_kdTree_1flat = new KD<pFlatPt*>(2);
//	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
//	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
//	//ofOneFlatKdTree.close();
//	//////////////////////////////////////////////////
//	// Build KD tree for 2Flat
//	if (g_kdTree_2flat){
//		g_kdTree_2flat->clear();
//		g_kdTree_2flat = NULL;
//	}
//	g_kdTree_2flat = new KD<pFlatPt*>(2);
//	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
//	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
//	//ofTwoFlatKdTree.close();
//	/////////////////////////////////////////////////
//}
//
//
//void XPCPdataAnalyzer::processTrajRawData(const vector<vector<float>>& rawData)
//{
//	if (rawData.empty())
//	{
//		cout << "Empty raw data!" << endl;
//		return;
//	}
//	cout << "Processing raw data as trajectory data!" << endl;
//	// Don't use nearest neighbors! Use vector values for p-flats directly!
//	//Set dimensions of the data
//	int dim_rawData = rawData[0].size();
//	_numSamples = UINT64(rawData.size());
//	_dimRawData = dim_rawData; // Raw data dimension
//	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
//	//Preperations: Clear list
//	g_xpcpData.resize(rawData.size(), XPCPSample(dim_rawData));
//	// build pcp if necessary
//	if (_pcp == NULL)
//	{
//		vector<string> attribs(dim_rawData);
//		for (size_t i = 0; i < attribs.size(); i++)
//			attribs[i] = string("attr_") + number2String(int(i));
//		_pcp = new PCPInselberg(attribs);
//	}
//
//	// Allocate space for Major eigen data
//	g_majEigData.resize(rawData.size());
//	// use global variables to keep the p-flat records
//	g_1flat_list.clear();
//	g_2flat_list.clear();
//
//	// Prepare variables to record min/max of p-flat strength in each subspace
//	_1flatsMinMaxPerSubspace.clear();
//	_1flatsMinMaxPerSubspace.resize(dim_rawData - 1);
//	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
//	{
//		_1flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_1flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//	_2flatsMinMaxPerSubspace.clear();
//	_2flatsMinMaxPerSubspace.resize(dim_rawData - 2);
//	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
//	{
//		_2flatsMinMaxPerSubspace[i].x = numeric_limits<float>::max();
//		_2flatsMinMaxPerSubspace[i].y = -numeric_limits<float>::max();
//	}
//
//	// Compute 1-flat and 2-flat for each and every item in pointCenterData
//	ofstream ofDebug("neighborFitTest.txt");
//	for (UINT64 i = 0; i < rawData.size(); i++)
//	{
//		// Just need to take the pointCenterData[i]
//		vector<float> p1 = rawData[i];     // current point
//		vector<float> p2 = (i == rawData.size() - 1) ? p1 : rawData[i + 1]; // next point
//		// Assume linear connection between p1 and p2 now!
//		// TODO: switch to quadratic connection
//		XPCPSample xpcp_tuple(dim_rawData);
//		calc_p_flats_fromTrajData(p1, p2, xpcp_tuple);
//		Eigen::VectorXf tangent(p1.size());
//		for (size_t j = 0; j < p1.size(); j++)
//			tangent(j) = p2[j] - p1[j];
//		g_xpcpData[i] = xpcp_tuple;
//		g_majEigData[i] = tangent; // Set Eigen vector data
//
//
//	}
//	ofDebug.close();
//	cout << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
//	/////////////////////////////////////////////////////
//	// !2 Sep 22, 2016: Added normalization func for p-flats 
//	// Do normalization after we compute all p-flats!
//	int normMethod = 0;
//	nomralizePFlatsInSubspaces(g_xpcpData, normMethod);
//	// Set g_1flat_list & g_2flat_list. Comment (!1) when we use normalization
//	for (size_t i = 0; i < g_xpcpData.size(); i++)
//	{
//		vector<pFlatPt*> oneFlats;
//		g_xpcpData[i].to_pFlatRec(oneFlats, 1, UINT64(i));
//		// Add to global list
//		g_1flat_list.insert(g_1flat_list.end(), oneFlats.begin(), oneFlats.end());
//		vector<pFlatPt*> twoFlats;
//		g_xpcpData[i].to_pFlatRec(twoFlats, 2, UINT64(i));
//		// Add to global list
//		//g_2flat_list.insert(g_2flat_list.end(), twoFlats.begin(), twoFlats.end());
//	}
//	////////////////////////////////////////////////
//	// Remove duplicate records
//	removeDuplicateRecords(rawData);
//
//	// Build KD trees for user query
//	g_0flat_list.resize(g_unique_sampl_raw.size());
//	for (size_t i = 0; i < g_unique_sampl_raw.size(); i++){
//		g_0flat_list[i] = new zeroFlatPt(g_unique_sampl_raw[i]);
//	}
//	if (g_kdTree_rawData){
//		g_kdTree_rawData->clear();
//		g_kdTree_rawData = NULL;
//	}
//	g_kdTree_rawData = new KD<zeroFlatPt*>(dim_rawData);
//	g_kdTree_rawData->buildBalanceTree(g_0flat_list);
//
//	// Build KD tree for 1Flat
//	if (g_kdTree_1flat){
//		g_kdTree_1flat->clear();
//		g_kdTree_1flat = NULL;
//	}
//	g_kdTree_1flat = new KD<pFlatPt*>(2);
//	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
//	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
//	//ofOneFlatKdTree.close();
//	//////////////////////////////////////////////////
//	// Build KD tree for 2Flat
//	if (g_kdTree_2flat){
//		g_kdTree_2flat->clear();
//		g_kdTree_2flat = NULL;
//	}
//	g_kdTree_2flat = new KD<pFlatPt*>(2);
//	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
//	//////////////////////////////////////////////////
//	// TODO: Uncomment to write to file!
//	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
//	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
//	//ofTwoFlatKdTree.close();
//	/////////////////////////////////////////////////
//}

void XPCPdataAnalyzer::processRawDataInSubSpaces(const vector<vector<float>>& pointCenterData, int num_nearest_neighbors)
{
	if (pointCenterData.empty())
	{
		cout << "Empty raw data!" << endl;
		return;
	}
	int k_0 = num_nearest_neighbors;
	int dim_rawData = pointCenterData[0].size();
	//Set dimensions of the data
	_numSamples = UINT64(pointCenterData.size());
	_dimRawData = dim_rawData; // Raw data dimension
	_dimXPCPdata = dim_rawData + 2 * (dim_rawData - 1) + 2 * (dim_rawData - 2); // XPCP data dimension
	// 2D-KD tree
	ANNkd_tree* kdTree2D = NULL;
	// 3D-KD tree
	ANNkd_tree* kdTree3D = NULL;
	// Prepare query vairables
	// Do procssing
	double eps = 0.0; // error bound
	ANNpoint     queryPt2D, queryPt3D; // query points

	ANNidxArray  nnIdx;   // near neighbor indices
	ANNdistArray dists;   // near neighbor distances   

	queryPt2D = annAllocPt(2);
	queryPt3D = annAllocPt(3);

	nnIdx = new ANNidx[k_0];
	dists = new ANNdist[k_0];
	//Preperations: Clear list
	g_xpcpData.clear();
	g_xpcpData = vector<XPCPSample>(pointCenterData.size(), XPCPSample(_dimRawData));

	// build pcp if necessary
	if (_pcp == NULL)
	{
		vector<string> attribs(dim_rawData);
		for (size_t i = 0; i < attribs.size(); i++)
			attribs[i] = string("attr_") + number2String(int(i));
		_pcp = new PCPInselberg(attribs);
	}
	// Allocate space for Major eigen data
	g_majEigData.resize(pointCenterData.size(), Eigen::VectorXf(_dimRawData));
	// use global variables to keep the p-flat records
	g_1flat_list.clear();
	g_2flat_list.clear();

	// Compute 1-flat and 2-flat for each and every item in pointCenterData
	ofstream ofDebug("neighborFitTest.txt");
	/////////////////////////////////////////////////
	// Calucluate 1-flats!!!
	// 1. Outer loop: subspaces
	int ssDim = 2;
	for (int ssId = 0; ssId < _dimRawData - 1; ssId++)
	{
		// build 2D tree
		buildKDTreeInSubSpaces(pointCenterData, ssId, ssDim, &kdTree2D);
		// 2. Inner loop: samples
		// 2.1 Compute 1-flats (2D subspace)
		int k = k_0;
		// neighboring samples
		vector<vector<float>> neighbors(k);
		// Search for nearest neighbors
		for (UINT64 i = 0; i < UINT64(pointCenterData.size()); i++)
		{
			// Compute 1-flat
			for (int dd = 0; dd < ssDim; dd++)
				queryPt2D[dd] = pointCenterData[i][ssId + dd];

			kdTree2D->annkSearch(// search
				queryPt2D,	// query point
				k,			// number of near neighbors
				nnIdx,		// nearest neighbors (returned)
				dists,		// distance (returned)
				eps			// error bound
				);

			// Get all points for computation
			//---------------------------------------------------------------
			// From ANN manual page 9:
			//https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf
			// the neighbors are sorted by their distances to the query point: i.e., nnIdx[0] is the neareast one, in our case, the query point itself!
			for (int nn = 0; nn < k; nn++)
			{
				int id = nnIdx[nn];
				neighbors[nn] = vector<float>(ssDim);
				for (int dd = 0; dd < ssDim; dd++)
					neighbors[nn][dd] = pointCenterData[id][ssId + dd];
			}
			Eigen::VectorXf mu;
			MyStatistics::ComputeMean<float>(neighbors, mu);
			if (i == 0)
			{
				cout << "Mean = " << endl << mu << endl;
				ofDebug << "Mean = " << endl << mu << endl;
			}
			// Compute the covariance matrix and then the SVD of that cov matrix
			Eigen::MatrixXf svdU, svdV;
			Eigen::VectorXf svdS;
			MyStatistics::SVDCovMatrix(neighbors, svdU, svdS, svdV);
			Eigen::VectorXf majSingV = svdU.col(0) * svdS(0);
			Eigen::VectorXf secMajSingV = svdU.col(1) * svdS(1);
			if (i == 0)
			{
				cout << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
				ofDebug << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
			}
			Eigen::VectorXf majEigV = majSingV;
			Eigen::VectorXf secMajEigV = secMajSingV;
			XPCPSample xpcp_tuple;
			// compute p-flat
			FLOATVECTOR2 fOneFlat;
			calc1Flat(mu, majEigV, ssId, fOneFlat);
			// Setup output
			g_xpcpData[i].pcp_1flat[ssId] = fOneFlat;
			g_majEigData[i][ssId] = majEigV(0); // Set the ssId dimension of the Eigen vector data
			g_majEigData[i][ssId + 1] = majEigV(1);
			// Add to global list
			g_1flat_list.insert(g_1flat_list.end(), new pFlatPt(fOneFlat.x, fOneFlat.y, ssId, i));
		}

	}
	////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////
	// Calculate 2-flats
	// 1. Outer loop: subspaces
	ssDim = 3;
	for (int ssId = 0; ssId < _dimRawData - 2; ssId++)
	{
		// build 3D tree
		buildKDTreeInSubSpaces(pointCenterData, ssId, ssDim, &kdTree3D);
		// 2. Inner loop: samples
		// 2.1 Compute 2-flats (3D subspace)
		int k = k_0;
		// neighboring samples
		vector<vector<float>> neighbors(k);
		// Search for nearest neighbors
		for (UINT64 i = 0; i < UINT64(pointCenterData.size()); i++)
		{
			// Compute 1-flat
			for (int dd = 0; dd < ssDim; dd++)
				queryPt3D[dd] = pointCenterData[i][ssId + dd];

			kdTree3D->annkSearch(// search
				queryPt3D,	// query point
				k,			// number of near neighbors
				nnIdx,		// nearest neighbors (returned)
				dists,		// distance (returned)
				eps			// error bound
				);

			// Get all points for computation
			//---------------------------------------------------------------
			// From ANN manual page 9:
			//https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ANNmanual_1.1.pdf
			// the neighbors are sorted by their distances to the query point: i.e., nnIdx[0] is the neareast one, in our case, the query point itself!
			for (int nn = 0; nn < k; nn++)
			{
				int id = nnIdx[nn];
				neighbors[nn] = vector<float>(ssDim);
				for (int dd = 0; dd < ssDim; dd++)
					neighbors[nn][dd] = pointCenterData[id][ssId + dd];
			}
			// Compute mean vector
			Eigen::VectorXf mu;
			MyStatistics::ComputeMean<float>(neighbors, mu);
			if (i == 0)
			{
				cout << "Mean = " << endl << mu << endl;
				ofDebug << "Mean = " << endl << mu << endl;
			}
			// Compute the covariance matrix and then the SVD of that cov matrix
			Eigen::MatrixXf svdU, svdV;
			Eigen::VectorXf svdS;
			MyStatistics::SVDCovMatrix(neighbors, svdU, svdS, svdV);
			Eigen::VectorXf majSingV = svdU.col(0) * svdS(0);
			Eigen::VectorXf secMajSingV = svdU.col(1) * svdS(1);
			if (i == 0)
			{
				cout << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
				ofDebug << "Maj Sing Vec: " << endl << majSingV << endl << "Sec Maj Sing Vec: " << endl << secMajSingV << endl;
			}
			Eigen::VectorXf majEigV = majSingV;
			Eigen::VectorXf secMajEigV = secMajSingV;
			XPCPSample xpcp_tuple;
			// compute p-flat
			FLOATVECTOR2 fTwoFlat;
			calc2Flat(mu, majEigV, secMajEigV, ssId, fTwoFlat);
			// Setup output
			g_xpcpData[i].pcp_2flat[ssId] = fTwoFlat;
			// Add to global list
			g_2flat_list.insert(g_2flat_list.end(), new pFlatPt(fTwoFlat.x, fTwoFlat.y, ssId, i));
			////////////////////////////////////////////////
		}

	}
	cout << "Number of entries in XPCP data: " << g_xpcpData.size() << endl;
	// Delete raw kd tree?
	if (_annKdTree_rawData)
	{
		SAFE_DELETE(_annKdTree_rawData);
		_annKdTree_rawData = NULL;
	}
	//// Output major eigen vectors
	//ofstream ofMajEig("majEig.txt");
	//for (size_t i = 0; i < g_majEigData.size(); i++)
	//	ofMajEig << g_majEigData[i] << endl;
	//ofMajEig.close();
	////////////////////////////////////////////////
	// TODO: comment out the following until the separators when test is done!
	//// Write result of 1-flats
	//ofstream ofOneFlatTest("oneFlatTest.txt");
	//for (size_t i = 0; i < g_1flat_list.size(); i++)
	//	ofOneFlatTest << *(g_1flat_list[i]) << endl;
	//ofOneFlatTest.close();
	////////////////////////////////////////////////
	// Remove duplicate records
	removeDuplicateRecords(pointCenterData);
	/////////////////////////////////////////////
	// TODO: comment out when test is done
	/*ofstream ofUniqueOneFlatTest("uniqOneFlatTest.txt");
	for (size_t i = 0; i < g_1flat_list.size(); i++)
	ofUniqueOneFlatTest << *(g_1flat_list[i]) << endl;
	ofUniqueOneFlatTest.close();*/
	////////////////////////////////////////////////

	// Build KD trees for p-Flats
	// Build KD tree for 1Flat
	if (g_kdTree_1flat){
		g_kdTree_1flat->clear();
		g_kdTree_1flat = NULL;
	}
	g_kdTree_1flat = new KD<pFlatPt*>(2);
	g_kdTree_1flat->buildBalanceTree(g_1flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofOneFlatKdTree("oneFlatKdtree.txt");
	//g_kdTree_1flat->printToFile(ofOneFlatKdTree);
	//ofOneFlatKdTree.close();
	//////////////////////////////////////////////////
	// Build KD tree for 2Flat
	if (g_kdTree_2flat){
		g_kdTree_2flat->clear();
		g_kdTree_2flat = NULL;
	}
	g_kdTree_2flat = new KD<pFlatPt*>(2);
	g_kdTree_2flat->buildBalanceTree(g_2flat_list);
	//////////////////////////////////////////////////
	// TODO: Uncomment to write to file!
	//ofstream ofTwoFlatKdTree("twoFlatKdtree.txt");
	//g_kdTree_2flat->printToFile(ofTwoFlatKdTree);
	//ofTwoFlatKdTree.close();
	/////////////////////////////////////////////////

	SAFE_DELETE_ARRAY(nnIdx);
	SAFE_DELETE_ARRAY(dists);
}

void XPCPdataAnalyzer::calc_p_flats_fromLowerDimMethod(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, XPCPSample& xpcp_tuple)
{
	//// The normal of the hyperplane
	//Eigen::VectorXf hyperPlaneNormal = majEigV.cross(secMajEigV);
	Eigen::VectorXf ptMajEig = majEig + mu;
	Eigen::VectorXf ptSecEig = secEig + mu;
	vector<float> vPtMajEig(ptMajEig.data(), ptMajEig.data() + ptMajEig.rows());
	vector<float> vPtSecEig(ptSecEig.data(), ptSecEig.data() + ptSecEig.rows());
	vector<float> vPtMu(mu.data(), mu.data() + mu.rows());
	// convert to pcp space

	vector<float> v_Eig_Cart0(majEig.data(), majEig.data() + majEig.rows());
	vector<float> v_Eig_Cart1(secEig.data(), secEig.data() + secEig.rows());
	// 0. 
	// Convert points to PCP coordinates: PCP coordinates with no repeat
	vector<Point> pcp_p0_no_repeat;
	vector<Point> pcp_p1_no_repeat;
	vector<Point> pcp_p2_no_repeat;
	_pcp->cartesian2PCP(vPtMu, pcp_p0_no_repeat);
	_pcp->cartesian2PCP(vPtMajEig, pcp_p1_no_repeat);
	_pcp->cartesian2PCP(vPtSecEig, pcp_p2_no_repeat);
	// 1. 
	// Get 1-flat for major and the second eigen vectors
	vector<Point> pcp_v0_no_repeat;
	vector<Point> pcp_v1_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p2_no_repeat, pcp_v1_no_repeat);
	// 2. compute 2-flat
	vector<Point> flat2;
	bool has2flat = true;
	if (!_pcp->calcFlatOneOrderHigher(pcp_v0_no_repeat, pcp_v1_no_repeat, flat2))
		has2flat = false; //  if something went wrong

	// set output
	xpcp_tuple.x = vPtMu;
	// 1-flat
	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(pcp_v0_no_repeat[i].x, pcp_v0_no_repeat[i].y);
	// store 2-flat
	if (has2flat){
		xpcp_tuple.pcp_2flat.resize(flat2.size());
		for (size_t i = 0; i < flat2.size(); i++)
			xpcp_tuple.pcp_2flat[i] = FLOATVECTOR2(flat2[i].x, flat2[i].y);
	}

}

void XPCPdataAnalyzer::calc_p_flats_fromLowFlats(const Eigen::VectorXf& mu, const vector<Point>& pcp_p0_no_repeat, const vector<Point>& pcp_p1_no_repeat, const vector<Point>& pcp_p2_no_repeat, XPCPSample& xpcp_tuple)
{
	// 1. 
	// Get 1-flat for major and the second eigen vectors
	vector<Point> pcp_v0_no_repeat;
	vector<Point> pcp_v1_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p2_no_repeat, pcp_v1_no_repeat);
	// 2. compute 2-flat
	vector<Point> flat2;
	bool has2flat = true;
	if (!_pcp->calcFlatOneOrderHigher(pcp_v0_no_repeat, pcp_v1_no_repeat, flat2))
		has2flat = false; //  if something went wrong

	// set output
	vector<float> vPtMu(mu.data(), mu.data() + mu.rows());
	xpcp_tuple.x = vPtMu;
	// 1-flat
	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(pcp_v0_no_repeat[i].x, pcp_v0_no_repeat[i].y);
	// store 2-flat
	if (has2flat) {
		xpcp_tuple.pcp_2flat.resize(flat2.size());
		for (size_t i = 0; i < flat2.size(); i++)
			xpcp_tuple.pcp_2flat[i] = FLOATVECTOR2(flat2[i].x, flat2[i].y);
	}
}

void XPCPdataAnalyzer::calc_p_flats_subSpaceMethod(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, XPCPSample& xpcp_tuple)
{
	int N = mu.size(); // dimension of the input data

	Eigen::VectorXf ptMajEig = majEig + mu;
	vector<float> vPtMajEig(ptMajEig.data(), ptMajEig.data() + ptMajEig.rows());
	vector<float> vPtMu(mu.data(), mu.data() + mu.rows());
	// 0. Set N-D value
	xpcp_tuple.x = vPtMu;
	// threshold the length of the major eigen vector
	float twoFlat_discard_thres = 1e-3f;
	float oneFlat_discard_thres = 1e-3f;

	// 1. Compute 1-flat
	vector<Point> pcp_v0_no_repeat;

	if (g_params.XformMethod() == XF_NGUYEN_ROSEN)
	{
		_pcp->calcOneFlatParametricForm(vPtMu, vPtMajEig, pcp_v0_no_repeat);
	}
	else
		_pcp->calcOneFlatGeneralForm(vPtMu, vPtMajEig, pcp_v0_no_repeat);
	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());

	double u = 0.0;
	double v = 0.0;
	// for nguyen & rosen's tvcg 2017
	////////////////////////////////////
	double uu = 0.0;
	double vv = 0.0;
	double x0 = 0.0;
	double y0 = 0.0;
	bool isRotated = false;
	///////////////////////////////////
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
	{
		Point idxPt = pcp_v0_no_repeat[i]; // 3D line coordinates0

		if (g_params.XformMethod() == XF_IDX_ONLY || g_params.XformMethod() == XF_IDX_LINE)
		{
			double c1 = idxPt.x;
			double c2 = idxPt.y;
			double c3 = idxPt.z;
			if (idxPt.y != 0.0)
			{
				double m = -idxPt.x / idxPt.y;
				double b = -idxPt.z / idxPt.y;

				c1 = -m;
				c2 = 1.0;
				c3 = -b;
			}

			if(g_params.IsScaleXform())
				xform_idxpt_complete_line_desc(c1, c2, c3, u, v);
			else
				xform_idxpt_complete_line_desc_noScaling(c1, c2, c3, u, v);
		}
		else if (g_params.XformMethod() == XF_NGUYEN_ROSEN)
		{
			x0 = vPtMu[i];
			y0 = vPtMu[i + 1];
			uu = idxPt.x;
			vv = idxPt.y;
			xform_nguyen_rosen_tvcg(uu, vv, x0, y0, u, v, isRotated);
		}
		else
			lineCoord2pointCoord2D(idxPt.x, idxPt.y, idxPt.z, u, v);
		// Adding the offset of the dimension
		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(float(u) + float(i), float(v));
		float eigVLen = Eigen::Vector2f(majEig[i], majEig[i + 1]).norm();
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
		xpcp_tuple.strength_1flat[i].x = eigVLen;
		if (g_params.XformMethod() == XF_NGUYEN_ROSEN) //HACK: use the y component of the strength to encode whether the point is rotated for Nguyen&Rosen!!!
			xpcp_tuple.isRotated_1flat[i] = isRotated;
		//**********************************************************
		// Alternatively, use the correlation value for 1-flat strength
		// NOTE: Uncomment to use correlation value as strength!
		//	xpcp_tuple.strength_1flat[i].x = abs(corr_mat(i,i+1));
		//**********************************************************
		float minEigVlen = MIN(_1flatsMinMaxPerSubspace[i].x, xpcp_tuple.strength_1flat[i].x);
		_1flatsMinMaxPerSubspace[i].x = minEigVlen;
		float maxEigVlen = MAX(_1flatsMinMaxPerSubspace[i].y, xpcp_tuple.strength_1flat[i].x);
		_1flatsMinMaxPerSubspace[i].y = maxEigVlen;
	}

	//=======================================================================			
	// 2. Compute 2-flat
	// For now, use 3-tuple in sequency of the axis ordering
	if (secEig.size() == 0)
		return; // no second eigen vector, we're done!
	vector<FLOATVECTOR2> indexPt;
	vector<FLOATVECTOR3> indexPtGf; // general form indexed points
	bool oneIndex = !g_use_repeat_pcp;
	int numIdxPerPlane = oneIndex ? 1 : 4;
	xpcp_tuple.pcp_2flat.resize((N - 2) * numIdxPerPlane);
	int normMethod = 1; // 0: length of the normal vector, 1: product of the lengths of two tangent vectors, 2: the larger of the tangents, 3: the smaller of the the tangents
	for (int j = 0; j < N - 2; j++)
	{
		// Get 3D subspace!!! Use projected eigenvectors as tangent vectors of the plane
		Eigen::Vector3f mu3d = mu.block(j, 0, 3, 1);
		Eigen::Vector3f majEig3d = majEig.block(j, 0, 3, 1);
		Eigen::Vector3f secEig3d = secEig.block(j, 0, 3, 1);
		//majEig3d.normalize(); // normalize these tangent vectors
		//secEig3d.normalize();

		// Do analytical computation to find 2-flat for the subspace
		Eigen::Vector3f norm3d = majEig3d.cross(secEig3d);
		float strength_2flat = 0.0f;
		if (norm3d == Eigen::Vector3f(0, 0, 0))
		{
			xpcp_tuple.strength_2flat[j].x = strength_2flat;
			xpcp_tuple.strength_2flat[j].y = strength_2flat;

			_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
			_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);
			continue;
		}

		/////////////////////////////////////

		// Normalize the normal vector
		strength_2flat = computeNormalVecStrength(majEig3d, secEig3d, norm3d, normMethod);
		//strength_2flat = norm3d.norm();
		//////////////////////////////////////////

		_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
		_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);

		PlaneCoeffs Pl;
		// The plane needs a normalized normal vector
		norm3d.normalize();
		Pl.c1 = norm3d.x(); // c1--c3 are components of the normal vector
		Pl.c2 = norm3d.y();
		Pl.c3 = norm3d.z();
		Pl.c0 = norm3d.x() * mu3d[0] + norm3d.y() * mu3d[1] + norm3d.z() * mu3d[2]; // c0 is the min distance from the origin to the plane

		// Compute indexed point(s)
		_pcp->calc2FlatFrom3DPlane(Pl, indexPt, oneIndex);
		//_pcp->calc2FlatFrom3DPlane_GeneralForm(Pl, indexPtGf, oneIndex);

		double u2 = 0.0;
		double v2 = 0.0;
		// Set 2-flat
		for (int k = 0; k < numIdxPerPlane; k++)
		{
			// Conduct the same transformation as 1-flats
			//FLOATVECTOR3 idxPt = indexPt[k]; 
			//FLOATVECTOR3 idxPt = indexPtGf[k];
			if (g_params.XformMethod() == XF_IDX_ONLY || g_params.XformMethod() == XF_IDX_LINE)
			{
				//double c1 = idxPt.x;
				//double c2 = idxPt.y;
				//double c3 = idxPt.z;
				//if (idxPt.y != 0.0)
				//{
				//	double m = -idxPt.x / idxPt.y;
				//	double b = -idxPt.z / idxPt.y;

				//	c1 = -m;
				//	c2 = 1.0;
				//	c3 = -b;
				//}

				//if (g_params.IsScaleXform())
				//	xform_idxpt_complete_line_desc(c1, c2, c3, u, v);
				//else
				//	xform_idxpt_complete_line_desc_noScaling(c1, c2, c3, u, v);

				xform(indexPt[k].x, indexPt[k].y, j, u2, v2);

				//if (g_params.IsScaleXform())
				//	xform_idxpt_complete_line_desc(indexPt[k].x, indexPt[k].y, j, u, v);
				//else
				//	xform_idxpt_complete_line_desc_noScaling(indexPt[k].x, indexPt[k].y, j, u, v);
			}

			// Shift to the starting dimension


			//indexPt[k].x += float(j) + 0.5f;
			u2 += float(j) + 0.5f;
			xpcp_tuple.pcp_2flat[numIdxPerPlane * j + k] = FLOATVECTOR2(u2, v2); // indexPt[k];
			xpcp_tuple.strength_2flat[j].x = strength_2flat; // set strength
		}

	}
}


void XPCPdataAnalyzer::xform(double xin, double yin, int subspaceid, double& xout, double& yout)
{
	int option = 0;
	static bool firstTime = true;

	if (xin == 0.0) // on the left axis
	{
		yout = yin;
		xout = xin;
		return;
	}
	double a = 1.0 - 1.0 / xin;
	double theta = atan(1.0 - 1.0 / xin);
	if (theta == M_PI / 4) // pi/4, the ideal point
	{
		xout = (g_randGen.rand() < 0.5) ? -0.5 : 1.5; // half chance to the left, half to the right
		yout = 0.0;
		return;
	}

	float scale = 1.0f;
	double delta = 1e-3;
	double s0 = delta / (atan(1.0 - 1.0 / delta) * 2.0 / M_PI + 1.0);
	double s1 = -delta / (atan(1.0 + 1.0 / delta) * 2.0 / M_PI - 1.0);
	if (firstTime)
	{
		cout << "S0 = " << s0 << ", S1 = " << s1 << endl;
		firstTime = false;
	}
	double a0 = 1.0 / delta;
	double t = 1.0;
	switch (option)
	{
	case 0:

		if (theta > M_PI / 4) // left of the axes pair
			xout = theta * 2 / M_PI - 1;
		else // For theta < 0 && theta < pi/4
			xout = theta * 2 / M_PI + 1.0;

		if (theta == 0.0)
			yout = yin; // if on either axes, use input y value
		else
		{
			// preserves linearity for y
			yout = (xout - 0.5) * yin / (xin - 0.5);
			//====================================
			//// NOTE:
			//// Uncomment for No scaling!!!
			return;
			/////////////////////////////////////

			// Scaling for y!


			if (theta < 0.0f){ // inside the pair
				t = xin;
				//scale = xin / (atan(1.0 - 1.0/xin)*2.0/M_PI + 1.0);
				scale = s0 * (1.0 - t) + t;
			}
			else if (theta > M_PI / 4) // left of the pair
			{

				t = -xout / 0.5; // normalized distance to the origin (left axis)
				float s1 = 1.1f * s0;
				float anchorpoint = 1.0f;
				if (t < anchorpoint)
				{
					t /= anchorpoint;
					scale = t * s1 + (1 - t) * s0; // left scale s1, right scale s0
				}

			}
			else if (theta > 0 && theta < M_PI / 4) // right of the pair
			{
				//scale = 1.0;
				t = (xout - 1.0) / 0.5;
				scale = 1.0 / s0 * t + (1.0 - t) * 1.0;
			}

			yout *= scale;
			//////////////////////////////////////
		}

		break;
	case 1:
		// 3. derived dx/dtheta = c transformation Dec.2.2016

		if (theta > M_PI / 4.0f) // left
			xout = (theta - M_PI / 2) * 2 / M_PI;
		else if (theta < 0.0f)// inside
		{
			xout = xin;
		}
		else // right
			xout = theta * 2 / M_PI + 1.0;
		if (yin < 0.0f)
			yout = tanh(double(yin)) * 0.5; // when y < 0, range of v is (-0.5,0)
		else if (yin > 1.0f)
			yout = tanh(double(yin) - 1.0) * 0.5 + 1.0; // when y > 1, range of v is (1,1.5)
		else
			yout = yin; // when 0<=y<=1, no changes
		break;
	}
}

void XPCPdataAnalyzer::xform_idxpt_complete_line_desc_origin0_noScaling(double c1, double c2, double c3, double& xout, double& yout)
{
	if (c2 == 0) // +/- pi/2; at x1 axis
	{
		xout = 0.0;
		yout = -c3 / c1;
		return;
	}

	if (c1 == -c2) // pi/4
	{
		xout = (g_randGen.rand() < 0.5) ? -0.5 : 1.5; // half chance to the left, half to the right
		yout = 0.0;
		return;
	}

	// parameters of the associated Cartesian line 
	double a = -c1 / c2;
	double b = -c3 / c2;

	double theta = atan(a);

	if (theta > M_PI / 4) // left of the axes pair
	{
		xout = theta * 2 / M_PI - 1;
	}
	else // For theta < 0 && theta < pi/4
	{
		xout = theta * 2 / M_PI + 1;
	}
	// Old transformation method with the origin at 0
	//yout = xout * b * scale; // output y
	// with new offset
	yout = xout * b; // output y
}

void XPCPdataAnalyzer::xform_nguyen_rosen_tvcg(double u, double v, double x0, double y0, double& q, double& r, bool& rotated)
{
	if (u*v < 0)
	{
		rotated = false;
		q = u / (u - v);
		r = x0 + (y0 - x0)*q;
	}
	else
	{
		rotated = true;
		q = -v / (-v - u);
		r = x0 + (y0 - x0)*q;
	}



}

bool compFirstElem(const DOUBLEVECTOR2& a, const DOUBLEVECTOR2& b)
{
	return a.x < b.x;
}

bool XPCPdataAnalyzer::loadScaleFactors(const string& fileName, vector<double>& scaleLocs, vector<double>& scaleFactors)
{
	ifstream ifSF(fileName.c_str());
	if (ifSF.is_open())
	{
		vector<DOUBLEVECTOR2> scaleVals;
		while (!ifSF.eof())
		{
			DOUBLEVECTOR2 tuple;
			ifSF >> tuple;
			scaleVals.push_back(tuple);
		}

		sort(scaleVals.begin(), scaleVals.end(), compFirstElem);
		scaleLocs.resize(scaleVals.size());
		scaleFactors.resize(scaleVals.size());
		for (size_t i = 0; i < scaleVals.size(); i++)
		{
			scaleLocs[i] = scaleVals[i].x;
			scaleFactors[i] = scaleVals[i].y;
		}
		ifSF.close();
		return true;
	}
	else
		return false;
}

// Transformation from traditional parallel coordinates to angle-uniform parallel coordinates
void XPCPdataAnalyzer::xform_idxpt_complete_line_desc(double c1, double c2, double c3, double& u, double& v)
{
	//=========================================================================
	// NOTE: for pi/4, we use the limit approach by taking an indexed point of a Cartesian line 
	// with the orientation angle close but not equal to pi/4, such that it has valid input parameters.
	//=========================================================================
	if (c2 == 0) // for +/- pi/2; the transformed point is at x1 axis
	{
		u = 0.0;
		v = -c3 / c1;
		return;
	}

	// parameters of the associated Cartesian line 
	double a = -c1 / c2;
	double b = 2 * c3 / (c1 - c2);  //-c3 / c2;

	// the angle of the orientation
	double theta = atan(a);
	// the scaling factor
	double scale = 1.0;

	if (theta > M_PI / 4) // left of the axes pair
	{
		u = theta * 2 / M_PI - 1;
		scale = spline_interp(u, spline_coefs, spline_breaks); // evaluate the cubic spline
		//scale = 1.0;
	}
	else // For theta < 0 && theta < pi/4
	{
		u = theta * 2 / M_PI + 1;
	
		if (theta < 0) // inside the pair
		{
			
			// Use the exact scaling for traditional parallel coordinates
			//scale = (1.0/(1.0 -a) - 0.5) / (u - 0.5) ;
			// Use the smooth scaling 
			scale = spline_interp(u, spline_coefs, spline_breaks); //(1.0/(1.0 -a) - 0.5) / (u - 0.5) ;
		}
		else
		{ 	// right of the pair
			//scale = 1.0;
			scale = spline_interp(u, spline_coefs, spline_breaks); // evaluate the cubic spline
		}
	}
	v = (u - 0.5) * b * scale; // the vertical coordinate
}

//
//void XPCPdataAnalyzer::xform_idxpt_complete_line_desc(double c1, double c2, double c3, double& xout, double& yout)
//{
//	if (scaleLocs.empty() || scaleFactors.empty()){
//		if (!loadScaleFactors(defautSFfileName, scaleLocs, scaleFactors))
//		{
//			xout = 0.0;
//			yout = 0.0;
//			return;
//		}
//	}
//		
//	if (c2 == 0) // +/- pi/2; at x1 axis
//	{
//		xout = 0.0;
//		yout = -c3 / c1;
//		return;
//	}
//
//	if ( c1 == -c2) // pi/4
//	{
//		xout = (g_randGen.rand() < 0.5) ? -0.5 : 1.5; // half chance to the left, half to the right
//		yout = 0.0;
//		return;
//	}
//	
//	// parameters of the associated Cartesian line 
//	double a = -c1 / c2;
//	double b = -c3 / c2;
//
//	double theta = atan(a);
//	double scale = 1.0;
//	//double delta = 1e-3;
//	//
//	//double s0 = delta / (atan(1.0 - 1.0 / delta) * 2.0 / M_PI + 1.0);
//	//double s1 = 1.1 * s0;
//	//double a0 = 1.0 / delta;
//	double t = 1.0;
//
//	double bS = 1.0;
//	if (theta > M_PI / 4) // left of the axes pair
//	{
//		xout = theta * 2 / M_PI - 1;
//		if (g_isScaleXform)
//		{
//			// get the scale from the right hand-side
//			t = (-xout) / 0.5;
//			bS = g_boundScale * t + (1.0 - t) * 1.0;
//
//			//scale = poly_func(polyScalingParams, xout);
//			scale = spline_interp(xout, -0.5, spline_coefs, spline_breaks);
//			scale *= bS;
//			////// scaling factor
//			////t = -xout / 0.5;
//			////scale = t * s1 + (1 - t)*s0;
//			//vector<double>::const_iterator Upit = std::lower_bound(scaleLocs.begin(), scaleLocs.end(), xout);
//			//// Lowerbound gives the first element that is not less than value!!!
//			//vector<double>::const_iterator Doit = std::upper_bound(scaleLocs.begin(), scaleLocs.end(), xout);
//			//// Upperbound gives the first element that is greater than value
//			//if (Upit == scaleLocs.end())
//			//{
//			//	scale = scaleFactors.back();
//			//}
//			//else if (Doit == scaleLocs.begin())
//			//{
//			//	scale = scaleFactors[0]; // 
//			//}
//			//else
//			//{
//			//	if (*Upit == xout)
//			//		scale = scaleFactors[Upit - scaleLocs.begin()];
//			//	else
//			//	{
//			//		// blend left and right elements
//			//		int up = Upit - scaleLocs.begin() - 1;
//			//		scale = scaleFactors[up];
//			//		//int down = Doit - scaleLocs.begin();
//			//		////scale = scaleFactors[down];
//			//		//double oneOver = 1.0 / abs(scaleLocs[up] - scaleLocs[down]);
//			//		//scale = Lerp((xout - scaleLocs[up])*oneOver, scaleFactors[down], scaleFactors[up]);
//			//	}
//			//}
//			//int left = std::distance(scaleLocs.begin(), std::lower_bound(scaleLocs.begin(), scaleLocs.end(), xout));
//			//scale = scaleFactors[left];
//	/*		int right = std::upper_bound(scaleLocs.begin(), scaleLocs.end(), xout) - scaleLocs.begin();
//			float oneOver = 1.0f / abs(scaleLocs[left] - scaleLocs[right]);
//			float tt = (xout - scaleLocs[left]) * oneOver;
//			scale = Lerp(tt, scaleFactors[left], scaleFactors[right]);*/
//		}
//	}
//	else // For theta < 0 && theta < pi/4
//	{
//		xout = theta * 2 / M_PI + 1;
//		if (g_isScaleXform)
//		{
//			if (theta < 0) // inside the pair
//			{
//				//double xin = 1.0 / (1 - a);
//				//t = xin;
//				//scale = s0 * (1.0 - t) + t;
//				//if (xout < 1.0 /*u_vertical_max*/)
//				{
//					//scale = poly_func(polyScalingParams, xout);
//					scale = spline_interp(xout, -0.5, spline_coefs, spline_breaks);
//					//vector<double>::const_iterator Upit = std::lower_bound(scaleLocs.begin(), scaleLocs.end(), xout);
//					//vector<double>::const_iterator Doit = std::upper_bound(scaleLocs.begin(), scaleLocs.end(), xout);
//					//if (Upit == scaleLocs.end())
//					//{
//					//	scale = scaleFactors.back();
//					//}
//					//else if (Doit == scaleLocs.begin())
//					//{
//					//	scale = scaleFactors[0]; // 
//					//}
//					//else
//					//{
//					//	if (*Upit == xout)
//					//		scale = scaleFactors[Upit - scaleLocs.begin()];
//					//	else
//					//	{
//					//		// blend left and right elements
//					//		int up = Upit - scaleLocs.begin() - 1;
//					//		scale = scaleFactors[up];
//					//		//int down = Doit - scaleLocs.begin();
//					//		////scale = scaleFactors[down];
//					//		//double oneOver = 1.0 / abs(scaleLocs[up] - scaleLocs[down]);
//					//		//scale = Lerp((xout - scaleLocs[up])*oneOver, scaleFactors[down], scaleFactors[up]);
//					//	}
//					//}
//				}
//			/*	else
//					scale = 1.0;*/
//			}
//			else
//			{ // right of the pair
//				t = (xout - 1.0) / 0.5;
//				scale = g_boundScale * t + (1.0 - t) * 1.0;
//
//				//scale = 1.0;
//			}
//		}
//	}
//	yout = xout * b * scale; // output y
//}

void XPCPdataAnalyzer::xform_idxpt_complete_line_desc_noScaling(double c1, double c2, double c3, double& xout, double& yout)
{

	if (c2 == 0) // +/- pi/2; at x1 axis
	{
		xout = 0.0;
		yout = -c3 / c1;
		if (c1 == 0.0)
			yout = 0.0;
		return;
	}

	if (c1 == -c2) // pi/4
	{
		xout = (g_randGen.rand() < 0.5) ? -0.5 : 1.5; // half chance to the left, half to the right
		yout = 0.0;
		return;
	}

	// parameters of the associated Cartesian line 
	double a = -c1 / c2;
	double b = 2 * c3 / (c1 - c2);  //-c3 / c2;

	double theta = atan(a);
	double scale = 1.0;

	if (theta > M_PI / 4) // left of the axes pair
	{
		xout = theta * 2 / M_PI - 1;
	}
	else // For theta < 0 && theta < pi/4
	{
		xout = theta * 2 / M_PI + 1;
	}
	// Old transformation method with the origin at 0
	//yout = xout * b * scale; // output y
	// with new offset
	yout = (xout - 0.5f) * b * scale; // output y
}

void XPCPdataAnalyzer::xform_idxpt_complete_line_desc_origin0(double c1, double c2, double c3, double& u, double& v)
{
	//=========================================================================
	// NOTE: for pi/4, we use the limit approach by taking an indexed point of a Cartesian line 
	// with the orientation angle close but not equal to pi/4, such that it has valid input parameters.
	//=========================================================================
	if (c2 == 0) // for +/- pi/2; the transformed point is at x1 axis
	{
		u = 0.0;
		v = -c3 / c1;
		return;
	}

	// parameters of the associated Cartesian line 
	double a = -c1 / c2;
	double b = -c3 / c2;

	// the angle of the orientation
	double theta = atan(a);
	// the scaling factor
	double scale = 1.0;

	if (theta > M_PI / 4) // left of the axes pair
	{
		u = theta * 2 / M_PI - 1;
		//scale = spline_interp(u, -0.5, spline_coefs, spline_breaks); // evaluate the cubic spline
		scale = 1.0;
	}
	else // For theta < 0 && theta < pi/4
	{
		u = theta * 2 / M_PI + 1;

		if (theta < 0) // inside the pair
		{
			scale = spline_interp(u, spline_coefs, spline_breaks);// evaluate the cubic spline
			//scale = poly_func(polyScalingParams, u);
		}
		else
		{ 	// right of the pair
			scale = 1.0;
		}
	}
	v = u * b * scale; // the vertical coordinate
}


void XPCPdataAnalyzer::lineCoord2pointCoord2D(double c1, double c2, double c3, double& x, double& y)
{
	if (c2 == 0.0)
	{
		x = 0;
		y = (c1 == 0)? 0 : -c3 / c1;
	}
	else
	{
		if (c1 == -c2)
		{
			x = INFINITY;
			y = INFINITY;
		}
		else
		{
			double m = -c1 / c2;
			double b = -c3 / c2;
			x = 1.0 / (1 - m);
			y = b * x;
		}
	}
}

void XPCPdataAnalyzer::calc_p_flats_fromTrajData(const vector<float>& p1, const vector<float>& p2, XPCPSample& xpcp_tuple)
{
	int N = p1.size(); // dimension of the input data

	// 0. Set N-D value 
	vector<float> pcpAxisVal = p1;
	xpcp_tuple.x = pcpAxisVal;
	// 1. Use the tangent as the line for 1-indexed point
	// Simply the difference between two neighboring points in time!
	vector<float> pStart = p1;
	vector<float> pEnd = p2;

	// convert to PCP
	vector<Point> pcp_p0_no_repeat;
	vector<Point> pcp_p1_no_repeat;
	_pcp->cartesian2PCP(pStart, pcp_p0_no_repeat);
	_pcp->cartesian2PCP(pEnd, pcp_p1_no_repeat);
	// Compute the 1-flat!
	vector<Point> pcp_v0_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);
	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
	{
		float x = pcp_v0_no_repeat[i].x - float(i); // org x-coordinate
		float y = pcp_v0_no_repeat[i].y; // org y-coordinate
		double u = x; // new x-coord
		double v = y; // new y-coord

		xform(x, y, i, u, v);
		/////////////////////////////////////////////////////////////////////////////////////////
		//===================================================================
		//**********************************************************************
		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(float(u) + float(i), float(v)/*pcp_v0_no_repeat[i].y*/);


		float vecLen = Eigen::Vector2f(p1[i], p1[i + 1]).norm();
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
		xpcp_tuple.strength_1flat[i].x = vecLen;
		_1flatsMinMaxPerSubspace[i].x = MIN(_1flatsMinMaxPerSubspace[i].x, vecLen);
		_1flatsMinMaxPerSubspace[i].y = MAX(_1flatsMinMaxPerSubspace[i].y, vecLen);
	}

	// 2. Compute 2-flat
	// We don't have 2-flats!
}

void XPCPdataAnalyzer::calc_p_flats_fromVecData(const vector<float>& vecVal, const vector<float>& pos, XPCPSample& xpcp_tuple)
{
	int N = vecVal.size(); // dimension of the input data

	// 0. Set N-D value by normalizing each channel
	vector<float> pcpAxisVal = vecVal;
	for (size_t i = 0; i < pcpAxisVal.size(); i++)
		pcpAxisVal[i] = (vecVal[i] - g_attrMinMax[i].x) / (g_attrMinMax[i].y - g_attrMinMax[i].x);
	xpcp_tuple.x = pcpAxisVal;
	// 1. Use the vector value as 1-flat!!!
	vector<float> pStart(N, 0.0f); // Use 0 as start point
	vector<float> pEnd = vecVal;   // Use vecVal as end point
	for (int i = 0; i < N; i++)
	{
		pStart[i] = g_randGen.rand() - 0.5f; // set a random starting point in [-0.5, 0.5] to have different "b" for the index point
		pEnd[i] = pStart[i] + vecVal[i];
	}
	// convert to PCP
	vector<Point> pcp_p0_no_repeat;
	vector<Point> pcp_p1_no_repeat;
	_pcp->cartesian2PCP(pStart, pcp_p0_no_repeat);
	_pcp->cartesian2PCP(pEnd, pcp_p1_no_repeat);
	// Compute the 1-flat!
	vector<Point> pcp_v0_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);
	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
	{
		float x = pcp_v0_no_repeat[i].x - float(i); // org x-coordinate
		float y = pcp_v0_no_repeat[i].y; // org y-coordinate
		double u = x; // new x-coord
		double v = y; // new y-coord

		xform(x, y, i, u, v);
		/////////////////////////////////////////////////////////////////////////////////////////
		//===================================================================
		//**********************************************************************
		xpcp_tuple.pcp_1flat[i] = FLOATVECTOR2(float(u) + float(i), float(v)/*pcp_v0_no_repeat[i].y*/);


		float vecLen = Eigen::Vector2f(vecVal[i], vecVal[i + 1]).norm();
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
		xpcp_tuple.strength_1flat[i].x = vecLen;
		_1flatsMinMaxPerSubspace[i].x = MIN(_1flatsMinMaxPerSubspace[i].x, vecLen);
		_1flatsMinMaxPerSubspace[i].y = MAX(_1flatsMinMaxPerSubspace[i].y, vecLen);
	}

	// 2. Compute 2-flat
	// We don't really have 2-flats for vector data!!!
	vector<FLOATVECTOR2> indexPt;
	bool oneIndex = true;

	//int normMethod = 1; // 0: length of the normal vector, 1: product of the lengths of two tangent vectors, 2: the larger of the tangents, 3: the smaller of the the tangents
	//for (int j = 0; j < N - 2; j++)
	//{
	//	// Get 3D subspace!!! Use projected eigenvectors as tangent vectors of the plane
	//	Eigen::Vector3f mu3d = mu.block(j, 0, 3, 1);
	//	Eigen::Vector3f majEig3d = majEig.block(j, 0, 3, 1);
	//	Eigen::Vector3f secEig3d = secEig.block(j, 0, 3, 1);
	//	//majEig3d.normalize(); // normalize these tangent vectors
	//	//secEig3d.normalize();

	//	// Do analytical computation to find 2-flat for the subspace
	//	Eigen::Vector3f norm3d = majEig3d.cross(secEig3d);
	//	float strength_2flat = 0.0f;
	//	if (norm3d == Eigen::Vector3f(0, 0, 0))
	//	{
	//		xpcp_tuple.strength_2flat[j].x = strength_2flat;
	//		xpcp_tuple.strength_2flat[j].y = strength_2flat;

	//		continue;
	//	}

	//	/////////////////////////////////////

	//	// Normalize the normal vector
	//	strength_2flat = computeNormalVecStrength(majEig3d, secEig3d, norm3d, normMethod);
	//	//strength_2flat = norm3d.norm();
	//	//////////////////////////////////////////

	//	_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
	//	_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);

	//	PlaneCoeffs Pl;
	//	Pl.c1 = norm3d.x();
	//	Pl.c2 = norm3d.y();
	//	Pl.c3 = norm3d.z();
	//	Pl.c0 = norm3d.x() * mu3d[0] + norm3d.y() * mu3d[1] + norm3d.z() * mu3d[2];

	//	// Compute only one index point
	//	_pcp->calc2FlatFrom3DPlane(Pl, indexPt, oneIndex);
	//	// Set 2-flat
	//	// Shift to the starting dimension
	//	indexPt[0].x += float(j);
	//	xpcp_tuple.pcp_2flat[j] = indexPt[0];
	//	xpcp_tuple.strength_2flat[j].x = strength_2flat; // set strength
	//}
}


void XPCPdataAnalyzer::calc1Flat(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, int ssId, FLOATVECTOR2& pFlat)
{
	Eigen::VectorXf ptMajEig = majEig + mu;
	vector<float> vPtMajEig(ptMajEig.data(), ptMajEig.data() + ptMajEig.rows());
	vector<float> vPtMu(mu.data(), mu.data() + mu.rows());
	// 1. Compute 1-flat
	// convert to PCP
	vector<Point> pcp_p0_no_repeat;
	vector<Point> pcp_p1_no_repeat;
	_pcp->cartesian2PCP(vPtMu, pcp_p0_no_repeat);
	_pcp->cartesian2PCP(vPtMajEig, pcp_p1_no_repeat);
	vector<Point> pcp_v0_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);

	// Convert to hyperbolic coordinates??
	// 1. https://en.wikipedia.org/wiki/Hyperbolic_coordinates
	if (pcp_v0_no_repeat[0].x < 0.0f) // x < 0 -> slope a > 1
		pcp_v0_no_repeat[0].x = -log(sqrt(-pcp_v0_no_repeat[0].x));
	else if (pcp_v0_no_repeat[0].x > 1.0f)
		pcp_v0_no_repeat[0].x = log(sqrt(pcp_v0_no_repeat[0].x));

	pFlat.x = pcp_v0_no_repeat[0].x + float(ssId); // offset by ssId 
	pFlat.y = pcp_v0_no_repeat[0].y;
}

void XPCPdataAnalyzer::calc2Flat(const Eigen::VectorXf& mu, const Eigen::VectorXf& majEig, const Eigen::VectorXf& secEig, int ssId, FLOATVECTOR2& pFlat)
{
	// 2. Compute 2-flat
	// For now, use 3-tuple in sequency of the axis ordering
	vector<FLOATVECTOR2> indexPt;
	bool oneIndex = true;
	int N = mu.size();
	pFlat = FLOATVECTOR2();

	int j = 0;
	// Get 3D subspace!!!
	Eigen::Vector3f mu3d = mu.block(j, 0, 3, 1);
	Eigen::Vector3f majEig3d = majEig.block(j, 0, 3, 1);
	Eigen::Vector3f secEig3d = secEig.block(j, 0, 3, 1);

	// Do analytical computation to find 2-flat for the subspace
	Eigen::Vector3f norm3d = majEig3d.cross(secEig3d);
	if (norm3d == Eigen::Vector3f(0, 0, 0))
		return; // if any of the eigen vectors are 0, quit
	PlaneCoeffs Pl;
	Pl.c1 = norm3d.x();
	Pl.c2 = norm3d.y();
	Pl.c3 = norm3d.z();
	Pl.c0 = norm3d.x() * mu3d[0] + norm3d.y() * mu3d[1] + norm3d.z() * mu3d[2];

	// Compute only one index point
	_pcp->calc2FlatFrom3DPlane(Pl, indexPt, oneIndex);
	// Set 2-flat
	// Shift to the starting dimension
	indexPt[0].x += float(ssId);
	pFlat = indexPt[0];

}

void XPCPdataAnalyzer::mapToRawData()
{}

bool vec_data_less(std::vector<float> v1, std::vector<float> v2) // For lower bound: return true if v1 < v2:
// We compare two vectors from the first element to the last
{
	bool lessThan = false;
	for (size_t i = 0; i < v1.size(); i++)
	{
		if (v1[i] > v2[i]){
			lessThan = false;
			break;
		}
		else if (v1[i] < v2[i]){
			lessThan = true;
			break;
		}
		else // Check next element
			continue;
	}
	return lessThan;
}

bool comp_pFlatPt(pFlatPt* v1, pFlatPt* v2)
{
	if (v1->subDimId() == v2->subDimId() &&
		v1->rawId() == v2->rawId() &&
		v1->x() == v2->x() &&
		v1->y() == v2->y())
		return true;
	else
		return false;
}

void XPCPdataAnalyzer::removeDuplicateRecords(const vector<vector<float>>& pointCenterData)
{
	// 1. Get unique records in the raw data
	set<vector<float>> uniqueRawDataSet(pointCenterData.begin(), pointCenterData.end());
	vector<vector<float>> uniqueRaw(uniqueRawDataSet.begin(), uniqueRawDataSet.end());
	std::sort(uniqueRaw.begin(), uniqueRaw.end(), vec_data_less); // sort 
	//ofstream ofUniqRaw("uniqueRaw.txt");
	//for (size_t i = 0; i < uniqueRaw.size(); i++)
	//{
	//	vector<float> raw = uniqueRaw[i];
	//	for (size_t j = 0; j < raw.size(); j++)
	//		ofUniqRaw << raw[j] << ",";
	//	ofUniqRaw << endl;
	//}
	//ofUniqRaw.close();

	// 2. Build the mapping from the original data to the unique data
	map<UINT64, UINT64> raw2uniqueMap;
	//ofstream ofR2Umap("raw2uniqueMap.txt");
	for (UINT64 i = 0; i < UINT64(pointCenterData.size()); i++)
	{
		// Locate all pointCenterData entries in uniqueRaw 
		vector<float> dataPoint = pointCenterData[i];
		vector<vector<float>>::iterator ITT = std::lower_bound(uniqueRaw.begin(), uniqueRaw.end(), dataPoint, vec_data_less); // find the first (only in our unique case) appearance of dataPoint in uniqueRaw
		if (ITT == uniqueRaw.end())
		{
			cout << "Something went wrong! Cannot find a value in pointCenterData in uniqueRaw" << endl;
			return;
		}
		raw2uniqueMap[i] = UINT64(ITT - uniqueRaw.begin());
		//ofR2Umap << i << "," << raw2uniqueMap[i] << endl;
	}
	//ofR2Umap.close();
	// 3. Re-assign unique rawIds to 1-flat array
	for (UINT64 i = 0; i < UINT64(g_1flat_list.size()); i++)
	{
		UINT64 newRawId = raw2uniqueMap[g_1flat_list[i]->rawId()];
		g_1flat_list[i]->setRawId(newRawId);
	}
	// 4. run unique on 1-flat list
	std::vector<pFlatPt*>::iterator it = std::unique(g_1flat_list.begin(), g_1flat_list.end(), comp_pFlatPt);
	g_1flat_list.erase(it, g_1flat_list.end());
	// 5. Repeat 3 for 2-flat
	for (UINT64 i = 0; i < UINT64(g_2flat_list.size()); i++)
	{
		UINT64 newRawId = raw2uniqueMap[g_2flat_list[i]->rawId()];
		g_2flat_list[i]->setRawId(newRawId);
	}
	// 4. run unique on 1-flat list
	g_2flat_list.erase(std::unique(g_2flat_list.begin(), g_2flat_list.end(), comp_pFlatPt), g_2flat_list.end());
	// save the unique, sorted raw data to a global variable!
	g_unique_sampl_raw = uniqueRaw;
}

float XPCPdataAnalyzer::computeNormalVecStrength(const Eigen::Vector3f& vT1, const Eigen::Vector3f& vT2, const Eigen::Vector3f& vNorm, int method)
{
	float g11 = 0.0f;
	float g12 = 0.0f;
	float g22 = 0.0f;
	float det = 1.0f;
	float len1 = 0.0f;
	float len2 = 0.0f;
	float avgLen = 1.0f;
	float maxLen = 1.0f;
	float area = 0.0f;
	float measure = 1.0f;
	Eigen::Matrix3f M;
	M(0, 0) = vNorm(0); M(0, 1) = vNorm(1); M(0, 2) = vNorm(2);
	M(1, 0) = vT1(0); M(1, 1) = vT1(1); M(1, 2) = vT1(2);
	M(2, 0) = vT2(0); M(2, 1) = vT2(1); M(2, 2) = vT2(2);
	static bool firstTime = true;
	if (firstTime){
		cout << M;
		firstTime = false;
	}
	switch (method)
	{
	case 0: // normalize by the length of the normal vector (considers the lengths of tangent vectors + the angle between them)
		measure = vNorm.norm();
		break;
	case 1: // normalize by the product of the lengths of two tangent vectors
		len1 = vT1.norm();
		len2 = vT2.norm();
		measure = len1 * len2;
		break;
	case 2: // normalize by the larger length of the tangent vectors
		maxLen = MAX(vT1.norm(), vT2.norm());
		measure = maxLen;
		break;
	case 3: // normalize by the smaller length of the tangents
		measure = MIN(vT1.norm(), vT2.norm()); //sqrtf(vT1.norm() * vT1.norm() + vT2.norm() * vT2.norm());
		break;
	default:
		measure = vNorm.norm();
		break;
	}
	return vNorm.norm();
}

void XPCPdataAnalyzer::nomralizePFlatsInSubspaces(vector<XPCPSample>& xpcpData, int normMethod)
{
	// For each subspace (2D for 1-flat, 3D for 2-flat), normalize the strength of p-flats
	// Test: output min/max of p-flat strength in each subspace
	for (size_t i = 0; i < _1flatsMinMaxPerSubspace.size(); i++)
		cout << "1-flat subspace " << i << ", min=" << _1flatsMinMaxPerSubspace[i].x << ", max=" << _1flatsMinMaxPerSubspace[i].y << endl;
	for (size_t i = 0; i < _2flatsMinMaxPerSubspace.size(); i++)
		cout << "2-flat subspace " << i << ", min=" << _2flatsMinMaxPerSubspace[i].x << ", max=" << _2flatsMinMaxPerSubspace[i].y << endl;
	// Actual normalization by population!!!
	if (normMethod == 0)
	{
		vector<FLOATVECTOR3> pdf1flats, pdf2flats, cdf1flats, cdf2flats;
		int mode = 0; // 1flat
		int numBin = MIN(256, xpcpData.size()); // number of bins
		for (int ssid = 0; ssid < _dimRawData - 1; ssid++) // Normalize by subspace
		{
			mode = 0; // 1flat
			// Compute pdf
			computeSubspacePDF(pdf1flats, ssid, mode, numBin, xpcpData);
			// compute cdf
			computeCDF(cdf1flats, pdf1flats);
			// find & set strength for p-flat using percentage
			findSetSubspaceSamplePercentage(xpcpData, ssid, mode, cdf1flats);
			//if (ssid < _dimRawData - 2)
			//{
			//	mode = 1;
			//	// Compute pdf
			//	computeSubspacePDF(pdf2flats, ssid, mode, numBin, xpcpData);
			//	// compute cdf
			//	computeCDF(cdf2flats, pdf2flats);
			//	// find & set strength for p-flat using percentage
			//	findSetSubspaceSamplePercentage(xpcpData, ssid, mode, cdf2flats);
			//}
		}
	}

}

void XPCPdataAnalyzer::computeSubspacePDF(vector<FLOATVECTOR3>& pdf, int subspaceId, int mode, int numBin, const vector<XPCPSample>& xpcpData)
{
	assert(mode >= 0 && mode < 2); // make sure mode is either 0 or 1!!!
	// prepare output
	pdf.resize(numBin);
	for (size_t i = 0; i < pdf.size(); i++)
	{
		pdf[i].x = numeric_limits<float>::max();
		pdf[i].y = -numeric_limits<float>::max();
		pdf[i].z = 0.0f;
	}
	// check mode
	if (mode == 0) // 1-flat
	{
		float oneOverRange = 1.0f / (_1flatsMinMaxPerSubspace[subspaceId].y - _1flatsMinMaxPerSubspace[subspaceId].x);
		if (_1flatsMinMaxPerSubspace[subspaceId].y == _1flatsMinMaxPerSubspace[subspaceId].x)
			oneOverRange = 0.0f;
		// TODO: what happens if min == max?? Should handle differently in the caller function
		if (_1flatsMinMaxPerSubspace[subspaceId].y == _1flatsMinMaxPerSubspace[subspaceId].x)
			oneOverRange = 0.0f;
		for (size_t i = 0; i < xpcpData.size(); i++)
		{
			float val = xpcpData[i].strength_1flat[subspaceId].x;
			float normVal = (val - _1flatsMinMaxPerSubspace[subspaceId].x) * oneOverRange;
			int bin = normVal * (numBin - 1);
			pdf[bin].x = MIN(val, pdf[bin].x);
			pdf[bin].y = MAX(val, pdf[bin].y);
			pdf[bin].z++;
		}

	}
	else // 2-flat
	{
		float oneOverRange = 1.0f / (_2flatsMinMaxPerSubspace[subspaceId].y - _2flatsMinMaxPerSubspace[subspaceId].x);
		// TODO: what happens if min == max?? Should handle differently in the caller function
		if (_2flatsMinMaxPerSubspace[subspaceId].y == _2flatsMinMaxPerSubspace[subspaceId].x)
			oneOverRange = 0.0f;
		for (size_t i = 0; i < xpcpData.size(); i++)
		{
			float val = xpcpData[i].strength_2flat[subspaceId].x;
			if (val == 0.0f)
				continue;
			float normVal = (val - _2flatsMinMaxPerSubspace[subspaceId].x) * oneOverRange;
			int bin = normVal * (numBin - 1);
			pdf[bin].x = MIN(val, pdf[bin].x);
			pdf[bin].y = MAX(val, pdf[bin].y);
			pdf[bin].z++;
		}
	}
}

void XPCPdataAnalyzer::computeCDF(vector<FLOATVECTOR3>& cdf, const vector<FLOATVECTOR3>& pdf)
{
	cdf.resize(pdf.size());
	for (size_t i = 0; i < pdf.size(); i++)
	{
		cdf[i].x = pdf[i].x;
		cdf[i].y = pdf[i].y;
		if (i == 0)
			cdf[i].z = pdf[i].z;
		else
			cdf[i].z = pdf[i].z + cdf[i - 1].z; // accumulate pdf
	}
	for (size_t i = 0; i < cdf.size(); i++)
		cdf[i].z /= float(_numSamples); // normalize cdf
}

void XPCPdataAnalyzer::findSetSubspaceSamplePercentage(std::vector<XPCPSample>& xpcpData, int subspaceId, int mode, std::vector<FLOATVECTOR3>& cdf)
{
	assert(mode >= 0 && mode < 2); // make sure mode is either 0 or 1!!!
	// check mode
	if (mode == 0) // 1-flat
	{
		float oneOverRange = 1.0f / (_1flatsMinMaxPerSubspace[subspaceId].y - _1flatsMinMaxPerSubspace[subspaceId].x);
		// TODO: what happens if min == max?? Should handle differently in the caller function
		if (_1flatsMinMaxPerSubspace[subspaceId].y == _1flatsMinMaxPerSubspace[subspaceId].x)
			oneOverRange = 0.0f;

		{
			for (size_t i = 0; i < xpcpData.size(); i++)
			{
				float val = xpcpData[i].strength_1flat[subspaceId].x;
				float normVal = (val - _1flatsMinMaxPerSubspace[subspaceId].x) * oneOverRange;
				int bin = normVal * (cdf.size() - 1);
				// look up the percentage
				float percentage = cdf[bin].z;
				// set 1-flat strength using this percentage
				xpcpData[i].strength_1flat[subspaceId].y = percentage;
			}
		}


	}
	else // 2-flat
	{
		float oneOverRange = 1.0f / (_2flatsMinMaxPerSubspace[subspaceId].y - _2flatsMinMaxPerSubspace[subspaceId].x);
		// TODO: what happens if min == max?? Should handle differently in the caller function
		if (_2flatsMinMaxPerSubspace[subspaceId].y == _2flatsMinMaxPerSubspace[subspaceId].x)
			oneOverRange = 0.0f;
		for (size_t i = 0; i < xpcpData.size(); i++)
		{
			float val = xpcpData[i].strength_2flat[subspaceId].x;
			float normVal = (val - _2flatsMinMaxPerSubspace[subspaceId].x) * oneOverRange;
			int bin = normVal * (cdf.size() - 1);
			// look up the percentage
			float percentage = cdf[bin].z;
			// set 2-flat percentile 
			xpcpData[i].strength_2flat[subspaceId].y = percentage;
		}
	}
}

void XPCPdataAnalyzer::uncertain_getPointCenters(const vector<GaussianXd>& distrData, vector<vector<float>>& pointCenterData)
{
	pointCenterData.resize(distrData.size());
	// Get the mean value of each Gaussian 
	for (vector<GaussianXd>::const_iterator IT = distrData.begin(); IT != distrData.end(); ++IT)
	{
		Eigen::VectorXd mu = IT->m_mu;
		vector<float>   pointCenter(mu.size());
		for (size_t i = 0; i < mu.size(); i++)
			pointCenter[i] = mu(i);
		pointCenterData[IT - distrData.begin()] = pointCenter;
	}
}

void XPCPdataAnalyzer::uncertain_sampleNeighborhood(const vector<GaussianXd>& distrData, int* nnIdx, int num_nn, vector<Eigen::MatrixXd>& pointNeighbors, int num_samples_per_distr)
{
	// Monte-Carlo sampling for the neighborhood
	if (pointNeighbors.size() != size_t(num_nn))
		pointNeighbors.resize(num_nn);
	for (int i = 0; i < num_nn; i++)
	{
		// Get the corresponding neighbor distribution
		int id = nnIdx[i];
		GaussianXd pointDistr = distrData[id];
		Eigen::MatrixXd point_Samples;
		// Sample the given multivariate Gaussian distribution
		uncertain_drawSamplesPointDistrib(pointDistr, num_samples_per_distr, point_Samples);
		// Save the result
		pointNeighbors[i] = point_Samples;
	}
}

void XPCPdataAnalyzer::uncertain_drawSamplesPointDistrib(const GaussianXd& pointDistr, int num_samples, Eigen::MatrixXd& samples)
{
	// the multivariate samples in the eigen::MatrixXd form
	drawSamplesFromGaussianXd(pointDistr, num_samples, samples);
}

void XPCPdataAnalyzer::uncertain_get_neighbors_one_sample(const vector<Eigen::MatrixXd>& sampled_neighbors, int sample_id, vector<Eigen::VectorXd>& neighbors_one_sample)
{
	if (neighbors_one_sample.size() != sampled_neighbors.size())
		neighbors_one_sample.resize(sampled_neighbors.size());
	for (vector<Eigen::MatrixXd>::const_iterator IT = sampled_neighbors.begin(); IT != sampled_neighbors.end(); ++IT)
	{
		Eigen::MatrixXd samples_one_distr = *IT;
		// Get a random combination of the samples!

		Eigen::VectorXd one_rec(samples_one_distr.rows());
		for (int ri = 0; ri < one_rec.rows(); ri++)
		{
			one_rec(ri) = samples_one_distr(ri, g_randGen.randInt() % samples_one_distr.cols());
		}

		neighbors_one_sample[IT - sampled_neighbors.begin()] = one_rec;
	}
}

void XPCPdataAnalyzer::uncertain_get_neighbors_from_all_samples(const vector<Eigen::MatrixXd>& all_samples, vector<Eigen::VectorXd>& neighbors)
{
	if (neighbors.size() != all_samples.size())
		neighbors.resize(all_samples.size());
	for (size_t i = 0; i < neighbors.size(); i++)
	{
		int ri = g_randGen.randInt() % (all_samples.size() - 1);
		Eigen::MatrixXd samples = all_samples[ri];
		Eigen::VectorXd one_rec(all_samples[0].rows());
		for (int rr = 0; rr < one_rec.rows(); rr++){
			one_rec(rr) = samples(rr, g_randGen.randInt() % (samples.cols() - 1));
		}
		neighbors[i] = one_rec;
	}

}
void XPCPdataAnalyzer::uncertain_calc_p_flats_subSpaceMethod(const Eigen::VectorXd& mu, const Eigen::VectorXd& majEig, const Eigen::VectorXd& secEig, UncertainXPCPSample& xpcp_tuple, int sampleId)
{
	int N = mu.size(); // dimension of the input data

	Eigen::VectorXd ptMajEig = majEig + mu;
	vector<float> vPtMajEig(ptMajEig.size());
	vector<float> vPtMu(ptMajEig.size());
	for (size_t i = 0; i < vPtMajEig.size(); i++)
	{
		vPtMajEig[i] = float(ptMajEig(i));
		vPtMu[i] = float(mu(i));
	}
	// 0. Set N-D value
	if (sampleId == 0) // set x for the first sample only
		xpcp_tuple.x = vPtMu;

	// 1. Compute 1-flat
	// convert to PCP
	vector<Point> pcp_p0_no_repeat;
	vector<Point> pcp_p1_no_repeat;
	_pcp->cartesian2PCP(vPtMu, pcp_p0_no_repeat);
	// threshold the length of the major eigen vector
	float twoFlat_discard_thres = 1e-3f;
	float oneFlat_discard_thres = 1e-3f;

	_pcp->cartesian2PCP(vPtMajEig, pcp_p1_no_repeat);
	vector<Point> pcp_v0_no_repeat;
	_pcp->calcFlatOneOrderHigher(pcp_p0_no_repeat, pcp_p1_no_repeat, pcp_v0_no_repeat);
	//if ( sampleId == 0 )
	//	xpcp_tuple.pcp_1flat.resize(pcp_v0_no_repeat.size());
	for (size_t i = 0; i < pcp_v0_no_repeat.size(); i++)
	{
		float x = pcp_v0_no_repeat[i].x - float(i); // org x-coordinate
		float y = pcp_v0_no_repeat[i].y; // org y-coordinate
		double u = x; // new x-coord
		double v = y; // new y-coord

		if (g_params.XformMethod() == XF_IDX_ONLY || g_params.XformMethod() == XF_IDX_LINE)
			xform(x, y, i, u, v);
		/////////////////////////////////////////////////////////////////////////////////////////
		//===================================================================
		//**********************************************************************
		int pflatid = sampleId * pcp_v0_no_repeat.size() + i;
		xpcp_tuple.pcp_1flat[pflatid] = FLOATVECTOR2(float(u) + float(i), float(v));

		float eigVLen = Eigen::Vector2f(majEig[i], majEig[i + 1]).norm();
		///////////////////////////////////////////////////////////
		// Record strength of the eigen vector in each subspace and get the min/max
		xpcp_tuple.strength_1flat[i].x = eigVLen;
		//**********************************************************
		// Alternatively, use the correlation value for 1-flat strength
		// NOTE: Uncomment to use correlation value as strength!
		//	xpcp_tuple.strength_1flat[i].x = abs(corr_mat(i,i+1));
		//**********************************************************
		float minEigVlen = MIN(_1flatsMinMaxPerSubspace[i].x, xpcp_tuple.strength_1flat[i].x);
		_1flatsMinMaxPerSubspace[i].x = minEigVlen;
		float maxEigVlen = MAX(_1flatsMinMaxPerSubspace[i].y, xpcp_tuple.strength_1flat[i].x);
		_1flatsMinMaxPerSubspace[i].y = maxEigVlen;
	}

	// 2. Compute 2-flat
	// For now, use 3-tuple in sequency of the axis ordering
	vector<FLOATVECTOR2> indexPt;
	bool oneIndex = !g_use_repeat_pcp;
	int numIdxPerPlane = oneIndex ? 1 : 4;

	if (xpcp_tuple.pcp_2flat.size() != xpcp_tuple.num_pts * (N - 2) * numIdxPerPlane)
		xpcp_tuple.pcp_2flat.resize(xpcp_tuple.num_pts * (N - 2) * numIdxPerPlane);
	int normMethod = 1; // 0: length of the normal vector, 1: product of the lengths of two tangent vectors, 2: the larger of the tangents, 3: the smaller of the the tangents
	for (int j = 0; j < N - 2; j++)
	{
		// Get 3D subspace!!! Use projected eigenvectors as tangent vectors of the plane
		Eigen::Vector3f mu3d = mu.cast<float>().block(j, 0, 3, 1);
		Eigen::Vector3f majEig3d = majEig.cast<float>().block(j, 0, 3, 1);
		Eigen::Vector3f secEig3d = secEig.cast<float>().block(j, 0, 3, 1);

		// Do analytical computation to find 2-flat for the subspace
		Eigen::Vector3f norm3d = majEig3d.cross(secEig3d);
		float strength_2flat = 0.0f;
		if (norm3d == Eigen::Vector3f(0, 0, 0))
		{
			xpcp_tuple.strength_2flat[j].x = strength_2flat;
			xpcp_tuple.strength_2flat[j].y = strength_2flat;

			_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
			_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);
			continue;
		}

		/////////////////////////////////////

		// Normalize the normal vector
		strength_2flat = computeNormalVecStrength(majEig3d, secEig3d, norm3d, normMethod);
		//strength_2flat = norm3d.norm();
		//////////////////////////////////////////

		_2flatsMinMaxPerSubspace[j].x = MIN(strength_2flat, _2flatsMinMaxPerSubspace[j].x);
		_2flatsMinMaxPerSubspace[j].y = MAX(strength_2flat, _2flatsMinMaxPerSubspace[j].y);

		PlaneCoeffs Pl;
		// The plane needs a normalized normal vector
		norm3d.normalize();
		Pl.c1 = norm3d.x();
		Pl.c2 = norm3d.y();
		Pl.c3 = norm3d.z();
		Pl.c0 = norm3d.x() * mu3d[0] + norm3d.y() * mu3d[1] + norm3d.z() * mu3d[2];

		// Compute indexed point(s)
		_pcp->calc2FlatFrom3DPlane(Pl, indexPt, oneIndex);
		// Set 2-flat
		int pflatid_base = sampleId * numIdxPerPlane * (N - 2);
		for (int k = 0; k < numIdxPerPlane; k++)
		{
			// Shift to the starting dimension
			indexPt[k].x += float(j);
			xpcp_tuple.pcp_2flat[pflatid_base + numIdxPerPlane * j + k] = indexPt[k];
			xpcp_tuple.strength_2flat[j].x = strength_2flat; // set strength
		}

	}
}

vector<UINT64> XPCPdataAnalyzer::detectOutliers(const std::vector<std::vector<float>>& rawData)
{
	
	vector<UINT64> outlierIdList;
	if (rawData.empty())
	{
		cout << "Input data is empty! Outlier detection skipped!" << endl;
		return outlierIdList;
	}

	if (rawData.size() < 100)
	{
		cout << "Input data is too few! Outlier detection skipped!" << endl;
		return outlierIdList;
	}
	// n: # of rows, d: # of columns
	UINT64 n = rawData.size();
	UINT64 d = rawData[0].size();
	UINT64 n_sample = MAX(20, UINT64(0.01 * double(n))); // try 10% of the sampled raw? 0.1 * n;
	UINT64 seed = 0; // skip seed, use random
	UINT64 k = 70; // keep only 50 outliers
	// compute qsps
	double  *X = new double[n * d];
	for (size_t i = 0; i < rawData.size(); i++)
	{
		for (size_t dd = 0; dd < d; dd++)
			X[i * d + dd] = double(rawData[i][dd]);
	}
	printf("computing qsp scores...");
	fflush(stdout);
	double* score = new double[n];
	qsp(X, n, d, n_sample, seed, score);
	printf("end\n");
	fflush(stdout);


	

	// index qsort for ranking
	UINT64* index = NULL;
	bool srt = true; // do sort
	if (srt) {
		index = new UINT64[n];
		outlierIdList.resize(k);
		for (UINT64 i = 0; i < n; i++) {
			index[i] = i;
		}
		iqsort(index, score, n);

		printf("\n");
		printf("indexes of top-%d ourliers:\n", k);
		printf("  ");
		for (UINT64 i = 0; i < k/*(k - 1)*/; i++) {
			printf("%ld, ", index[i] + 1);
			outlierIdList[i] = index[i];
		}
		//printf("%ld\n", index[i] + 1);
	}
	// write resulting scores to an output file
	string output = "outlier_scores.txt";
	printf("writing scores to the file \"%s\"...", output.c_str());
	fflush(stdout);
	FILE* fp = fopen(output.c_str(), "w");
	if (fp == NULL) {
		printf("%s: cannot open the file\n", output);
		exit(1);
	}
	for (UINT64 i = 0; i < n; i++) {
		fprintf(fp, "%f\n", score[i]);
	}
	fclose(fp);
	printf("end\n");
	fflush(stdout);

	delete[] X;
	X = NULL;
	delete[] score;
	score = NULL;
	delete[] index;
	index = NULL;
	return outlierIdList;
}

void XPCPdataAnalyzer::testEigenVectorsInterp()
{
	vector<UINT64VECTOR3> offset;
	offset.push_back(UINT64VECTOR3(1, 0, 0));
	offset.push_back(UINT64VECTOR3(-1, 0, 0));
	offset.push_back(UINT64VECTOR3(0, 1, 0));
	offset.push_back(UINT64VECTOR3(0, -1, 0));
	offset.push_back(UINT64VECTOR3(0, 0, 1));
	offset.push_back(UINT64VECTOR3(0, 0, -1));


	vector<FLOATVECTOR3> eigVList;
	vector<float> eigvalList;
	vector<float> tList;

	UINT64VECTOR3 testPos = UINT64VECTOR3(10, 10, 10);
	for (size_t m = 0; m < g_params.majEigVols().size(); m++)
	{
		VolumeData* vol = g_params.majEigVols()[m];


		for (size_t nn = 0; nn < offset.size(); nn++)
		{
			UINT64VECTOR3 pos = testPos + offset[nn];
			if (pos.x < 0 || pos.y < 0 || pos.z < 0 ||
				pos.x >= vol->getDim().x || pos.y >= vol->getDim().y || pos.z >= vol->getDim().z)
			{
				float vx = vol->getVoxel(pos, 0);
				float vy = vol->getVoxel(pos, 1);
				float vz = vol->getVoxel(pos, 2);
				FLOATVECTOR3 eigVorg = FLOATVECTOR3(vx, vy, vz);
				float eigVal = eigVorg.length();
				eigVList.push_back(eigVorg);
				eigvalList.push_back(eigVal);
				tList.push_back(0.5f);
			}

		}

		FLOATVECTOR3 vt;
		float lamt;
		eigenVecInterpND(vt, eigVList, tList);
		eigenValInterpND(lamt, eigvalList, tList);
	}
}

void XPCPdataAnalyzer::cluster(const std::vector<std::vector<float>>& rawData, std::vector<int>& labelData, int numK)
{
	if (rawData.empty())
		return;
	int rawDataDim = rawData[0].size();
	int rawDataCount = rawData.size();
	cv::Mat matInput(rawDataCount, rawDataDim, CV_32FC1);
	cv::Mat matLabel;

	int attempts = 8;
	//cv::Mat centers;
	//cv::kmeans(matInput, numK, matLabel, cv::TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1), attempts, cv::KMEANS_RANDOM_CENTERS);
	cv::kmeans(matInput, numK, matLabel, cv::TermCriteria(cv::TermCriteria::COUNT, 10, 1), attempts, cv::KMEANS_RANDOM_CENTERS);

	// set label data
	labelData.resize(rawDataCount, 0);
	for (int i = 0; i < rawDataCount; i++)
		labelData[i] = matLabel.at<int>(i);
}

void XPCPdataAnalyzer::eigenVecInterp1D(FLOATVECTOR3& vt, const FLOATVECTOR3& v1, const FLOATVECTOR3& v2, float t)
{
	float thetCos1 = v1 ^ v2;
	float theta1 = acosf(thetCos1);
	vt = sin((1.0f - t) * theta1) / sin(theta1) * v1 + sin(t * theta1) / sin(theta1) * v2;
}

void XPCPdataAnalyzer::eigenVecInterp1D(Eigen::VectorXf& vt, const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float t)
{
	float thetaCos = v1.dot(v2);
	float theta = acosf(thetaCos);
	vt = sin((1.0f - t) * theta) / sin(theta) * v1 + sin(t * theta) / sin(theta) * v2;
}

void XPCPdataAnalyzer::eigenVecInterpND(FLOATVECTOR3& vt, const std::vector<FLOATVECTOR3>& v, const std::vector<float>& t)
{
	//U\left(t\right)=\left({\nu }_{1}\left(t\right),{\nu }_{2}\left(t\right),{\nu }_{3}\left(t\right)\right)=\sum_{i=1}^{4}{U}_{i}diags\left(\frac{sin{t}_{i}{\theta }_{1}}{sin{\theta }_{1}},\frac{sin{t}_{i}{\theta }_{2}}{sin{\theta }_{2}},\frac{sin{t}_{i}{\theta }_{3}}{sin{\theta }_{3}}\right)
	// For the case of 2D, we have 4 items of v
	int n = v.size();
	float theta = -FLT_MAX;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			float thetai = acosf(v[i] ^ v[j]);
			theta = MAX(thetai, theta);
		}
	}

	vt = FLOATVECTOR3(0,0,0);
	for (int i = 0; i < n; i++)
	{
		float w = sin(t[i] * theta) / sin(theta);
		vt += w * v[i];
	}
}

void XPCPdataAnalyzer::eigenVecInterpND(Eigen::VectorXf& vt, const std::vector<Eigen::VectorXf>& v, const std::vector<float>& t)
{
	int n = v.size();
	float theta = -FLT_MAX;
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			float thetai = acosf(v[i].dot(v[j]));
			theta = MAX(thetai, theta);
		}
	}

	vt = Eigen::VectorXf::Zero(3);
	for (int i = 0; i < n; i++)
	{
		float w = sin(t[i] * theta) / sin(theta);
		vt += w * v[i];
	}
}

void XPCPdataAnalyzer::eigenValInterpND(float& lamt, const std::vector<float>& lambda, const std::vector<float>& t)
{
	//\Lambda \left(t\right)=exp\left(\sum_{i=1}^{4}{t}_{i}log{\Lambda }_{i}\right)
	lamt = 0.0f;
	for (size_t i = 0; i < lambda.size(); i++)
	{
		lamt += t[i] * logf(lambda[i]);
	}
	lamt = expf(lamt);
}

void XPCPdataAnalyzer::matEigInterp1D(std::vector<Eigen::VectorXf>& Mt, Eigen::VectorXf& Lambt, const Eigen::MatrixXf& M1, const Eigen::MatrixXf& M2, const Eigen::VectorXf& Lamb1, const Eigen::VectorXf& Lamb2, float t)
{
	int d = 3;
	Mt.clear();
	vector<float> ts;
	ts.push_back(1.0f - t);
	ts.push_back(t);
	for (int i = 0; i < d; i++)
	{
		vector<float> lambs;
		lambs.push_back(Lamb1[i]);
		lambs.push_back(Lamb2[i]);

		// Eigen value interpolation
		eigenValInterpND(Lambt[i], lambs, ts);
		// Eigen vector interpolation
		//majEigV = eigVec.col(eigVec.cols() - 1);
		Eigen::VectorXf v = M1.col(i);
		eigenVecInterp1D(v, M1.col(i), M2.col(i), t);
		Mt.push_back(v);
	}
}

void XPCPdataAnalyzer::matEigInterpND(std::vector<Eigen::VectorXf>& Mt, Eigen::VectorXf& Lambt, const std::vector<Eigen::MatrixXf>& Ms, const std::vector<Eigen::VectorXf>& Lambs, const std::vector<float>& t)
{
	int d = 3;
	Mt.clear();
	int n = Lambs.size();
	for (int i = 0; i < d; i++)
	{
		vector<float> lambda;
		Eigen::VectorXf v = Eigen::VectorXf::Zero(3);
		vector<Eigen::VectorXf> vs;
		for (int j = 0; j < n; j++) {
			lambda.push_back(Lambs[j][i]);
			vs.push_back(Ms[j].col(i));
		}
		// Eigen value interpolation
		eigenValInterpND(Lambt[i], lambda, t);
		// Eigen vector interpolation
		//majEigV = eigVec.col(eigVec.cols() - 1);
		eigenVecInterpND(v, vs, t);
		Mt.push_back(v);
	}
}


