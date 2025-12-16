#include "PointMonteCarloSampler.h"
#include "ScatterPlotBuilder.h"
#include "PCPbuilder.h"
#include "ScatterPlot3DBuilder.h"
#include "SplomPlotBuilder.h"
#include "MyStatistics.h"
#include "Eigen/SVD"
#include "ANN/ANN.h"
using namespace std;
// Seeding parameters

PointMonteCarloSampler::PointMonteCarloSampler(const vector<VolumeData*>& volData, UINT64VECTOR3 volDim, string& strOutputFile)
:DirectionalMonteCarloSampler(volData, volDim, strOutputFile)
{
	_sample_portion = 0.0f;
	_neighborSize = 0;
	_numTotalSamples = 0;
	_sampleCnt = 0;

	//  Test 1: verify mean, covariance computation
	vector<vector<float>> testSamples;
	testSamples.resize(5);
	testSamples[0].push_back(90);	testSamples[0].push_back(60); 	testSamples[0].push_back(90);
	testSamples[1].push_back(90);	testSamples[1].push_back(90);	testSamples[1].push_back(30);
	testSamples[2].push_back(60);	testSamples[2].push_back(60);	testSamples[2].push_back(60);
	testSamples[3].push_back(60);	testSamples[3].push_back(60);	testSamples[3].push_back(90);
	testSamples[4].push_back(30);	testSamples[4].push_back(30);	testSamples[4].push_back(30);
	Eigen::VectorXd mu;
	MyStatistics::ComputeMean<float>(testSamples, mu);
	cout << "test mean=" << endl << mu << endl;
	Eigen::MatrixXd cov;
	MyStatistics::ComputeCovariance<float>(testSamples, mu, cov);
	cout << "test cov= " << endl << cov << endl;

	// Test 2: SVD
	Eigen::MatrixXd A(3, 3);
	A << 1.5, -0.5, 0.0,
		-0.5, 1.5, 0.0,
		0.0, 0.0, 3.0;
	cout << "mat A for SVD ="<<endl<<A<<endl;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	cout << "U = " << endl << svd.matrixU() << endl;
	cout << "D = " << endl << svd.singularValues() << endl;
	cout << "V = " << endl << svd.matrixV() << endl;
}

PointMonteCarloSampler::~PointMonteCarloSampler()
{}

void PointMonteCarloSampler::sample(float portion_voxels, int neighborSize)
{
	if (!init(portion_voxels, neighborSize))
		return; 

	for (UINT64 i = 0; i < _numTotalSamples; i++)
	{
		_sampleCnt = i;
		//// Sample the integer location with neighborhood
		//UINT64VECTOR3 loc = UINT64VECTOR3(_rng->randInt(UINT32(_volDim.x)), _rng->randInt(UINT32(_volDim.y)), _rng->randInt(UINT32(_volDim.z)));
		//sampleCubeNeighborhood(loc, neighborSize);

		// Sample the normalized continuous location with neighborhood
		FLOATVECTOR3 loc = FLOATVECTOR3(_rng->rand(), _rng->rand(), _rng->rand());
		// Rejection sampling in the scatterplot
		// Compare to a ratio of the dense of the scatterplot


		// Sample only on every 10 z slices
		//loc.z = float(_rng->randInt(9) * 10) / float(_volDim.z);
		sampleCubeNeighborhoodNormCoords(loc, neighborSize);
	}
}

vector<vector<float>> PointMonteCarloSampler::buildValueDomain(float portion_voxels, vector<vector<float>>& sampledData)
{
	vector<vector<float>> sampleLocs;
	cout << "Build value domain!" << endl;
	if (!init(portion_voxels, 1))
	{
		cout << "Failed to initialize sampling parameters!" << endl;
		return sampleLocs;
	}

	sampleLocs.resize(_numTotalSamples);
	sampledData.resize(_numTotalSamples);
	for (UINT64 i = 0; i < _numTotalSamples; i++)
	{
		// Sample the normalized continuous location
		FLOATVECTOR3 loc = FLOATVECTOR3(_rng->rand(), _rng->rand(), _rng->rand());
		//sample
		sampledData[i] = sampleVoxelWithInterpNormCoords(loc);
		// record sample location
		sampleLocs[i].resize(3);
		sampleLocs[i][0] = loc.x;
		sampleLocs[i][1] = loc.y;
		sampleLocs[i][2] = loc.z;
	}
	return sampleLocs;
}

std::vector<std::vector<float>> PointMonteCarloSampler::sampleSpatialDomain(float portion_voxels, std::vector<std::vector<float>>& sampledData)
{
	cout << "Sample the spatial domain!" << endl;
	return buildValueDomain(portion_voxels, sampledData);
}

vector<vector<float>> PointMonteCarloSampler::buildValueDomain_regDownsample_volDimSamplePos(FLOATVECTOR3 samplePortion, vector<vector<float>>& sampledData)
{
	INTVECTOR3 sampleDim = INTVECTOR3(samplePortion.x * _volDim.x, samplePortion.y * _volDim.y, MAX(1, samplePortion.z * _volDim.z));
	return buildValueDomain_regDownsample_volDimSamplePos(sampleDim, sampledData);
}

vector<vector<float>> PointMonteCarloSampler::buildValueDomain_regDownsample_volDimSamplePos(INTVECTOR3 sampleDim, vector<vector<float>>& sampledData)
{
	vector<vector<float>> sampleLocs;
	cout << "Build value domain by downsampling of a target size =" << sampleDim << "!"<<endl;
	if (!init(sampleDim, 1))
	{
		cout << "Failed to initialize sampling parameters!" << endl;
		return sampleLocs;
	}

	sampleLocs.resize(_numTotalSamples);
	sampledData.resize(_numTotalSamples);

	INTVECTOR3 inc = INTVECTOR3(_volDim.x / sampleDim.x, _volDim.y / sampleDim.y, _volDim.z / sampleDim.z);
	UINT64 cnt = 0;
	vector<float> val(_volData.size());
	for (int z = int(inc.z/2); z < _volDim.z; z+=inc.z)
	{
		for (int y = int(inc.y/2); y < _volDim.y; y+=inc.y)
		{
			for (int x = int(inc.x/2); x < _volDim.x; x+=inc.x)
			{
				FLOATVECTOR3 loc = FLOATVECTOR3(x, y, z);
				
				for (int i = 0; i < _volData.size(); i++)
					val[i] = _volData[i]->getVoxel(UINT64VECTOR3(x, y, z)) / _volValNormalizer;
				if (cnt >= _numTotalSamples)
					break;
				sampledData[cnt]= val;

				sampleLocs[cnt].resize(3);
				sampleLocs[cnt][0] = loc.x;
				sampleLocs[cnt][1] = loc.y;
				sampleLocs[cnt][2] = loc.z;
				cnt++;
			}
		}
	}
	return sampleLocs;
}

vector<vector<float>> PointMonteCarloSampler::buildVecValueDomain_regDownsample_volDimSamplePos(FLOATVECTOR3 samplePortion, vector<vector<float>>& sampledData)
{
	INTVECTOR3 sampleDim = INTVECTOR3(samplePortion.x * _volDim.x, samplePortion.y * _volDim.y, MAX(1, samplePortion.z * _volDim.z));
	return buildVecValueDomain_regDownsample_volDimSamplePos(sampleDim, sampledData);
}


vector<vector<float>> PointMonteCarloSampler::buildVecValueDomain_regDownsample_volDimSamplePos(INTVECTOR3 sampleDim, vector<vector<float>>& sampledData)
{
	vector<vector<float>> sampleLocs;
	cout << "Build value domain by downsampling of a target size =" << sampleDim << "!" << endl;
	if (!init(sampleDim, 1))
	{
		cout << "Failed to initialize sampling parameters!" << endl;
		return sampleLocs;
	}

	sampleLocs.resize(_numTotalSamples);
	sampledData.resize(_numTotalSamples);

	INTVECTOR3 inc = INTVECTOR3(_volDim.x / sampleDim.x, _volDim.y / sampleDim.y, _volDim.z / sampleDim.z);
	UINT64 cnt = 0;
	int totalValDim = 0;
	for (int i = 0; i < _volData.size(); i++)
		totalValDim += _volData[i]->numChannels();
	vector<float> val(totalValDim);

	for (int z = int(inc.z / 2); z < _volDim.z; z += inc.z)
	{
		for (int y = int(inc.y / 2); y < _volDim.y; y += inc.y)
		{
			for (int x = int(inc.x / 2); x < _volDim.x; x += inc.x)
			{
				FLOATVECTOR3 loc = FLOATVECTOR3(x, y, z);

				int dv = 0;
				for (int i = 0; i < _volData.size(); i++) {
					int nch = _volData[i]->numChannels();
					for (int ch = 0; ch < nch; ch++)
						val[dv++] = _volData[i]->getVoxel(UINT64VECTOR3(x, y, z),ch) / _volValNormalizer;

				}
				if (cnt >= _numTotalSamples)
					break;
				sampledData[cnt] = val;

				sampleLocs[cnt].resize(3);
				sampleLocs[cnt][0] = loc.x;
				sampleLocs[cnt][1] = loc.y;
				sampleLocs[cnt][2] = loc.z;
				cnt++;
			}
		}
	}
	return sampleLocs;
}

bool PointMonteCarloSampler::init(float portion_voxels, int neighborSize)
{
	if (neighborSize == 0){
		cout << "Neighborhood size cannot be 0" << endl;
		return false;
	}
	if (neighborSize % 2 == 0)
		neighborSize++; // make sure the neighbor size is odd
	_neighborSize = neighborSize;
	_sample_portion = portion_voxels;
	UINT64VECTOR3 dim = _volDim;
	_numTotalSamples = MIN(dim.volume(), UINT64(ceil(double(portion_voxels) * double(dim.volume()))));
	return true; 
}

bool PointMonteCarloSampler::init(INTVECTOR3 sample_dim, int neighborSize)
{
	if (neighborSize == 0){
		cout << "Neighborhood size cannot be 0" << endl;
		return false;
	}
	if (neighborSize % 2 == 0)
		neighborSize++; // make sure the neighbor size is odd
	_neighborSize = neighborSize;
	_numTotalSamples = sample_dim.volume();
	_sample_portion = float(_numTotalSamples) / float(_volDim.volume());
	return true;
}

bool PointMonteCarloSampler::updatePlot()
{

	if (!g_sampled_raw.empty()) // We have sampled the data
		return PointDataHandler::updatePlotFromLayerData(); 
	else // No sampling on the data is available
		return updatePlotFromVolumeSampling();
}

bool PointMonteCarloSampler::updatePlotFromVolumeSampling()
{
	if (_sampleCnt >= _numTotalSamples)
		return true; // we are done here!
	int startCnt = _sampleCnt;
	const int num_sample_per_update = 2000;
	// Otherwise, do sampling the volume 
	for (; _sampleCnt < MIN(_numTotalSamples, startCnt + num_sample_per_update); _sampleCnt++){

		//UINT64VECTOR3 loc = UINT64VECTOR3(_rng->randInt(UINT32(_volDim.x)), _rng->randInt(UINT32(_volDim.y)), _rng->randInt(UINT32(_volDim.z)));
		// Sample the neighborhood
		//sampleCubeNeighborhood(loc, _neighborSize);
		// Sample the normalized continuous location with neighborhood
		FLOATVECTOR3 loc = FLOATVECTOR3(_rng->rand(), _rng->rand(), _rng->rand());
		// Sample only on every 10 z slices
		//loc.z = float(_rng->randInt(9) * 10) / float(_volDim.z);
		sampleCubeNeighborhoodNormCoords(loc, _neighborSize);
	}
	return false;
}

bool PointMonteCarloSampler::updatePCPlot()
{
	const vector<vector<FLOATVECTOR2>>& streamlines = reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->getStreamLines();

	//for (vector<vector<FLOATVECTOR2>>::const_iterator IT = streamlines.begin(); IT != streamlines.end(); IT++)
	static vector<vector<FLOATVECTOR2>>::const_iterator IT = streamlines.begin();
	if (IT == streamlines.end())
		return true;

	{
		const vector<FLOATVECTOR2>& lines = *IT;
		IT++;
		if (lines.size() <= 2)
			return false; 

		vector<float> s1, s2;
		s1.resize(2, 0);
		s2.resize(2, 0);
		s1[0] = lines[0].x;
		s1[1] = lines[0].y;

		_pcpPlotBuilder->processSample(s1);
		for (size_t i = 1; i < lines.size(); i++){
			s2[0] = lines[i].x;
			s2[1] = lines[i].y;

			_pcpPlotBuilder->processSampleLineTo(s1, s2);
			//_pcpPlotBuilder->processSample(s2);
			s1 = s2;
		}
		
		return false;
	}
}

bool PointMonteCarloSampler::updateSCnPCPlots_from_streamline()
{

	// 1. SCPlotBuilder detects streamline
	static int seedId = 0;
	vector<FLOATVECTOR2> streamline;
	bool flipY = true; // Flip Y axis for the second attribute and so on...
	bool drawMore = reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->buildOneStreamline(streamline, seedId, g_scVec2DParam.max_pts_per_curve);
	if (!drawMore)
		return true; // we are done here
	else
	{ 

		if (streamline.size() <= 2)
			return false;

		vector<float> s1, s2;
		s1.resize(2, 0);
		s2.resize(2, 0);
		s1[0] = streamline[0].x;
		s1[1] = streamline[0].y;

		_pcpPlotBuilder->processSample(s1);
		for (size_t i = 1; i < streamline.size(); i++)
		{
			s2[0] = streamline[i].x;
			s2[1] = streamline[i].y;

			//_pcpPlotBuilder->processSampleLineTo(s1, s2);
			if (reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->processSampleLineToCheckPCPpointLoc(s1, s2))
			{
				_scPlotBuilder->setPenMode(false); // not using background pen
				// Redraw s1->s2 in the scatterplot
				bool emph = true;
				_scPlotBuilder->processSampleLineTo(s1, s2, emph);
				// Redraw s1->s2 in the PCP

			}
			//_pcpPlotBuilder->processSample(s2);
			s1 = s2;
		}
	}
	seedId++;
	return false; // need to draw more
}

void PointMonteCarloSampler::finishPlot(bool buildAllStreamlines)
{
	//// Write out the vector field
	//reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->writeOutVectorField("vecField");
	//// finish the plot
	//if ( buildAllStreamlines ) // build all streamlines
	//	reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->buildStreamlines(g_scVec2DParam.num_seeds, g_scVec2DParam.max_pts_per_curve);
	//else // Simply generate streamline seeds
	//	reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->genStreamlineSeeds(g_scVec2DParam.num_seeds);


}

void PointMonteCarloSampler::sampleCubeNeighborhood(UINT64VECTOR3 c, int neighborSize)
{
	int halfSize = (neighborSize - 1) / 2; 
	int zmin = MAX(0, int(c.z) - halfSize);
	int zmax = MIN(int(_volDim.z - 1), int(c.z) + halfSize);
	int ymin = MAX(0, int(c.y) - halfSize);
	int ymax = MIN(int(_volDim.y - 1), int(c.y) + halfSize);
	int xmin = MAX(0, int(c.x) - halfSize);
	int xmax = MIN(int(_volDim.x - 1), int(c.x) + halfSize);

	FLOATVECTOR2 centerVal; 
	for (int z = zmin; z <= zmax; z++)
	{
		for (int y = ymin; y <= ymax; y++)
		{
			for (int x = xmin; x <= xmax; x++)
			{
				// Sample randomly inside the grids
				FLOATVECTOR3 p;
				p.x = float(x) + _rng->rand(); // random inside the voxel
				p.y = float(y) + _rng->rand();
				p.z = float(z) + _rng->rand();
				vector<float> val = sampleVoxelWithInterp(p);
				if (x == c.x && y == c.y && z == c.z)
					centerVal = FLOATVECTOR2(val[0], val[1]);
			}
		}
	}
	// compute linear regression for these samples
	calcLinearRegScatterplot(centerVal);
}

void PointMonteCarloSampler::sampleCubeNeighborhoodNormCoords(FLOATVECTOR3 c_nl, int neighborSize1D)
{
	int halfSize1D = (neighborSize1D - 1) / 2; 
	FLOATVECTOR3 halfNeighborSize = FLOATVECTOR3(float(halfSize1D) / float(_volDim.x), float(halfSize1D) / float(_volDim.y), float(halfSize1D) / float(_volDim.z));
	FLOATVECTOR3 oneVoxelSize = FLOATVECTOR3(1.0f / float(_volDim.x), 1.0f / float(_volDim.y), 1.0f / float(_volDim.z));

	FLOATVECTOR2 centerVal;
	vector<float> val = sampleVoxelWithInterp(FLOATVECTOR3(c_nl.x * float(_volDim.x - 1), c_nl.y * float(_volDim.y - 1), c_nl.z * float(_volDim.z - 1)));
	// records all values inside the neighborhood
	UINT64 cnt = 0; // a counter keeps tracking where we are
	vector< vector<float> > valList; 
	valList.resize(neighborSize1D * neighborSize1D * neighborSize1D);

	// For now, we just use the first 2 components of val
	// TODO: make it to full neighboring pairs!
	centerVal = FLOATVECTOR2(val[0], val[1]);
	vector<float> vCenterVal = val;
	accumSampleDensity(centerVal);
	// Sample neighborhood randomly
	for (float z = MAX(c_nl.z - halfNeighborSize.z, 0.0f); z <= MIN(c_nl.z + halfNeighborSize.z, 1.0f); z += oneVoxelSize.z)
	{
		for (float y = MAX(c_nl.y - halfNeighborSize.y, 0.0f); y <= MIN(c_nl.y + halfNeighborSize.y, 1.0f); y += oneVoxelSize.y)
		{
			for (float x = MAX(c_nl.x - halfNeighborSize.x, 0.0f); x <= MIN(c_nl.x + halfNeighborSize.x, 1.0f); x += oneVoxelSize.x)
			{
				FLOATVECTOR3 p;
				p.x = x + _rng->rand() * oneVoxelSize.x;
				p.y = y + _rng->rand() * oneVoxelSize.y;
				p.z = z + _rng->rand() * oneVoxelSize.z; 
				val = sampleVoxelWithInterp(FLOATVECTOR3(p.x * float(_volDim.x - 1), p.y * float(_volDim.y - 1), p.z * float(_volDim.z - 1)));
				valList[cnt++] = val;
			}
		}
	}
	// Remove unused elements
	valList.erase(valList.begin() + cnt, valList.end());
	// compute linear regression for these samples
	calcCovMatrix(vCenterVal, valList);
	//calcLinearRegScatterplot(centerVal);
}

void PointMonteCarloSampler::accumSampleDensity(const FLOATVECTOR2& val)
{
	reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->accumSampleOnBuffer(val);
}

vector<float> PointMonteCarloSampler::sampleVoxelWithInterp(FLOATVECTOR3 p)
{
	vector<float> vals = trilinear(p);
	processSampleScatterplot(vals);
	return vals; 
}

vector<float> PointMonteCarloSampler::sampleVoxelWithInterpNormCoords(FLOATVECTOR3 c_nl)
{
	vector<float> val = trilinear(FLOATVECTOR3(c_nl.x * float(_volDim.x - 1), c_nl.y * float(_volDim.y - 1), c_nl.z * float(_volDim.z - 1)));
	return val;
}

vector<float> PointMonteCarloSampler::sampleVoxelNormCoords(FLOATVECTOR3 c_nl, const std::vector<VolumeData*>& volList)
{
	assert(!volList.empty());
	UINT64VECTOR3 volDim = volList[0]->getDim();
	vector<float> val = DirectionalMonteCarloSampler::trilinear(FLOATVECTOR3(c_nl.x * float(volDim.x - 1), c_nl.y * float(volDim.y - 1), c_nl.z * float(volDim.z - 1)), volList);
	return val;
}

vector<float> PointMonteCarloSampler::readVoxelNormalCoords(FLOATVECTOR3 c_nl, const std::vector<VolumeData*>& volList)
{
	assert(!volList.empty());
	UINT64VECTOR3 volDim = volList[0]->getDim();
	//UINT64VECTOR3 pos = UINT64VECTOR3(std::roundf(c_nl.x * float(volDim.x - 1)),
	//	std::roundf(c_nl.y * float(volDim.y - 1)), 
	//	std::roundf(c_nl.z * float(volDim.z - 1)));

	UINT64VECTOR3 pos = UINT64VECTOR3(
		MIN(volDim.x - 1, c_nl.x * float(volDim.x)),
		MIN(volDim.y - 1, c_nl.y * float(volDim.y)),
		MIN(volDim.z - 1, c_nl.z * float(volDim.z)));
#ifdef IDX_PT_DEBUG
	cout << c_nl << "ix ";
	cout << pos << ": ";
#endif
	//vector<float> val = DirectionalMonteCarloSampler::trilinear(, volList);
	int numCh = volList[0]->numChannels();
	vector<float> val(volList.size()*numCh);
	for (size_t i = 0; i < val.size(); i++) {
		for (int ch = 0; ch < numCh; ch++) {
			val[i * numCh + ch] = volList[i]->getVoxel(pos, ch);
#ifdef IDX_PT_DEBUG
			cout << val[i * numCh + ch] << " ";
#endif
		}	
	}
#ifdef IDX_PT_DEBUG
	cout << endl;
#endif
	return val;
}


void PointMonteCarloSampler::calcCovMatrix(const vector<float>& centralVal, const vector<vector<float>>& valList)
{
	// NOTE: for Eigen SVD decomposition
	// http://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
	Eigen::VectorXd mu;
	MyStatistics::ComputeMean<float>(valList, mu);
	//cout << "mu = " << endl<< mu << endl;
	Eigen::MatrixXd cov;
	MyStatistics::ComputeCovariance<float>(valList, mu, cov);
	//cout << "Cov = " << endl<< cov << endl;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
	// Left singular vectors
	//cout << "U= " <<endl<< svd.matrixU() << endl;
	// Get the major eigen vector
	//cout << "Eigen vectors = " << endl;

	////////////////////
	// SVD gives us A = USV(^T), where
	// U: each row of U is an eigen vector
	// S: each diagonal element is the eigen value
	
	int maxEigId = -1;
	float maxEigVal = -numeric_limits<float>::max();

	// NOTE: no need to sort eigen values!!!
	// Eigen make sure singular values are always sorted in decreasing order!!!
	//For the SVD decomposition of a n - by - p matrix, letting m be the minimum of n and p, the returned vector has size m.Singular values are always sorted in decreasing order.

	// For test, take the first 2-component of the major eigen vector
	Eigen::VectorXd majEigV = svd.matrixU().col(0) * svd.singularValues().row(0);
	Eigen::VectorXd secMajEigV = svd.matrixU().col(1) * svd.singularValues().row(1);

	Eigen::VectorXd v0 = majEigV + mu;
	Eigen::VectorXd v1 = secMajEigV + mu;
	//cout << "min eigvec =" << svd.matrixU().col(2) << endl;
	//reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->set2FlatNormDir(FLOATVECTOR2(centralVal[0], centralVal[1]), FLOATVECTOR2(majEigV(0), majEigV(1)));
	vector<double> v_EigV0d(majEigV.data(), majEigV.data() + majEigV.rows());
	vector<double> v_EigV1d(secMajEigV.data(), secMajEigV.data() + secMajEigV.rows());
	vector<float>  v_EigV0(v_EigV0d.size(), 0.0f);
	vector<float>  v_EigV1(v_EigV1d.size(), 0.0f);
	// c + v0, c + v1
	vector<float>  v_cplusEigV0(v_EigV0d.size(), 0.0f);
	vector<float>  v_cplusEigV1(v_EigV0d.size(), 0.0f);
	
	for (size_t i = 0; i < v_EigV0d.size(); i++)
	{
	
		v_EigV0[i] = float(v_EigV0d[i]);
		v_EigV1[i] = float(v_EigV1d[i]);

		v_cplusEigV0[i] = float(v0(i));
		v_cplusEigV1[i] = float(v1(i));
	}


	// Draw the original plot for centralVal
	// 1. Draw 0-flat (point in cartesian)
	reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->processSampleRepeatAxes(centralVal);
	// 2. Draw 1-flats (lines in cartesian)
	// Vec 1
	reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->processSampleLineToRepeatAxes(centralVal, v_cplusEigV0);
	// Vec 2
	reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->processSampleLineToRepeatAxes(centralVal, v_cplusEigV1);

	// Update other views
	if (_splomPlotBuilder)
		reinterpret_cast<SplomPlotBuilder*>(_splomPlotBuilder)->processSample(centralVal);
	if (_sc3DPlotBuilder)
		reinterpret_cast<ScatterPlot3DBuilder*>(_sc3DPlotBuilder)->processSample(centralVal);
	// right singular vectors
	//cout<<"V= "<<endl<<svd.matrixV()<<endl;
	// singular values
	//cout << "SingValues= " << endl << svd.singularValues() << endl;

	// 3. Draw 2-flat (plane in cartesian)
	// Create the plane
	PlaneCoeffs Pl;
	if (calcPlane(v_EigV0, v_EigV1, centralVal, Pl)){
		//cout << "Failed to create a plane!" << endl;
		reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->calc2FlatIndexPoints(Pl);
	}


}

bool PointMonteCarloSampler::calcPlane(const vector<float>& v1, const vector<float>& v2, const vector<float>& p, PlaneCoeffs& Pl)
{
	float delta = 1e-8f;
	// Compute the cross product of two vectors
	// TODO: need to extend to N-D
	// For now, use only 3-D vectors!!!
	if (v1.size() != 3 || v1.size() != v2.size())
		return false;

	Vector vv1(v1[0], v1[1], v1[2]);
	Vector vv2(v2[0], v2[1], v2[2]);
	Vector N = Cross(vv1, vv2);
	if (N.Length() <= delta)
		return false;
	Pl.c1 = N.x;
	Pl.c2 = N.y;
	Pl.c3 = N.z;
	
	Pl.c0 = N.x * p[0] + N.y * p[1] + N.z * p[2];
	return true;
}

void PointMonteCarloSampler::calcLinearRegScatterplot(FLOATVECTOR2 cVal)
{
	reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->calcLinearRegForNeighborhood(cVal);
}

void PointMonteCarloSampler::processSampleScatterplot(const vector<float>& vals)
{
	reinterpret_cast<ScatterPlotBuilder*>(_scPlotBuilder)->processAndStoreSample(vals);
}

void PointMonteCarloSampler::testPCP2Flat()
{
	// Test the construction of 2-flat in the PCP
	PlaneCoeffs test2flat;
	test2flat.c1 = /*0.25f; */6.0f;
	test2flat.c2 = /*0.532f;*/  5.0f;
	test2flat.c3 = /*0.216f;*/  8.0f;
	test2flat.c0 = /*0.631f;*/  8.0f; // coefficients for the plane
	vector<float> s(3, 0.0f);
	reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->init(s);
	reinterpret_cast<PCPbuilder*>(_pcpPlotBuilder)->construct2Flat(test2flat);
}

