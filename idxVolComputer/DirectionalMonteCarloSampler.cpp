#include "DirectionalMonteCarloSampler.h"
#include "MersenneTwister.h"
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "film.h"
#include "OrthoCamera.h"
#include "PlotBuilder.h"
#include "ScatterPlotBuilder.h"
#include "PCPbuilder.h"
#include "PCPLineDensityBuilder.h"
#include "ScatterPlot3DBuilder.h"
#include "SplomPlotBuilder.h"

using namespace std;

std::ostream& operator<<(std::ostream& os, const DirectionalSample& obj)
{
	os<<obj._pos<<"/"<<obj._dir<<"/"<<obj._len;
	return os;
}

std::istream& operator>>(std::istream& is, DirectionalSample& obj)
{
	is>>obj._pos>>obj._dir>>obj._len;
	return is;
}

DirectionalMonteCarloSampler::DirectionalMonteCarloSampler(const vector<VolumeData*>& volData, UINT64VECTOR3 volDim, const string& strOutputFile):
PointDataHandler(),
_scatterplotDim(g_params.getPCPDispBufSize().x, g_params.getPCPDispBufSize().y),
_volData(volData),
_volDim(volDim),
_strOutputFile(strOutputFile),
_sampleStepSize(0.1f),
_rayCheckpointSteps(10),
_film(NULL)
{
	_rng = new MTRand();
	_rng->seed();

	_scatterplot = new UINT64[_scatterplotDim.area()];
	memset(_scatterplot, 0, _scatterplotDim.area() * sizeof(UINT64));
	
	
	// Build world to texture transformation
	FLOATMATRIX4 w2tTrans, w2tScale;
	w2tTrans.Translation(0.5f, 0.5f, 0.5f);
	w2tScale.Scaling(float(_volDim.x) - 1, float(_volDim.y) - 1, float(_volDim.z) - 1);
	_worldToTex = w2tScale * w2tTrans; 

	// Scatter plot OR PCP?
	//_pcpPlotBuilder = new ScatterPlotBuilder(_scatterplotDim.y, _scatterplotDim.x);
	int plotDispBufW = g_params.getPCPDispBufSize().x;// g_pcpDispBufWidth;
	int plotDispBufH = g_params.getPCPDispBufSize().y;// g_pcpDispBufHeight
	////////////////////////////////////////////////////////////
	// PCP plot 
	switch (g_params.Pcp_RenderMode())
	{
	case PR_LINE_IDX_PT:
		_pcpPlotBuilder = new PCPbuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->flipY(true);
		break;
	case PR_DS_LINE_DENSITY:
		_pcpPlotBuilder = new PCPLineDensityBuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPLineDensityBuilder*>(_pcpPlotBuilder)->flipY(true);
		cout << "PCPLineDensityBuilder created!" << endl;
		break;
	default:
		_pcpPlotBuilder = new PCPbuilder(plotDispBufW, plotDispBufH);
		// Shall we flip Y????
		dynamic_cast<PCPbuilder*>(_pcpPlotBuilder)->flipY(true);
		break;
	}
	////////////////////////////////////////////////////////
	_scPlotBuilder = new ScatterPlotBuilder(g_splomDispBufWidth, g_splomDispBufHeight);


	// Get value range 
	int dataRange = 256;
	// If the type of input data is not unsigned char, it has already been normalized!!!
	// We now supports only unsigned data!!!
	if (_volData[0]->getVoxelByteSize() == 1) 
	{
		if (_volData[0]->isSigned())
		{
			_volValNormalizer = double(MAX_VOL_VAL_CHAR);
		}
		else
			_volValNormalizer = double(MAX_VOL_VAL_UCHAR);
		//dataRange = 256;
	}
	else if (_volData[0]->getVoxelByteSize() == 2)
	{
		if (_volData[0]->isSigned())
			_volValNormalizer = double(MAX_VOL_VAL_SHORT);
		else
			_volValNormalizer = double(MAX_VOL_VAL_USHORT);
		//dataRange = 512;
	}
	else if (_volData[0]->getVoxelByteSize() == 4 )
	{
		if (_volData[0]->isFloat())
			_volValNormalizer = 1.0;
		else
		{
			if (_volData[0]->isSigned())
				_volValNormalizer = double(MAX_VOL_VAL_INT);
			else
				_volValNormalizer = double(MAX_VOL_VAL_UINT);
		}
		//dataRange = 1024;
	}
	else
		_volValNormalizer = 1.0;

	// Create sc3d
	_sc3DPlotBuilder = NULL;
	//if (volData.size() == 3)
		_sc3DPlotBuilder = new ScatterPlot3DBuilder(dataRange, dataRange, dataRange, dataRange, dataRange, dataRange);
	// create splom
	_splomPlotBuilder = NULL;
	if (volData.size() >= 2 || g_datVolInfo._components >= 3) // either multiple scalar volumes or a vector volume is eligible for a splom plot builder
		_splomPlotBuilder = new SplomPlotBuilder(_scatterplotDim.x, _scatterplotDim.y);
	_isRedraw = false; // by default, we are not redrawing

	_valueDim = int(volData.size());
}


DirectionalMonteCarloSampler::~DirectionalMonteCarloSampler(void)
{
	SAFE_DELETE(_rng);
	SAFE_DELETE(_film);
}

void DirectionalMonteCarloSampler::sample(float portion_voxels, int lines_per_voxel)
{
	_samples.clear();
	generatesamples(_samples, portion_voxels, lines_per_voxel);
	sampleVolume(_samples);
	// output scatterplot
	writeResult();
}

void DirectionalMonteCarloSampler::generatesamples(vector<DirectionalSample>& samples, float portion_voxels, int lines_per_voxel)
{
	UINT64VECTOR3 dim = _volDim;
	FLOATVECTOR3  fdim = FLOATVECTOR3(dim.x, dim.y, dim.z);
	float diagonal = float(fdim.length());
	UINT64 num_samples = MIN(dim.volume(), UINT64(ceil(double(portion_voxels) * double(dim.volume()))));
	samples.resize(num_samples * lines_per_voxel);
	UINT64 cnt = 0;
	for(UINT64 i = 0; i < num_samples; i++ )
	{
		
		FLOATVECTOR3 loc = FLOATVECTOR3(_rng->rand(), _rng->rand(), _rng->rand());
		for(int d = 0; d < lines_per_voxel; d++ )
		{
			FLOATVECTOR3 dir = FLOATVECTOR3(_rng->rand(), _rng->rand(), _rng->rand());
			DirectionalSample ss;
			ss._dir = dir;
			ss._pos = UINT64VECTOR3(UINT64(float(dim.x - 1) * loc.x), UINT64(float(dim.y - 1) * loc.y), UINT64(float(dim.z - 1) * loc.z));
			ss._len = _rng->rand() * diagonal;
			samples[cnt++] = ss;
		}
	}
	// print out 
	//ofstream ofSamples("samples.txt");
	//for(vector<DirectionalSample>::iterator IT = samples.begin(); IT != samples.end(); ++IT)
	//	ofSamples<<*IT<<endl;
	//ofSamples.close();
}

void DirectionalMonteCarloSampler::sampleVolume(const vector<DirectionalSample>& samples)
{

	for(vector<DirectionalSample>::const_iterator IT = samples.begin(); IT != samples.end(); ++IT)
	{
		//sampleLine(IT->_pos, IT->_dir, IT->_len);
		sampleLineWithInterp(FLOATVECTOR3(IT->_pos.x, IT->_pos.y, IT->_pos.z), IT->_dir, IT->_len);
	//	temp = _matScatterplot + 1;
	//	cv::log(temp, logScatterplot);
	///*	logScatterplot.convertTo(ucLogScatterplot, CV_8UC1);
	//	cv::imshow("Scatterplot", ucLogScatterplot * 255);*/
	//	cv::imshow("scatterplot", logScatterplot/10.0);
	//	cv::waitKey(1);
	}
	cout<<"Done!"<<endl;
	//cv::waitKey(0);
}

void DirectionalMonteCarloSampler::sampleLine(UINT64VECTOR3 start, FLOATVECTOR3 dir, float len)
{
	FLOATVECTOR3 fstart = FLOATVECTOR3(start.x, start.y, start.z);
	FLOATVECTOR3 fend = fstart + dir * len;
	FLOATVECTOR3 dim = FLOATVECTOR3(_volDim.x, _volDim.y, _volDim.z);

	fend.x = MAX(0, MIN(dim.x-1, fend.x));
	fend.y = MAX(0, MIN(dim.y-1, fend.y));
	fend.z = MAX(0, MIN(dim.z-1, fend.z));

	FLOATVECTOR3 vec = fend - fstart;
	int N = vec.abs().maxVal();
	FLOATVECTOR3 s = vec / N;
	FLOATVECTOR3 p = FLOATVECTOR3(start.x, start.y, start.z);
	for(int i = 0; i < N; i++ )
	{
		UINT64VECTOR3 destVoxel = UINT64VECTOR3(round<float>(p.x), round<float>(p.y), round<float>(p.z));
		if( i == 0 )
			sampleVoxel(destVoxel);
		else
			lineToSampleVoxel(destVoxel);
		p = p + s;
	}
}


void DirectionalMonteCarloSampler::sampleLineWithInterp(FLOATVECTOR3 start, FLOATVECTOR3 dir, float len)
{
	FLOATVECTOR3 fstart = FLOATVECTOR3(start.x, start.y, start.z);
	FLOATVECTOR3 fend = fstart + dir * len;
	FLOATVECTOR3 dim = FLOATVECTOR3(_volDim.x, _volDim.y, _volDim.z);

	fend.x = MAX(0, MIN(dim.x-1, fend.x));
	fend.y = MAX(0, MIN(dim.y-1, fend.y));
	fend.z = MAX(0, MIN(dim.z-1, fend.z));

	FLOATVECTOR3 vec = fend - fstart;
	int N = vec.length() / _sampleStepSize; //vec.abs().maxVal();
	cout<<"Ray length = "<<vec.length()<<" voxels."<<endl;
	FLOATVECTOR3 s = vec / N; // incremental vector for each step
	FLOATVECTOR3 p = FLOATVECTOR3(start.x, start.y, start.z);

	bool drawWithEmph = false;
	_isRedraw = false; 
	for(int i = 0; i < N; i++ )
	{
		FLOATVECTOR3 destVoxel = p;
		//UINT64VECTOR3 destVoxel = UINT64VECTOR3(round<float>(p.x), round<float>(p.y), round<float>(p.z));
		if( i == 0 )
			sampleVoxelWithInterp(destVoxel);
		else
			lineToSampleVoxelWithInterp(destVoxel);

		if (i % _rayCheckpointSteps == 0 && i > 0)
		{
			if (_isRedraw)
				_isRedraw = false; // when we are redrawing, set flag to false, and keep going
			else
			{
				drawWithEmph = processCheckpoint();
				if (drawWithEmph)
				{
					p = p - float(_rayCheckpointSteps) * s;  // If we need to redraw, then step back to the last check point and redraw!
					_isRedraw = true;
				}
			}
		}
		else
			p = p + s;
	}
}

bool DirectionalMonteCarloSampler::processCheckpoint()
{
	bool needToRedraw = _pcpPlotBuilder->processCheckpoint();
	bool isBkgrPen = !needToRedraw;
	_scPlotBuilder->setPenMode(isBkgrPen);
	return needToRedraw;
}

void DirectionalMonteCarloSampler::lineTo2D(UINT64VECTOR2 start, UINT64VECTOR2 end)
{
	FLOATVECTOR2 fstart = FLOATVECTOR2(start.x, start.y);
	FLOATVECTOR2 fend = FLOATVECTOR2(end.x, end.y);
	FLOATVECTOR2 vec = fend - fstart;
	int N = vec.abs().maxVal();
	FLOATVECTOR2 s = vec / N;
	FLOATVECTOR2 p = fstart;
	for(int i = 0; i < N; i++ )
	{
		UINT64VECTOR2 ip = UINT64VECTOR2(round<float>(p.x), round<float>(p.y));
		incrScatterplotElemCnt(ip);
		p = p + s;
	}

}

void DirectionalMonteCarloSampler::incrScatterplotElemCnt(UINT64VECTOR2 p)
{
	UINT64 id = p.y * _scatterplotDim.x + p.x;
	_scatterplot[id]++;

}

void DirectionalMonteCarloSampler::incrScatterplotElemCnt(FLOATVECTOR2 p)
{
	UINT64 id = (p.y * (_scatterplotDim.y - 1)) * _scatterplotDim.x + p.x * (_scatterplotDim.x - 1);
	_scatterplot[id]++;
}

void DirectionalMonteCarloSampler::sampleVoxel(UINT64VECTOR3 p)
{
	// For now we just use two volumes...later add more to build parallel coord
	unsigned int val = _volData[0]->getVoxel(p); 
	unsigned int gm = _volData[1]->getVoxel(p);

	float x, y;
	if( _volData[0]->getVoxelByteSize() == 1)
	{
		x = float(val) / float(numeric_limits<unsigned char>::max());
		y = float(gm) / float(numeric_limits<unsigned char>::max());
	}
	else if( _volData[0]->getVoxelByteSize() == 2 )
	{
		x = float(val) / float(numeric_limits<unsigned int>::max());
		y = float(gm) / float(numeric_limits<unsigned int>::max());
	}
	UINT64 iy = UINT64(y * float(_scatterplotDim.y - 1));
	UINT64 ix = UINT64(x * float(_scatterplotDim.x - 1));
	incrScatterplotElemCnt(UINT64VECTOR2(ix, iy));
	
	_lastPosScatterplot = UINT64VECTOR2(ix, iy);
	_fLastPosScatterplot = FLOATVECTOR2(x,y);
}

void DirectionalMonteCarloSampler::sampleVoxelWithInterp(FLOATVECTOR3 p)
{
	vector<float> vals = trilinear(p);

	// Let the plot builder do its job
	processSample(vals);


	//// Take only the first two elements for now
	//float x = vals[0];
	//float y = vals[1];
	//if( _volData[0]->getVoxelByteSize() == 1)
	//{
	//	x = float(x) / float(numeric_limits<unsigned char>::max());
	//	y = float(y) / float(numeric_limits<unsigned char>::max());
	//}
	//else if( _volData[0]->getVoxelByteSize() == 2 )
	//{
	//	x = float(x) / float(numeric_limits<unsigned int>::max());
	//	y = float(y) / float(numeric_limits<unsigned int>::max());
	//}
	//UINT64 iy = UINT64(y * float(_scatterplotDim.y - 1));
	//UINT64 ix = UINT64(x * float(_scatterplotDim.x - 1));


	//incrScatterplotElemCnt(UINT64VECTOR2(ix, iy));

	//_lastPosScatterplot = UINT64VECTOR2(ix, iy);
	//_fLastPosScatterplot = FLOATVECTOR2(x,y);

}

void DirectionalMonteCarloSampler::processSample(const vector<float>& vals)
{
	_pcpPlotBuilder->processSample(vals);
	_scPlotBuilder->processSample(vals);
	_lastSampleVal = vals;
}

void DirectionalMonteCarloSampler::processSampleLineTo(const vector<float>& vals)
{
	_pcpPlotBuilder->processSampleLineTo(_lastSampleVal, vals);
	_scPlotBuilder->processSampleLineTo(_lastSampleVal, vals);
	_lastSampleVal = vals;
}

void DirectionalMonteCarloSampler::lineToSampleVoxel(UINT64VECTOR3 p)
{
	// For now we just use two volumes...later add more to build parallel coord
	unsigned int val = _volData[0]->getVoxel(p); 
	unsigned int gm = _volData[1]->getVoxel(p);

	float x, y;
	if( _volData[0]->getVoxelByteSize() == 1)
	{
		x = float(val) / float(numeric_limits<unsigned char>::max());
		y = float(gm) / float(numeric_limits<unsigned char>::max());
	}
	else if( _volData[0]->getVoxelByteSize() == 2 )
	{
		x = float(val) / float(numeric_limits<unsigned int>::max());
		y = float(gm) / float(numeric_limits<unsigned int>::max());
	}
	UINT64 iy = UINT64(y * float(_scatterplotDim.y - 1));
	UINT64 ix = UINT64(x * float(_scatterplotDim.x - 1));
	UINT64VECTOR2 dest2D(ix, iy);

	lineTo2D(_lastPosScatterplot, dest2D);

	_lastPosScatterplot = dest2D;
}

void DirectionalMonteCarloSampler::lineToSampleVoxelWithInterp(FLOATVECTOR3 p)
{
	vector<float> vals = trilinear(p);
	// Let the plot builder do its job
	processSampleLineTo(vals);

	//// Take only the first two elements for now
	//float x = vals[0];
	//float y = vals[1];
	//if( _volData[0]->getVoxelByteSize() == 1)
	//{
	//	x = float(x) / float(numeric_limits<unsigned char>::max());
	//	y = float(y) / float(numeric_limits<unsigned char>::max());
	//}
	//else if( _volData[0]->getVoxelByteSize() == 2 )
	//{
	//	x = float(x) / float(numeric_limits<unsigned int>::max());
	//	y = float(y) / float(numeric_limits<unsigned int>::max());
	//}
	//UINT64 iy = UINT64(y * float(_scatterplotDim.y - 1));
	//UINT64 ix = UINT64(x * float(_scatterplotDim.x - 1));
	//UINT64VECTOR2 dest2D(ix, iy);

	//lineTo2D(_lastPosScatterplot, dest2D);

	//_lastPosScatterplot = dest2D;
}

void DirectionalMonteCarloSampler::showPlot()
{
	_pcpPlotBuilder->show();
	_scPlotBuilder->show();
}

void DirectionalMonteCarloSampler::clearWorkingPlotBuf()
{
	_pcpPlotBuilder->clearTempBuf();
	_scPlotBuilder->clearTempBuf();
}

void DirectionalMonteCarloSampler::writeResult()
{
	// writeout scatterplot
	ofstream ofScatterplot(_strOutputFile.c_str());
	for(UINT64 y = 0; y < _scatterplotDim.y; y++ )
	{
		for(UINT64 x = 0; x < _scatterplotDim.x; x++ )
			ofScatterplot<<_scatterplot[y * _scatterplotDim.x + x]<<" ";
		ofScatterplot<<endl;
	}
	ofScatterplot.close();
}

vector<float> DirectionalMonteCarloSampler::trilinear(FLOATVECTOR3 p)
{
	vector<float> normval =  trilinear(p, _volData);
	for(size_t i = 0; i < normval.size(); i++)
		normval[i] /= _volValNormalizer;
	return normval;
}

std::vector<float> DirectionalMonteCarloSampler::trilinear(FLOATVECTOR3 p, const std::vector<VolumeData*>& volList)
{
	// Also supports vector volume data
	assert(volList.size() > 0);
	int numchannels = volList[0]->numChannels();
	vector<float> val(volList.size() * numchannels, 0.0f);
	UINT64VECTOR3 volDim = volList[0]->getDim();
	float gx, gy, gz, tx, ty, tz;
	UINT64 gxi, gyi, gzi;
	gx = p.x;
	gy = p.y;
	gz = p.z;

	gxi = UINT64(p.x);
	gyi = UINT64(p.y);
	gzi = UINT64(p.z);

	tx = gx - gxi;
	ty = gy - gyi;
	tz = gz - gzi;

	float c000, c100, c010, c110, c001, c101, c011, c111;

	UINT64VECTOR3 p000 = UINT64VECTOR3(gxi, gyi, gzi);
	UINT64VECTOR3 p100 = UINT64VECTOR3(MIN(volDim.x - 1, gxi + 1), gyi, gzi);
	UINT64VECTOR3 p010 = UINT64VECTOR3(gxi, MIN(volDim.y - 1, gyi + 1), gzi);
	UINT64VECTOR3 p110 = UINT64VECTOR3(MIN(volDim.x - 1, gxi + 1), MIN(volDim.y - 1, gyi + 1), gzi);
	UINT64VECTOR3 p001 = UINT64VECTOR3(gxi, gyi, MIN(volDim.z - 1, gzi));
	UINT64VECTOR3 p101 = UINT64VECTOR3(MIN(volDim.x - 1, gxi + 1), gyi, MIN(volDim.z - 1, gzi + 1));
	UINT64VECTOR3 p011 = UINT64VECTOR3(gxi, MIN(volDim.y - 1, gyi + 1), MIN(volDim.z - 1, gzi + 1));
	UINT64VECTOR3 p111 = UINT64VECTOR3(MIN(volDim.x - 1, gxi + 1), MIN(volDim.y - 1, gyi + 1), MIN(volDim.z - 1, gzi + 1));
	for (size_t i = 0; i < volList.size(); i++)
	{
		for (int d = 0; d < numchannels; d++)
		{
			int j = i * numchannels + d;
			c000 = volList[i]->getVoxel(p000, d);
			c100 = volList[i]->getVoxel(p100, d);

			c010 = volList[i]->getVoxel(p010, d);
			c110 = volList[i]->getVoxel(p110, d);

			c001 = volList[i]->getVoxel(p001, d);
			c101 = volList[i]->getVoxel(p101, d);

			c011 = volList[i]->getVoxel(p011, d);
			c111 = volList[i]->getVoxel(p111, d);

			// interpolate
			val[j] = (1.0f - tx)*(1.0f - ty)*(1.0f - tz)*c000 +
				tx*(1.0f - ty)*(1.0f - tz)*c100 +
				(1.0f - tx)*ty*(1.0f - tz)*c010 +
				tx*ty*(1.0f - tz)*c110 +
				(1.0f - tx)*(1.0f - ty)*tz*c001 +
				tx*(1.0f - ty)*tz*c101 +
				(1.0f - tx)*ty*tz*c011 +
				tx*ty*tz*c111;

		}
	}

	return val;
}

void DirectionalMonteCarloSampler::generateSphereDirectionalSamples(std::vector<FLOATVECTOR3>& samples, UINT64 num_samples)
{

	UINT64VECTOR3 dim = _volDim;
	FLOATVECTOR3  fdim = FLOATVECTOR3(dim.x, dim.y, dim.z);

	samples.clear();
	samples.resize(num_samples);
	UINT64 cnt = 0;
	for(UINT64 i = 0; i < num_samples; i++ )
	{
		float u1 = _rng->rand();
		float u2 = _rng->rand();

		FLOATVECTOR3 sample = uniformSampleFullsphere(u1, u2);
		samples[i] = sample;
	}
}

void DirectionalMonteCarloSampler::raySampleVolume(FLOATVECTOR3 eyePos)
{

}

FLOATVECTOR3 DirectionalMonteCarloSampler::uniformSampleHemisphere(float u1, float u2)
{
	float z = u1;
	float r = sqrtf(MAX(0.0f, 1.0f - z * z));
	float phi = 2 * M_PI * u2;
	float x = r * cosf(phi);
	float y = r * sinf(phi);
	return FLOATVECTOR3(x, y, z);
}

FLOATVECTOR3 DirectionalMonteCarloSampler::uniformSampleFullsphere(float u1, float u2)
{
	float z = 1.0f - 2.0f * u1;
	float r = sqrtf(MAX(0.0f, 1.0f - z * z));
	float phi = 2 * M_PI * u2;
	float x = r * cosf(phi);
	float y = r * sinf(phi);
	return FLOATVECTOR3(x, y, z);
}

float DirectionalMonteCarloSampler::uniformHemispherePdf()
{
	return 1.0f / (2.0f * M_PI);
}

float DirectionalMonteCarloSampler::uniformFullspherePdf()
{
	return 1.0f / (4.0f * M_PI);
}

// Ray tracing based sampling...
void DirectionalMonteCarloSampler::rayTracingSampling(int num_dir, int sample_img_width, int sample_img_height)
{
	// Prepare film
	if (_film == NULL)
		_film = new Film(sample_img_width, sample_img_height);
	else
	{
		if (_film->xResolution != sample_img_width ||
			_film->yResolution != sample_img_height)
		{
			delete _film;
			_film = new Film(sample_img_width, sample_img_height);
		}
	}
	// generate samples on the film
	vector<Sample> filmSamples;
	bool jitter = false;
	_film->generateFilmSamples(filmSamples, jitter);
	// Generate camera locations
	vector<FLOATVECTOR3> cameraLocs;
	generateSphereDirectionalSamples(cameraLocs, UINT64(num_dir));
	// Get bounding box of the volume
	Point aabbMin(-0.5f, -0.5f, -0.5f);
	Point aabbMax(0.5f, 0.5f, 0.5f);
	BBox aabb(aabbMin, aabbMax);


	FLOATMATRIX4 w2tTrans, w2tScale;
	w2tTrans.Translation(-aabbMin.x, -aabbMin.y, -aabbMin.z);
	cout << "Trans" << endl << w2tTrans << endl;
	w2tScale.Scaling(float(_volDim.x) - 1, float(_volDim.y) - 1, float(_volDim.z) - 1);
	cout << "Scale" << endl << w2tScale << endl;
	FLOATVECTOR3 worldOrg(0.0f, 0.0f, 0.0f);
	FLOATVECTOR3 up(0.0f, 1.0f, 0.0f);
	float znear = 0.0f;
	float zfar = 1e3f;
	float aspectRatio = _film->xResolution / _film->yResolution;
	Ray   ray;
	float rayStep = _sampleStepSize / _volDim.length();
	// For each camera location, do ray casting
	for(size_t i = 0; i < cameraLocs.size(); i++ )
	{
		FLOATMATRIX4 worldToCamera;
		worldToCamera.BuildLookAt(cameraLocs[i], worldOrg, up);
		cout<<"World2Cam"<<endl<<worldToCamera<<endl;
		FLOATVECTOR4 worldOrgInCam = worldToCamera * FLOATVECTOR4(worldOrg.x, worldOrg.y, worldOrg.z, 1.0f);
		cout<<"World Org in Cam"<<endl<<worldOrgInCam<<endl;
		// get inverse 
		FLOATMATRIX4 cameraToWorld = worldToCamera.inverse();
		cout<<"World Org :"<<cameraToWorld * worldOrgInCam<<endl;
		OrthoCamera camera(worldToCamera, aspectRatio, znear, zfar, _film);
		for(vector<Sample>::iterator IT = filmSamples.begin(); IT != filmSamples.end(); ++IT)
		{
			camera.GenerateRay(*IT, &ray);
			// Test ray-aabb intersection
			float hitt0 = -INFINITY;
			float hitt1 = INFINITY;
			if( aabb.IntersectP(ray, &hitt0, &hitt1) )
			{

				// sampling from hitt0 to hitt1 (These are already in the world coordinates!!!)
				Point worldStart = ray(hitt0);
				Point worldEnd = ray(hitt1);
			
				FLOATVECTOR4 texStart = _worldToTex * FLOATVECTOR4(worldStart.x, worldStart.y, worldStart.z, 1.0f);
				FLOATVECTOR4 texEnd = _worldToTex * FLOATVECTOR4(worldEnd.x, worldEnd.y, worldEnd.z, 1.0f); 
				FLOATVECTOR3 texDir = texEnd.xyz()/texEnd.w - texStart.xyz()/texStart.w;
				float len = texDir.length();
				texDir.normalize();

				//clearWorkingPlotBuf();
				//for(float t = hitt0; t <= hitt1; t += rayStep)
				sampleLineWithInterp(texStart.xyz()/texStart.w, texDir, len);
			}

		}
	}
	cout<<"Done!"<<endl;
}
