#ifndef VOLUME_DATA_H
#define VOLUME_DATA_H
#include "MyVectors.h"
#include <limits>
using namespace std;

const int MAX_VOL_VAL_INT = numeric_limits<int>::max();
const int MAX_VOL_VAL_UINT = numeric_limits<unsigned int>::max();
const int MAX_VOL_VAL_SHORT = numeric_limits<short>::max();
const int MAX_VOL_VAL_USHORT = numeric_limits<unsigned short>::max();
const int MAX_VOL_VAL_CHAR = numeric_limits<char>::max();
const int MAX_VOL_VAL_UCHAR = numeric_limits<unsigned char>::max();


struct voxelTypeDesc{
	int byteSize;
	bool isSigned;
	bool isFloat;
};

class VolumeInfo
{
public:
	VolumeInfo() :
		_dim(0, 0, 0)
	{
		_voxType.byteSize = 1;
		_voxType.isSigned = false;
		_voxType.isFloat = false;
	}

	VolumeInfo(voxelTypeDesc voxType, UINT64VECTOR3 dim) :
		_voxType(voxType),
		_dim(dim)
	{}

	voxelTypeDesc _voxType;
	UINT64VECTOR3 _dim;
};

class VolumeData
{
public:
	VolumeData();
	VolumeData(UINT64VECTOR3 dim, int voxelByteSize, int value, bool isSigned = false, bool isFloat = false, int numChannels = 1);
	VolumeData(UINT64VECTOR3 dim, int voxelByteSize, void* data, bool isSigned = false, bool isFloat = false, int numChannels = 1);
	VolumeData(UINT64VECTOR3 dim, int voxelByteSize, const char* filename, bool isSigned = false, bool isFloat = false, int numChannels = 1);
	virtual ~VolumeData();

	void setDataPtr(void* ptr) { _data = ptr; }
	void setDataPtrNull() {_data = NULL;}
	// Added from the ColorPlayground project
	void setValueAllVoxels(int val); // set all voxels 
	void castDataToFloat(float** dstPtr);

	void getSubVol(void* buf, UINT64VECTOR3 pos, UINT64VECTOR3 sz) const;
	void setSubVol(const void* buf, UINT64VECTOR3 pos, UINT64VECTOR3 sz);
	// Accessors
	int getVoxelByteSize() const { return _voxelByteSize; }
	bool isSigned() const { return _isSigned; }
	bool isFloat() const { return _isFloat; }
	int  numChannels() const { return _numChannels; }

	void setSigned(bool isSigned) {_isSigned = isSigned; }
	void setFloat(bool isFloat) {_isFloat = isFloat; }

	float        getVoxel(UINT64VECTOR3 pos, int channel = 0) const;  // Get a scalar at channel "channel" at voxel "pos"
	void         getVec(UINT64VECTOR3 pos, std::vector<float>& vecVal) const;
	void         setVoxel(UINT64VECTOR3 pos, float value, int channel = 0); // Set a normalized scalar at channel "channel" at voxel pos. Note that value is a normalized float, in [0,1]
	void         getLayer(void* buf, UINT64 z) const; // get entire Z layer	
	void         setLayer(const void* buf, UINT64 z); // set entire Z layer

	void         getSlice(void* buf, UINT64 sliceId, int sliceDir) const; // Get slice given slice Id and slice Dir

	UINT64VECTOR3 getDim() const { return _dim; }

	//////////////////////////////////////////
	// Normalization
	void normalize(float& minVal, float& maxVal);
	void findMinMax(float& minVal, float& maxVal);
	void normalizeToRange(const float minVal, const float maxVal);
	// Acquire min/max for each channel
	void getMinMax(std::vector<float>& minVals, std::vector<float>& maxVals);
	void setMinMax(float minVal, float maxVal) { _minval = minVal; _maxval = maxVal; };
	//////////////////////////////////////////
	// Differential operations
	/////////////////////////////////////////
	// compute gradient magnitude:
	VolumeData* computeGradientMagnitude();
	// compute second derivative magnitude
	VolumeData* computeSecondDerivMag(); 
	// compute gradient, gradient magnitude and also directional second derivative
	void computeGradMagandDirSecondDeriv(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool normalize = true);
	void genSyntheticVols(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool normalize = true);
	void genSynGoodPatternVols(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool normalize = true);
	float computeDataSampleOnPlane(const PlaneCoeffs& P, const float& x, const float& y); // Given coefficients of a plane and two elements of a value, and compute the third element
	void  rotateDataSampleToPlane(const PlaneCoeffs& P, float& x, float& y, float& z);
	/////////////////////////////////////////

	// create a 1/8 volume
	VolumeData* downSampleOneLevel();

	void* getDataBuf() const { return _data; }
	// histograms
	void computeHistogram1D();
	void computeHistogram2D(VolumeData* gradMagVol);

	int  findBackgrndVal(); // Finds the background value in integer
	// Converte to float (by default normalize values as well!)
	void convert2Float(float** ppFloatVol, bool normalize = true) const; 

	bool writeToFile(const char* filename, char* finalName);
	void writePPMFiles(const char* filename);
	bool writeToNrrdFile(const char* filename);
	void writeLayerToPPMFile(const char* filename, UINT64 z);

	bool loadFromFile(const char* filename);
	bool loadFromNrrdFile(const char* filename, bool convertToFloat= true, bool doNormalize = false);
private:
	std::string            _inputFileName;
	void*                  _data;
	UINT64VECTOR3          _dim;
	std::vector<UINT64>    _histogram1D;
	int                    _voxelByteSize;
	int                    _backgrndVal;
	UINT64*                _histogram2D;
	int                      _nrrdType;
	bool                   _isSigned;
	bool                   _isFloat;

	// if it's vector data
	int                    _numChannels;
	std::vector<float>     _minVals;
	std::vector<float>     _maxVals;

	// For scalar data only...
	float                  _minval;
	float                  _maxval;
};

#endif