#include "VolumeData.h"
#include "Image.h"
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include "SysTools.h"
#include "Eigen\Dense"
#include "geometry.h"
#include "teem\nrrd.h"
#pragma warning(disable : 4996)
using namespace std;
using namespace Eigen;

VolumeData::VolumeData() :
	_dim(UINT64VECTOR3(0, 0, 0)),
	_voxelByteSize(0),
	_data(NULL),
	_isSigned(false),
	_isFloat(false),
	_backgrndVal(-1),
	_histogram2D(NULL),
	_maxval(-FLT_MAX),
	_minval(FLT_MAX),
	_nrrdType(nrrdTypeUChar),
	_numChannels(1)
{}


VolumeData::VolumeData(UINT64VECTOR3 dim, int voxelByteSize, int value, bool isSigned, bool isFloat, int numChannels):
_dim(dim),
_voxelByteSize(voxelByteSize),
_isSigned(isSigned),
_isFloat(isFloat),
_backgrndVal(-1),
_histogram2D(NULL),
_numChannels(numChannels)
{
	UINT64 allocSize = UINT64(_voxelByteSize) * _dim.volume() * UINT64(_numChannels);
	_data = new unsigned char[allocSize];
	memset(_data, value, allocSize);
}

VolumeData::VolumeData(UINT64VECTOR3 dim, int voxelByteSize, const char* filename, bool isSigned, bool isFloat, int numChannels):
_dim(dim),
_voxelByteSize(voxelByteSize),
_isSigned(isSigned),
_isFloat(isFloat),
_backgrndVal(-1),
_histogram2D(NULL),
_numChannels(numChannels)
{
	ifstream ifOr(filename, ios::in | ios::binary);
	assert(ifOr.is_open());
	//_data = malloc(_voxelByteSize*_dimX*_dimY*_dimZ);
	UINT64 allocSize = UINT64(_voxelByteSize) * _dim.volume() * UINT64(_numChannels);
	_data = malloc(allocSize);
	if( _data == NULL )
	{
		cout<<"cannot allocate data buffer!"<<endl;
	}
	else
	{
		memset(_data, 0, allocSize);
		ifOr.read(reinterpret_cast<char*>(_data), allocSize);
	}
	ifOr.close();
	assert(_data);
}

VolumeData::VolumeData(UINT64VECTOR3 dim, int voxelByteSize, void* data, bool isSigned, bool isFloat, int numChannels):
_dim(dim),
_voxelByteSize(voxelByteSize),
_isSigned(isSigned),
_isFloat(isFloat),
_backgrndVal(-1),
_histogram2D(NULL),
_numChannels(numChannels)
{
	UINT64 dataSize = UINT64(_voxelByteSize) * _dim.volume() * UINT64(_numChannels);
	_data = new unsigned char[dataSize];
#ifdef WIN32
	memcpy_s(_data, dataSize, data, dataSize);
#else
	memcpy(_data, data, dataSize);
#endif
}

VolumeData::~VolumeData()
{
	SAFE_DELETE_ARRAY(_data);
	SAFE_DELETE_ARRAY(_histogram2D);
}

void VolumeData::normalize(float& minVal, float& maxVal)
{
	minVal = numeric_limits<float>::max();
	maxVal = -numeric_limits<float>::max();
	// For vector data, simply normalize all to unit vector
	if (_numChannels > 1)
	{
		for (UINT64 z = 0; z < _dim.z; z++)
		{
			for (UINT64 y = 0; y < _dim.y; y++)
			{
				for (UINT64 x = 0; x < _dim.x; x++)
				{
					float vecLen = 0.0f;
					UINT64VECTOR3 loc(x, y, z);
					
					vector<float> vecVal;
					getVec(loc, vecVal);
					for (size_t c = 0; c < vecVal.size(); c++)
						vecLen += vecVal[c] * vecVal[c];
					vecLen = sqrtf(vecLen);
					// normalize to unit vector
					if (vecLen > 0.0f){
						for (size_t c = 0; c < vecVal.size(); c++)
							setVoxel(loc, vecVal[c] / vecLen, c);
					}
	
					minVal = MIN(minVal, vecLen);
					maxVal = MAX(maxVal, vecLen);
				}
			}
		}
	}
	else
	{
		// For scalar data, normalize based on min,max
		// First pass find min,max
		for (UINT64 z = 0; z < _dim.z; z++)
		{
			for (UINT64 y = 0; y < _dim.y; y++)
			{
				for (UINT64 x = 0; x < _dim.x; x++)
				{
					float voxVal = 0.0f;
					UINT64VECTOR3 loc(x, y, z);
					voxVal = getVoxel(loc);
					minVal = MIN(minVal, voxVal);
					maxVal = MAX(maxVal, voxVal);
				}
			}
		}
		// Second pass, normalization
		float oneOverRange = (minVal == maxVal) ? 0.0f : 1.0f / (maxVal - minVal);
		for (UINT64 z = 0; z < _dim.z; z++)
		{
			for (UINT64 y = 0; y < _dim.y; y++)
			{
				for (UINT64 x = 0; x < _dim.x; x++)
				{
					float voxVal = getVoxel(UINT64VECTOR3(x, y, z));
					float normVal = (voxVal - minVal) * oneOverRange;
					setVoxel(UINT64VECTOR3(x, y, z), normVal);
				}
			}
		}
	}

}

void VolumeData::getMinMax(vector<float>& minVals, vector<float>& maxVals)
{
	minVals.resize(_numChannels, numeric_limits<float>::max());
	maxVals.resize(_numChannels,-numeric_limits<float>::max());
	// First pass find min,max
	for (UINT64 z = 0; z < _dim.z; z++)
	{
		for (UINT64 y = 0; y < _dim.y; y++)
		{
			for (UINT64 x = 0; x < _dim.x; x++)
			{
				for (int ch = 0; ch < _numChannels; ch++)
				{
					float voxVal = getVoxel(UINT64VECTOR3(x, y, z), ch);
					minVals[ch] = MIN(minVals[ch], voxVal);
					maxVals[ch] = MAX(maxVals[ch], voxVal);
				}
	
			}
		}
	}

	_minVals = minVals;
	_maxVals = maxVals;
}

float VolumeData::getVoxel(UINT64VECTOR3 pos, int channel) const
{
	assert(pos.x < _dim.x && pos.y < _dim.y && pos.z < _dim.z);
	assert(_voxelByteSize == 1 || _voxelByteSize == 2 || _voxelByteSize == 4);
	assert(channel >= 0 && channel < _numChannels);
	
	UINT64 ind = UINT64(_numChannels) * ((pos.z * _dim.y + pos.y) * _dim.x + pos.x) + UINT64(channel);
	if (_voxelByteSize == 4)
	{
		if (_isFloat)
			return reinterpret_cast<float*>(_data)[ind];
		else
		{
			if (_isSigned)
				return float(reinterpret_cast<int*>(_data)[ind]);
			else
				return float(reinterpret_cast<unsigned int*>(_data)[ind]);
		}
	}
	else if (_voxelByteSize == 2)
	{
		if (_isSigned)
			return float(reinterpret_cast<short*>(_data)[ind]);
		else
			return float(reinterpret_cast<unsigned short*>(_data)[ind]);
	}

	else if (_voxelByteSize == 1)
	{
		if (_isSigned)
			return float(reinterpret_cast<char*>(_data)[ind]);
		else
			return float(reinterpret_cast<unsigned char*>(_data)[ind]);
	}
	else
		return 0.0f;
}

void VolumeData::getVec(UINT64VECTOR3 pos, vector<float>& vecVal) const
{
	assert(pos.x < _dim.x && pos.y < _dim.y && pos.z < _dim.z);
	assert(_voxelByteSize == 1 || _voxelByteSize == 2 || _voxelByteSize == 4);

	UINT64 ind = UINT64(_numChannels) * ((pos.z * _dim.y + pos.y) * _dim.x + pos.x);
	vecVal.clear();
	vecVal.resize(_numChannels);
	if (_voxelByteSize == 4)
	{
		if (_isFloat)
		{
			for (int i = 0; i < _numChannels; i++)
				vecVal[i] = reinterpret_cast<float*>(_data)[ind + i];			
		}
		else
		{
			if (_isSigned)
			{
				for (int i = 0; i < _numChannels; i++)
					vecVal[i] = float(reinterpret_cast<int*>(_data)[ind + i]);
			}
			else
			{
				for (int i = 0; i < _numChannels; i++)
					vecVal[i] = float(reinterpret_cast<unsigned int*>(_data)[ind + i]);
			}
		}
	}
	else if (_voxelByteSize == 2)
	{
		if (_isSigned)
		{
			for (int i = 0; i < _numChannels; i++)
				vecVal[i] = float(reinterpret_cast<short*>(_data)[ind + i]);
		}
		else
		{
			for (int i = 0; i < _numChannels; i++)
				vecVal[i] = float(reinterpret_cast<unsigned short*>(_data)[ind + i]);
		}
	}

	else if (_voxelByteSize == 1)
	{
		if (_isSigned)
		{
			for (int i = 0; i < _numChannels; i++)
				vecVal[i] = float(reinterpret_cast<char*>(_data)[ind + i]);
		}
		else
		{
			for (int i = 0; i < _numChannels; i++)
				vecVal[i] = float(reinterpret_cast<unsigned char*>(_data)[ind + i]);
		}
	}
	else
		cout << "Format not supported!" << endl;
}

void VolumeData::setValueAllVoxels(int val)
{
	assert(_voxelByteSize == 1 || _voxelByteSize == 2 || _voxelByteSize == 4 || _voxelByteSize == 8);

	memset(_data, val, _voxelByteSize * _dim.volume());
}


void VolumeData::castDataToFloat(float** dstPtr)
{
	SAFE_DELETE_ARRAY(*dstPtr);
	*dstPtr = new float[_dim.volume()];

	for (UINT64 k = 0; k < _dim.z; k++)
		for (UINT64 j = 0; j < _dim.y; j++)
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;
				float val = getVoxel(UINT64VECTOR3(i, j, k));
				(*dstPtr)[ind] = val;
			}
}

void VolumeData::setVoxel(UINT64VECTOR3 pos, float value, int channel) 
{
	assert(pos.x < _dim.x && pos.y < _dim.y && pos.z < _dim.z);
	assert(_voxelByteSize == 1 || _voxelByteSize == 2 || _voxelByteSize == 4);
	assert(channel >= 0 && channel < _numChannels);

	UINT64 ind = _numChannels * ((pos.z * _dim.y + pos.y) * _dim.x + pos.x) + channel;
	if (_voxelByteSize == 4)
	{
		if (_isFloat)
			reinterpret_cast<float*>(_data)[ind] = value;
		else
		{
			if (_isSigned)
				reinterpret_cast<int*>(_data)[ind] = int(value * MAX_VOL_VAL_INT);
			else
				reinterpret_cast<unsigned int*>(_data)[ind] = unsigned int(value * MAX_VOL_VAL_UINT);
		}
	}
	else if (_voxelByteSize == 2)
	{
		if (_isSigned)
			reinterpret_cast<short*>(_data)[ind] = short(value * MAX_VOL_VAL_SHORT);
		else
			reinterpret_cast<unsigned short*>(_data)[ind] = unsigned short(value* MAX_VOL_VAL_USHORT);
	}

	else if (_voxelByteSize == 1)
	{
		if (_isSigned)
			reinterpret_cast<char*>(_data)[ind] = char(value * MAX_VOL_VAL_CHAR);
		else
			reinterpret_cast<unsigned char*>(_data)[ind] = unsigned char(value * MAX_VOL_VAL_UCHAR);
	}
}

void VolumeData::setLayer(const void* buf, UINT64 z)
{
	assert(buf != NULL && z < _dim.z && _voxelByteSize <= 4);
	if( _voxelByteSize == 1 )
		memcpy_s( &(reinterpret_cast<unsigned char*>(_data)[z * _dim.x * _dim.y]), _dim.x * _dim.y * _voxelByteSize, buf, _dim.x * _dim.y * _voxelByteSize);
	else if(_voxelByteSize == 2)
	{
		memcpy_s( &(reinterpret_cast<unsigned short*>(_data)[z * _dim.x * _dim.y]), _dim.x * _dim.y * _voxelByteSize, buf, _dim.x * _dim.y * _voxelByteSize);
	}
	else
	{
		if(_isSigned)
		{
			if(_isFloat){
				// float
				memcpy_s( &(reinterpret_cast<float*>(_data)[z * _dim.x * _dim.y]), _dim.x * _dim.y * _voxelByteSize, buf, _dim.x * _dim.y * _voxelByteSize);
			}else{
				// int
				memcpy_s( &(reinterpret_cast<int*>(_data)[z * _dim.x * _dim.y]), _dim.x * _dim.y * _voxelByteSize, buf, _dim.x * _dim.y * _voxelByteSize);
			}
		}
		else
		{
			// uint
			memcpy_s( &(reinterpret_cast<unsigned int*>(_data)[z * _dim.x * _dim.y]), _dim.x * _dim.y * _voxelByteSize, buf, _dim.x * _dim.y * _voxelByteSize);

		}
	}
}

void VolumeData::getLayer(void* buf, UINT64 z) const
{
	assert( buf!= NULL && z < _dim.z && _voxelByteSize <= 4);
	if( _voxelByteSize == 1)
		memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<unsigned char*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
	else if(_voxelByteSize == 2)
	{
		if(_isSigned)
			memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<short*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
		else
			memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<unsigned short*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
	}
	else
	{
		if(_isSigned)
		{
			if(_isFloat){
				// float
				memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<float*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
			}else{
				// int
				memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<int*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
			}
		}
		else
		{
			// uint
			memcpy_s(buf, _dim.x * _dim.y * _voxelByteSize, &(reinterpret_cast<unsigned int*>(_data)[z * _dim.x * _dim.y]),  _dim.x * _dim.y * _voxelByteSize);
		}
	}
}

void VolumeData::getSlice(void* buf, UINT64 sliceId, int sliceDir) const
{
	assert( buf!= NULL && sliceDir <= 2);
	switch(sliceDir)
	{
	case 0: // get yz plane
		{
			if( _voxelByteSize == 1 )
			{
				// char...
				if( _isSigned )
				{
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 y = 0; y < _dim.y; y++ )
						{
							UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
							UINT64 sliceInd = z * _dim.y + y;
							reinterpret_cast<char*>(buf)[sliceInd] = reinterpret_cast<char*>(_data)[ind];
						}
					}
				}
				else
				{
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 y = 0; y < _dim.y; y++ )
						{
							UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
							UINT64 sliceInd = z * _dim.y + y;
							reinterpret_cast<unsigned char*>(buf)[sliceInd] = reinterpret_cast<unsigned char*>(_data)[ind];
						}
					}
				}
			}
			else if( _voxelByteSize == 2)
			{
				
				if( _isSigned )
				{// short
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 y = 0; y < _dim.y; y++ )
						{
							UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
							UINT64 sliceInd = z * _dim.y + y;
							reinterpret_cast<float*>(buf)[sliceInd] = reinterpret_cast<short*>(_data)[ind];
						}
					}
				}
				else
				{// unsigned short...
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 y = 0; y < _dim.y; y++ )
						{
							UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
							UINT64 sliceInd = z * _dim.y + y;
							reinterpret_cast<unsigned short*>(buf)[sliceInd] = reinterpret_cast<unsigned short*>(_data)[ind];
						}
					}
				}
				
			}
			else
			{
				if( _isSigned )
				{
					if( _isFloat )
					{
						// Float
						for(UINT64 z = 0; z < _dim.z; z++ )
						{
							for(UINT64 y = 0; y < _dim.y; y++ )
							{
								UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
								UINT64 sliceInd = z * _dim.y + y;
								reinterpret_cast<float*>(buf)[sliceInd] = reinterpret_cast<float*>(_data)[ind];
							}
						}
					}
					else
					{
						// Int
						for(UINT64 z = 0; z < _dim.z; z++ )
						{
							for(UINT64 y = 0; y < _dim.y; y++ )
							{
								UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
								UINT64 sliceInd = z * _dim.y + y;
								reinterpret_cast<int*>(buf)[sliceInd] = reinterpret_cast<int*>(_data)[ind];
							}
						}
					}
				}
				else
				{
					// unsigned int
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 y = 0; y < _dim.y; y++ )
						{
							UINT64 ind =  (z * _dim.y + y) * _dim.x + sliceId;
							UINT64 sliceInd = z * _dim.y + y;
							reinterpret_cast<unsigned int*>(buf)[sliceInd] = reinterpret_cast<unsigned int*>(_data)[ind];
						}
					}	
				}
			}
		}
		break;
	case 1: // get xz plane
		{
			if( _voxelByteSize == 1 )
			{
				// char...
				if( _isSigned )
				{
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 x = 0; x < _dim.x; x++ )
						{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
							reinterpret_cast<char*>(buf)[sliceInd] = reinterpret_cast<char*>(_data)[ind];
						}
					}
				}
				else
				{
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 x = 0; x < _dim.x; x++ )
						{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
							reinterpret_cast<unsigned char*>(buf)[sliceInd] = reinterpret_cast<unsigned char*>(_data)[ind];
						}
					}
				}
			}
			else if( _voxelByteSize == 2)
			{
				
				if( _isSigned )
				{// short
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 x = 0; x < _dim.x; x++ )
						{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
							reinterpret_cast<float*>(buf)[sliceInd] = reinterpret_cast<short*>(_data)[ind];
						}
					}
				}
				else
				{// unsigned short...
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 x = 0; x < _dim.x; x++ )
						{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
							reinterpret_cast<unsigned short*>(buf)[sliceInd] = reinterpret_cast<unsigned short*>(_data)[ind];
						}
					}
				}
				
			}
			else
			{
				if( _isSigned )
				{
					if( _isFloat )
					{
						// Float
						for(UINT64 z = 0; z < _dim.z; z++ )
						{
							for(UINT64 x = 0; x < _dim.x; x++ )
							{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
								reinterpret_cast<float*>(buf)[sliceInd] = reinterpret_cast<float*>(_data)[ind];
							}
						}
					}
					else
					{
						// Int
						for(UINT64 z = 0; z < _dim.z; z++ )
						{
							for(UINT64 x = 0; x < _dim.x; x++ )
							{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
								reinterpret_cast<int*>(buf)[sliceInd] = reinterpret_cast<int*>(_data)[ind];
							}
						}
					}
				}
				else
				{
					// unsigned int
					for(UINT64 z = 0; z < _dim.z; z++ )
					{
						for(UINT64 x = 0; x < _dim.x; x++ )
						{
							UINT64 ind = (z * _dim.y + sliceId) * _dim.x + x;
							UINT64 sliceInd = z * _dim.x + x;
							reinterpret_cast<unsigned int*>(buf)[sliceInd] = reinterpret_cast<unsigned int*>(_data)[ind];
						}
					}	
				}
			}
		}
		break;
	case 2: // get xy plane
		getLayer(buf, sliceId);
		break;
	}

}

void VolumeData::getSubVol(void* buf, UINT64VECTOR3 pos, UINT64VECTOR3 sz) const
{
	bool allZero = false;
	assert( buf != NULL ); // caller has to allocate the buffer
	for(UINT64 k = 0; k < sz.z; k++ )
		for(UINT64 j = 0; j< sz.y; j++ )
			for(UINT64 i = 0; i < sz.x; i++ )
			{
				unsigned int value = 0;
				if( (i + pos.x < _dim.x) && (j + pos.y < _dim.y) && (k + pos.z < _dim.z) )
					value = this->getVoxel(UINT64VECTOR3(i + pos.x, j + pos.y, k + pos.z));
				else
					value = 0;

				if( value != 0 )
					allZero = false;

				if (_voxelByteSize == 2)
					reinterpret_cast<unsigned short*>(buf)[ (k * sz.y + j) * sz.x + i ] = (unsigned short)value;
				else if(_voxelByteSize == 1) 
					reinterpret_cast<unsigned char*>(buf)[ (k * sz.y + j) * sz.x + i] = (unsigned char)value;
			}
	if( allZero )
		SAFE_DELETE_ARRAY(buf); // return a null buffer
}

void VolumeData::setSubVol(const void* buf, UINT64VECTOR3 pos, UINT64VECTOR3 sz)
{
	assert( buf != NULL ); // caller has to allocate the buffer
	for(UINT64 k = 0; k < sz.z; k++ )
		for(UINT64 j = 0; j< sz.y; j++ )
			for(UINT64 i = 0; i < sz.x; i++ )
			{
				unsigned int value = 0;
				if(_voxelByteSize == 2)
					value = float(((unsigned short*)buf)[ (k * sz.y + j) * sz.x + i ]); 
				else if( _voxelByteSize == 1)
					value = float(((unsigned char*)buf)[ (k * sz.y + j) * sz.x + i ]); 
				setVoxel(UINT64VECTOR3(i + pos.x, j + pos.y, k + pos.z), value);
			}
}

void VolumeData::computeHistogram1D()
{
	// simply traverse the volume 
	int histoLen = 1<<(_voxelByteSize * 8);
	_histogram1D.clear();
	_histogram1D.resize(histoLen, 0);
	UINT64 totalSize = _dim.volume();
	for(UINT64 i = 0; i < totalSize; i++ )
	{
		UINT64 destPos = 0;
		switch(_voxelByteSize)
		{
		case 1:
			destPos = reinterpret_cast<unsigned char*>(_data)[i];
		    break;
		case 2:
			destPos = reinterpret_cast<unsigned short*>(_data)[i];
			break;
		default:
			cout<<"voxel byte size not supported yet!"<<endl;
			break;
		}
		_histogram1D[destPos]++;
	}
}

void VolumeData::computeHistogram2D(VolumeData* gradMagVol)
{
	if( gradMagVol == NULL )
		gradMagVol = computeGradientMagnitude();

	int histoLen = 1<<(_voxelByteSize * 8);
	if( _voxelByteSize > 1 )
		histoLen = 1024;
	SAFE_DELETE_ARRAY(_histogram2D);
	_histogram2D = new UINT64[histoLen * histoLen];
	memset(_histogram2D, 0, sizeof(UINT64) * histoLen * histoLen);
	if( _voxelByteSize == 1 ){
		for(UINT64 i = 0; i < _dim.volume(); i++ )
		{
			float val = float(((unsigned char*)_data)[i]);
			float gm = float(((unsigned char*)gradMagVol->getDataBuf())[i]);

			int ix = val / 255.0f * (histoLen - 1);
			int iy = gm / 255.0f * (histoLen - 1);
			_histogram2D[UINT64(iy * histoLen + ix)]++;
		}
	}
	else if(_voxelByteSize == 2 )
	{
		for(UINT64 i = 0; i < _dim.volume(); i++ )
		{
			float val = float(((unsigned char*)_data)[i]);
			float gm = float(((unsigned char*)gradMagVol->getDataBuf())[i]);

			int ix = val / float(numeric_limits<unsigned short>::max()) * (histoLen - 1);
			int iy = gm / float(numeric_limits<unsigned short>::max()) * (histoLen - 1);
			_histogram2D[UINT64(iy * histoLen + ix)]++;
		}
	}

	// writeout scatterplot
	ofstream ofScatterplot("regularHisto2D.txt");
	for(UINT64 y = 0; y < histoLen; y++ )
	{
		for(UINT64 x = 0; x < histoLen; x++ )
			ofScatterplot<<_histogram2D[y * histoLen + x]<<" ";
		ofScatterplot<<endl;
	}
	ofScatterplot.close();
}

bool VolumeData::writeToFile(const char* filename, char* finalName)
{
	FILE *pFile; 

	std::string fname=std::string(filename);
	fname += std::string("-xyz-");
	fname += (number2String(_dim.x) + string("X")+ number2String(_dim.y) + string("X")+number2String(_dim.z));
	fname += std::string(".raw");

	pFile= fopen(fname.c_str(), "wb");
	if(finalName)
		memcpy(finalName, fname.c_str(), fname.size());
	if (NULL == pFile) {
		cout<<"can't open file!"<<endl;
		return false;
	}

	size_t size = fwrite(_data, UINT64(_voxelByteSize), _dim.volume() * _numChannels, pFile);
	fclose(pFile);
	if (size !=  _dim.volume()*_numChannels)
		return false;
	cout <<"write:"<<fname.c_str()<< endl;
	return true;
}

void VolumeData::writePPMFiles(const char* filename)
{
	Image img(int(_dim.x), int(_dim.y));
	for(UINT64 z = 0; z < _dim.z; z++ )
	{
		for(UINT64 y = 0; y < _dim.y; y++ )
			for(UINT64 x = 0; x < _dim.x; x++ )
			{
				double val = 0.0;
				if( _voxelByteSize == 1)
					val = (double)reinterpret_cast<unsigned char*>(_data)[(z * _dim.y + y) * _dim.x + x];
				else if( _voxelByteSize == 2)
					val = (double)reinterpret_cast<unsigned short*>(_data)[(z * _dim.y + y) * _dim.x + x];
				else if( _voxelByteSize == 4)
					val = (double)reinterpret_cast<unsigned short*>(_data)[(z * _dim.y + y) * _dim.x + x];
				val /= 255.0;
				Color C;
				C.set_r(val);
				C.set_g(val);
				C.set_b(val);
				img.set(int(x), int(y), C);
			}
			string sliceName = string(filename) + "-" + number2String(z) + ".ppm";
			img.write(sliceName);
	}
}

bool VolumeData::writeToNrrdFile(const char* filename)
{
	if (_numChannels > 1)
	{
		cout << "Failed to write data. NRRD supports only data with one channel!" << endl;
		return false;
	}

	Nrrd* nval = nrrdNew();

	size_t v_size[3];
	v_size[0] = _dim.x;
	v_size[1] = _dim.y;
	v_size[2] = _dim.z;

	int dataType = nrrdTypeUChar;
	switch (_voxelByteSize)
	{
	case 1:
		if (_isSigned)
			dataType = nrrdTypeChar;
		else
			dataType = nrrdTypeUChar;
		break;
	case 2:
		if (_isSigned)
			dataType = nrrdTypeShort;
		else
			dataType = nrrdTypeUShort;
		break;
	case 4:
		if (_isFloat)
		{
			dataType = nrrdTypeFloat;
		}
		else
		{
			if (_isSigned)
				dataType = nrrdTypeInt;
			else
				dataType = nrrdTypeUInt;
		}
		break;
	case 8:
		if (_isFloat)
			dataType = nrrdTypeDouble;
		else
		{
			if (_isSigned)
				dataType = nrrdTypeLLong;
			else
				dataType = nrrdTypeULLong;
		}
		break;
	default:
		break;
	}

	////////////////////
	// WARNING: fSize must be EXPLICITYLY converted to const int*, otherwise, we got weird errors!!!!
	nrrdWrap_nva(nval, getDataBuf(), dataType, 3, (const size_t*)v_size);
	//nrrdWrap(nval, _data, dataType, 3, fSize[0], fSize[1], fSize[2]);
	if (nrrdSave(filename, nval, NULL))
	{
		char* err = biffGetDone(NRRD);
		fprintf(stderr, "%s\n", err);
		free(err);

		nrrdNix(nval);
		return false;
	}
	nrrdNix(nval);
	return true;
}

void VolumeData::writeLayerToPPMFile(const char* filename, UINT64 z)
{
	Image img(int(_dim.x), int(_dim.y));
	cout << "Min max = " << _minval << ", " << _maxval << endl;
	for (UINT64 y = 0; y < _dim.y; y++) 
	{
		for (UINT64 x = 0; x < _dim.x; x++)
		{
			double val = 0.0;
			if (!_isSigned) {
				if (_voxelByteSize == 1) {
					val = (double)reinterpret_cast<unsigned char*>(_data)[(z * _dim.y + y) * _dim.x + x];
					val /= 255.0;
				}
				else if (_voxelByteSize == 2)
					val = (double)reinterpret_cast<unsigned short*>(_data)[(z * _dim.y + y) * _dim.x + x];
				else if (_voxelByteSize == 4)
					val = (double)reinterpret_cast<unsigned int*>(_data)[(z * _dim.y + y) * _dim.x + x];
			}
			else
			{
				if (_voxelByteSize == 1)
				{
					val = (double)reinterpret_cast<char*>(_data)[(z * _dim.y + y) * _dim.x + x];
					val = (val - 127.0) / 255.0;
				}
				else if (_voxelByteSize == 2)
					val = (double)reinterpret_cast<short*>(_data)[(z * _dim.y + y) * _dim.x + x];
				else if (_voxelByteSize == 4)
				{
					if (_isFloat)
						val = (double)reinterpret_cast<float*>(_data)[(z * _dim.y + y) * _dim.x + x];
					else
						val = (double)reinterpret_cast<int*>(_data)[(z * _dim.y + y) * _dim.x + x];
				}
			}
			val = (val - _minval)/(_maxval - _minval);
			Color C;
			C.set_r(val);
			C.set_g(val);
			C.set_b(val);
			img.set(int(x), int(y), C);
		}
	}
	string sliceName = string(filename) + "-" + number2String(z) + ".ppm";
	img.write(sliceName);
}

VolumeData* VolumeData::computeGradientMagnitude()
{
	VolumeData* vol = 0;
	if (_voxelByteSize != 1 && _voxelByteSize != 2 && _voxelByteSize != 4){
		cout << "Unsupported format!" << endl;
		return NULL; // data format not supported
	}
	// Always allocate a temproray volume of float first
	void* fTemp = malloc(_dim.volume() * sizeof(float));
	memset(fTemp, 0, _dim.volume() * sizeof(float));


	float maxgradmag = -numeric_limits<float>::max();
	for(UINT64 k = 0; k < _dim.z; k++ )
	{
		for(UINT64 j = 0; j < _dim.y; j++ )
		{
			for(UINT64 i = 0; i < _dim.x; i++ )
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;
				UINT64VECTOR3 px = UINT64VECTOR3(MIN(_dim.x-1, i+1), j, k);
				UINT64VECTOR3 nx = UINT64VECTOR3(MAX(0, int(i)-1), j, k);

				UINT64VECTOR3 py = UINT64VECTOR3(i, MIN(_dim.y-1,j+1), k);
				UINT64VECTOR3 ny = UINT64VECTOR3(i, MAX(0, int(j)-1), k);

				UINT64VECTOR3 pz = UINT64VECTOR3(i, j, MIN(_dim.z-1, k+1));
				UINT64VECTOR3 nz = UINT64VECTOR3(i, j, MAX(0, int(k)-1));

				FLOATVECTOR3 grad = FLOATVECTOR3(float(getVoxel(px)) - float(getVoxel(nx)), 
					float(getVoxel(py)) - float(getVoxel(ny)),
					float(getVoxel(pz)) - float(getVoxel(nz)));

				float gradMag = grad.length();
				maxgradmag = MAX(gradMag, maxgradmag);
				((float*)fTemp)[ind] = gradMag;
			}
		}
	}

	if (maxgradmag == 0.0f){
		free(fTemp);
		return NULL;
	}
	else
	{
		void* temp = malloc(_dim.volume() * UINT64(_voxelByteSize));
		memset(temp, 0, _dim.volume() * UINT64(_voxelByteSize));
		// Do normalization of gradmag
		for (UINT64 i = 0; i < _dim.volume(); i++)
		{
			float normGradMag = ((float*)fTemp)[i] / maxgradmag; 

			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)temp)[i] = char(normGradMag * 128.0f);
				else
					((unsigned char*)temp)[i] = unsigned char(normGradMag * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)temp)[i] = short(normGradMag * numeric_limits<short>::max());
				else
					((unsigned short*)temp)[i] = unsigned short(normGradMag * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)temp)[i] = float(normGradMag);
				else
				{
					if (_isSigned)
						((int*)temp)[i] = int(normGradMag * numeric_limits<int>::max());
					else
						((unsigned int*)temp)[i] = int(normGradMag * numeric_limits<unsigned int>::max());
				}
			}
		}

		vol = new VolumeData(_dim, _voxelByteSize, temp);
		free(temp);
		free(fTemp);
		return vol;
	}
}



VolumeData* VolumeData::computeSecondDerivMag()
{
	if (_voxelByteSize != 1 && _voxelByteSize != 2 && _voxelByteSize != 4){
		cout << "Unsupported format!" << endl;
		return NULL; // data format not supported
	}
	VolumeData* vol = 0;
	void* temp = malloc(_dim.volume() * UINT64(_voxelByteSize));
	memset(temp, 0, _dim.volume() * UINT64(_voxelByteSize));

	// Always allocate a temproray volume of float first
	float* fTemp = new float[_dim.volume()];
	memset(fTemp, 0, _dim.volume() * sizeof(float));

	float max2ndDerivMag = -numeric_limits<float>::max();
	for (UINT64 k = 0; k < _dim.z; k++)
	{
		for (UINT64 j = 0; j < _dim.y; j++)
		{
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;

				UINT64VECTOR3 c = UINT64VECTOR3(i, j, k); // The center

				UINT64VECTOR3 px = UINT64VECTOR3(MIN(_dim.x - 1, i + 1), j, k);
				UINT64VECTOR3 nx = UINT64VECTOR3(MAX(0, double(i) - 1), j, k);

				UINT64VECTOR3 py = UINT64VECTOR3(i, MIN(_dim.y - 1, j + 1), k);
				UINT64VECTOR3 ny = UINT64VECTOR3(i, MAX(0, double(j) - 1), k);

				UINT64VECTOR3 pz = UINT64VECTOR3(i, j, MIN(_dim.z - 1, k + 1));
				UINT64VECTOR3 nz = UINT64VECTOR3(i, j, MAX(0, double(k) - 1));

				// Compute Laplacian
				float lapX = -2.0f * float(getVoxel(c)) + float(getVoxel(px)) + float(getVoxel(nx));
				float lapY = -2.0f * float(getVoxel(c)) + float(getVoxel(py)) + float(getVoxel(ny));
				float lapZ = -2.0f * float(getVoxel(c)) + float(getVoxel(pz)) + float(getVoxel(nz));

				FLOATVECTOR3 lap(lapX, lapY, lapZ);
				float lapMag = lap.length();
				max2ndDerivMag = MAX(max2ndDerivMag, lapMag);
				fTemp[ind] = lapMag;
			}
		}
	}
	cout << "Max 2nd derivative magnitude = " << max2ndDerivMag << endl;
	// Do normalization
	if (max2ndDerivMag == 0.0f){
		SAFE_DELETE_ARRAY(fTemp); 
		return NULL;
	}
	else
	{
		float oneOverMaxMag = 1.0f / max2ndDerivMag;
		void* temp = malloc(_dim.volume() * UINT64(_voxelByteSize));
		memset(temp, 0, _dim.volume() * UINT64(_voxelByteSize));
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			float normGradMag = fTemp[ind] * oneOverMaxMag;
			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)temp)[ind] = char(normGradMag * 128.0f);
				else
					((unsigned char*)temp)[ind] = unsigned char(normGradMag * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)temp)[ind] = short(normGradMag * numeric_limits<short>::max());
				else
					((unsigned short*)temp)[ind] = unsigned short(normGradMag * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)temp)[ind] = float(normGradMag);
				else
				{
					if (_isSigned)
						((int*)temp)[ind] = int(normGradMag * numeric_limits<int>::max());
					else
						((unsigned int*)temp)[ind] = int(normGradMag * numeric_limits<unsigned int>::max());
				}
			}
		}
		vol = new VolumeData(_dim, _voxelByteSize, temp);
		free(temp);
		free(fTemp);
		return vol;
	}
}

void VolumeData::computeGradMagandDirSecondDeriv(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool  normalize)
{
	if (_voxelByteSize != 1 && _voxelByteSize != 2 && _voxelByteSize != 4){
		cout << "Unsupported format!" << endl;
		return; // data format not supported
	}
	// Always allocate a temproray volume of float first
	float* fGradMag = new float[_dim.volume()];
	memset(fGradMag, 0, _dim.volume() * sizeof(float));
	float maxgradmag = -numeric_limits<float>::max();
	// Store the actual gradient volume
	Vector3f* fGrad = new Vector3f[_dim.volume()];
	memset(fGrad, 0, _dim.volume() * sizeof(Vector3f));
	// First, compute gradient volume (vec3!)
	for (UINT64 k = 0; k < _dim.z; k++)
	{
		for (UINT64 j = 0; j < _dim.y; j++)
		{
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;

				UINT64VECTOR3 c = UINT64VECTOR3(i, j, k); // The center

				UINT64VECTOR3 px = UINT64VECTOR3(MIN(_dim.x - 1, i + 1), j, k);
				UINT64VECTOR3 nx = UINT64VECTOR3(MAX(0, double(i) - 1), j, k);

				UINT64VECTOR3 py = UINT64VECTOR3(i, MIN(_dim.y - 1, j + 1), k);
				UINT64VECTOR3 ny = UINT64VECTOR3(i, MAX(0, double(j) - 1), k);

				UINT64VECTOR3 pz = UINT64VECTOR3(i, j, MIN(_dim.z - 1, k + 1));
				UINT64VECTOR3 nz = UINT64VECTOR3(i, j, MAX(0, double(k) - 1));
				// Compute gradient
				Vector3f grad = Vector3f(float(getVoxel(px)) - float(getVoxel(nx)),
					float(getVoxel(py)) - float(getVoxel(ny)),
					float(getVoxel(pz)) - float(getVoxel(nz)));

				fGrad[ind] = grad;

				float gradMag = grad.norm(); // Second norm hopefully
				maxgradmag = MAX(gradMag, maxgradmag);
				fGradMag[ind] = gradMag;
			}
		}
	}
	cout << "Max gradient magnitude = " << maxgradmag << endl;
	// Second, compute directional second-derivative
	// Store the second derivative magnitude
	float maxSecDeriv = -numeric_limits<float>::max();
	float minSecDeriv = numeric_limits<float>::max();
	float* fSecDeriv = new float[_dim.volume()];
	memset(fSecDeriv, 0, _dim.volume() * sizeof(float));

	///////////////////
	// Test plane:
	PlaneCoeffs P_test;
	// Plane 1
	P_test.c1 = 6.0f;// -0.11f;
	P_test.c2 = 5.0f;// 1.89f;
	P_test.c3 = 8.0f;// -0.78f;
	P_test.c0 = 8.0f; //2.40f;
	// Plane 2
	//P_test.c1 = -0.11f;
	//P_test.c2 =  1.89f;
	//P_test.c3 = -0.78f;
	//P_test.c0 =  2.40f;

	for (UINT64 k = 0; k < _dim.z; k++)
	{
		for (UINT64 j = 0; j < _dim.y; j++)
		{
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;
				// integer indices
				UINT64 ipx = (k * _dim.y + j) * _dim.x + MIN(_dim.x - 1, i + 1);
				UINT64 inx = (k * _dim.y + j) * _dim.x + MAX(0, double(i) - 1);

				UINT64 ipy = (k * _dim.y + MIN(_dim.y - 1, j + 1)) * _dim.x + i;
				UINT64 iny = (k * _dim.y + MAX(0, double(j) - 1)) * _dim.x + i;

				UINT64 ipz = (MIN(_dim.z - 1, k + 1) * _dim.y + j) * _dim.x + i;
				UINT64 inz = (MAX(0, double(k) - 1) * _dim.y + j) * _dim.x + i;

				// Get the gradient at center
				Vector3f g = fGrad[ind];
				float g_len2 = g.squaredNorm();
				if (g_len2 == 0.0f){
					fSecDeriv[ind] = 0.0f;
					continue;
				}
				// compute Hessian matrix
				Matrix3f H;
				H.row(0) = fGrad[ipx] - fGrad[inx]; // ddf/ddx, ddf/(dxdy), ddf/(dxdz)
				H.row(1) = fGrad[ipy] - fGrad[iny]; // ddf/(dxdy), ddf/ddy, ddf/(dydz)
				H.row(2) = fGrad[ipz] - fGrad[inz]; // ddf(dxdz), ddf/(dydz), ddf/ddz

				float secDeriv = g.transpose() * H * g;
				secDeriv /= g_len2;

				///////////////////////////////////////////////
				//// TODO: Disable comments to test samples on a plane
				//// TEST: set value based on the z value
				//secDeriv = float(int(k/128));
				//secDeriv = 0.5f; // all values lie on the same plane
				//secDeriv = 0.5f;
				// //Test: all values are on a plane in the 3D data domain
				//secDeriv = computeDataSampleOnPlane(P_test, getVoxel(UINT64VECTOR3(i,j,k)), g.norm());

				fSecDeriv[ind] = secDeriv;


				maxSecDeriv = MAX(secDeriv, maxSecDeriv);
				minSecDeriv = MIN(secDeriv, minSecDeriv);
			}
		}
	}

	cout << "Max 2nd derivative = " << maxSecDeriv << " Min 2nd derivative = " << minSecDeriv << endl;
	// Do normalization and quantization
	void* tempGM = 0;
	if (maxgradmag != 0.0f)
	{
		// No normalization
		float oneOverMaxMag = 1.0f;
		// With normalization
		//if ( normalize )
		oneOverMaxMag = 1.0f / maxgradmag;
		tempGM = malloc(_dim.volume() * UINT64(_voxelByteSize));
		memset(tempGM, 0, _dim.volume() * UINT64(_voxelByteSize));
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			float normGradMag = fGradMag[ind] * oneOverMaxMag;
			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempGM)[ind] = char(normGradMag * 128.0f);
				else
					((unsigned char*)tempGM)[ind] = unsigned char(normGradMag * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempGM)[ind] = short(normGradMag * /*numeric_limits<short>::max()*/MAX_VOL_VAL_SHORT);
				else
					((unsigned short*)tempGM)[ind] = unsigned short(normGradMag * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempGM)[ind] = float(normGradMag);
				else
				{
					if (_isSigned)
						((int*)tempGM)[ind] = int(normGradMag * numeric_limits<int>::max());
					else
						((unsigned int*)tempGM)[ind] = int(normGradMag * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volGradMag = new VolumeData(_dim, _voxelByteSize, tempGM);
	}
	// For second derivative
	void* tempSD = 0;
	tempSD = malloc(_dim.volume() * UINT64(_voxelByteSize));
	memset(tempSD, 0, _dim.volume() * UINT64(_voxelByteSize));
	//if (minSecDeriv != maxSecDeriv)
	{

		float oneOverRange = 1.0f;
		//if (normalize)
		{
			if (minSecDeriv != maxSecDeriv)
				oneOverRange /= (maxSecDeriv - minSecDeriv);
		}
		float normSecDeriv = 0.5f;
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			if (minSecDeriv != maxSecDeriv) // do normalization only when the min max are not the same
				normSecDeriv = (fSecDeriv[ind] - minSecDeriv) * oneOverRange;

			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempSD)[ind] = char(normSecDeriv * 128.0f);
				else
					((unsigned char*)tempSD)[ind] = unsigned char(normSecDeriv * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempSD)[ind] = short(normSecDeriv * numeric_limits<short>::max());
				else
					((unsigned short*)tempSD)[ind] = unsigned short(normSecDeriv * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempSD)[ind] = float(normSecDeriv);
				else
				{
					if (_isSigned)
						((int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<int>::max());
					else
						((unsigned int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volDirSecDeriv = new VolumeData(_dim, _voxelByteSize, tempSD);
	}
	free(tempGM);
	free(tempSD);
	SAFE_DELETE_ARRAY(fGrad);
	SAFE_DELETE_ARRAY(fGradMag);
	SAFE_DELETE_ARRAY(fSecDeriv);

}


void VolumeData::genSyntheticVols(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool  normalize)
{
	if (_voxelByteSize != 1 && _voxelByteSize != 2 && _voxelByteSize != 4){
		cout << "Unsupported format!" << endl;
		return; // data format not supported
	}
	// Always allocate a temproray volume of float first
	float* fGradMag = new float[_dim.volume()];
	memset(fGradMag, 0, _dim.volume() * sizeof(float));
	float maxgradmag = -numeric_limits<float>::max();
	// Second, compute directional second-derivative
	// Store the second derivative magnitude
	float maxSecDeriv = -numeric_limits<float>::max();
	float minSecDeriv = numeric_limits<float>::max();
	float* fSecDeriv = new float[_dim.volume()];
	memset(fSecDeriv, 0, _dim.volume() * sizeof(float));
	// Store the actual gradient volume
	Vector3f* fGrad = new Vector3f[_dim.volume()];
	memset(fGrad, 0, _dim.volume() * sizeof(Vector3f));
	///////////////////
	// Test plane:
	PlaneCoeffs P_test;
	// Plane 1
	P_test.c1 = 8.0f;// -0.11f;
	P_test.c2 = 8.0f;// 1.89f;
	P_test.c3 = 8.0f;// -0.78f;
	P_test.c0 = 8.0f; //2.40f;
	// Plane 2
	//P_test.c1 = -0.11f;
	//P_test.c2 =  1.89f;
	//P_test.c3 = -0.78f;
	//P_test.c0 =  2.40f;
	// First, compute gradient volume (vec3!)
	for (UINT64 k = 0; k < _dim.z; k++)
	{
		for (UINT64 j = 0; j < _dim.y; j++)
		{
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;

				UINT64VECTOR3 c = UINT64VECTOR3(i, j, k); // The center
				//float scalarV = float(getVoxel(c));
				//float gradMag = grad.norm(); // Second norm hopefully
				float secDeriv = 0.0f;
				float scalarV = g_randGen.rand();
				float gradMag = g_randGen.rand();
				// 1. Make v and gradmag to be positively correlated
				float a = g_randGen.rand();
				float b = g_randGen.rand();
				gradMag = scalarV * a + b;


				// 2. Generate samples on the plane

				secDeriv = computeDataSampleOnPlane(P_test, scalarV, gradMag);
				//rotateDataSampleToPlane(P_test, scalarV, gradMag, secDeriv);
				setVoxel(c, scalarV);
				maxgradmag = MAX(gradMag, maxgradmag);
				fGradMag[ind] = gradMag;
				fSecDeriv[ind] = secDeriv;
			

				maxSecDeriv = MAX(secDeriv, maxSecDeriv);
				minSecDeriv = MIN(secDeriv, minSecDeriv);
			}
		}
	}
	cout << "Max gradient magnitude = " << maxgradmag << endl;


	cout << "Max 2nd derivative = " << maxSecDeriv<<" Min 2nd derivative = "<< minSecDeriv << endl;
	// Do normalization and quantization
	void* tempGM = 0;
	if ( maxgradmag != 0.0f)	
	{
		// No normalization
		float oneOverMaxMag = 1.0f;
		// With normalization
		//if ( normalize )
			oneOverMaxMag = 1.0f / maxgradmag;
		tempGM = malloc(_dim.volume() * UINT64(_voxelByteSize));
		memset(tempGM, 0, _dim.volume() * UINT64(_voxelByteSize));
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			float normGradMag = fGradMag[ind] * oneOverMaxMag;
			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempGM)[ind] = char(normGradMag * 128.0f);
				else
					((unsigned char*)tempGM)[ind] = unsigned char(normGradMag * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempGM)[ind] = short(normGradMag * numeric_limits<short>::max());
				else
					((unsigned short*)tempGM)[ind] = unsigned short(normGradMag * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempGM)[ind] = float(normGradMag);
				else
				{
					if (_isSigned)
						((int*)tempGM)[ind] = int(normGradMag * numeric_limits<int>::max());
					else
						((unsigned int*)tempGM)[ind] = int(normGradMag * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volGradMag = new VolumeData(_dim, _voxelByteSize, tempGM);
	}
	// For second derivative
	void* tempSD = 0;
	tempSD = malloc(_dim.volume() * UINT64(_voxelByteSize));
	memset(tempSD, 0, _dim.volume() * UINT64(_voxelByteSize));
	//if (minSecDeriv != maxSecDeriv)
	{
		
		float oneOverRange = 1.0f;
		//if (normalize)
		{
			if (minSecDeriv != maxSecDeriv)
				oneOverRange /= (maxSecDeriv - minSecDeriv);
		}
		float normSecDeriv = 0.5f;
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			if ( minSecDeriv != maxSecDeriv ) // do normalization only when the min max are not the same
				normSecDeriv = (fSecDeriv[ind] - minSecDeriv) * oneOverRange;

			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempSD)[ind] = char(normSecDeriv * 128.0f);
				else
					((unsigned char*)tempSD)[ind] = unsigned char(normSecDeriv * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempSD)[ind] = short(normSecDeriv * numeric_limits<short>::max());
				else
					((unsigned short*)tempSD)[ind] = unsigned short(normSecDeriv * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempSD)[ind] = float(normSecDeriv);
				else
				{
					if (_isSigned)
						((int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<int>::max());
					else
						((unsigned int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volDirSecDeriv = new VolumeData(_dim, _voxelByteSize, tempSD);
	}
	free(tempGM);
	free(tempSD);
	SAFE_DELETE_ARRAY(fGrad);
	SAFE_DELETE_ARRAY(fGradMag);
	SAFE_DELETE_ARRAY(fSecDeriv); 

}

void VolumeData::genSynGoodPatternVols(VolumeData** volGradMag, VolumeData** volDirSecDeriv, bool  normalize)
{
	if (_voxelByteSize != 1 && _voxelByteSize != 2 && _voxelByteSize != 4){
		cout << "Unsupported format!" << endl;
		return; // data format not supported
	}
	// Always allocate a temproray volume of float first
	float* fGradMag = new float[_dim.volume()];
	memset(fGradMag, 0, _dim.volume() * sizeof(float));
	float maxgradmag = -numeric_limits<float>::max();
	// Second, compute directional second-derivative
	// Store the second derivative magnitude
	float maxSecDeriv = -numeric_limits<float>::max();
	float minSecDeriv = numeric_limits<float>::max();
	float* fSecDeriv = new float[_dim.volume()];
	memset(fSecDeriv, 0, _dim.volume() * sizeof(float));
	// Store the actual gradient volume
	Vector3f* fGrad = new Vector3f[_dim.volume()];
	memset(fGrad, 0, _dim.volume() * sizeof(Vector3f));
	///////////////////
	// Test plane:
	PlaneCoeffs P_test;
	// Plane 1
	P_test.c1 = 8.0f;// -0.11f;
	P_test.c2 = 8.0f;// 1.89f;
	P_test.c3 = 8.0f;// -0.78f;
	P_test.c0 = 8.0f; //2.40f;
	// Plane 2
	//P_test.c1 = -0.11f;
	//P_test.c2 =  1.89f;
	//P_test.c3 = -0.78f;
	//P_test.c0 =  2.40f;
	// First, compute gradient volume (vec3!)

	bool isNegCorr = false;
	float a = -0.6f; // negative corr								
	float b = 1.0f;
	float a2 = 0.4f;
	float b2 = 0.0f;
	if (!isNegCorr)
	{
		a = 0.4f;
		b = 0.f;
	}

	for (UINT64 k = 0; k < _dim.z; k++)
	{
		for (UINT64 j = 0; j < _dim.y; j++)
		{
			for (UINT64 i = 0; i < _dim.x; i++)
			{
				UINT64 ind = (k * _dim.y + j) * _dim.x + i;

				UINT64VECTOR3 c = UINT64VECTOR3(i, j, k); // The center

				//float scalarV = float(getVoxel(c));
				//float gradMag = grad.norm(); // Second norm hopefully
				float secDeriv = secDeriv = g_randGen.rand53();// g_randGen.rand();
				float scalarV = g_randGen.rand53();
				float gradMag = 0.0f;// g_randGen.rand();

				
				// 2. Generate samples on the plane
				
				//secDeriv = computeDataSampleOnPlane(P_test, scalarV, gradMag);
				//rotateDataSampleToPlane(P_test, scalarV, gradMag, secDeriv);
				setVoxel(c, scalarV);

				// 1. Make v and gradmag to be positively correlated in the centeral area in data domain
				// 2. Generate samples on the plane
				if (scalarV > 0.05f && scalarV < 0.95f)
				{

					gradMag = scalarV * a + b;
				}
				else
				{
					gradMag = g_randGen.rand53();
				
				}
				
				if (scalarV > 0.4f && scalarV < 0.9f)
					secDeriv = computeDataSampleOnPlane(P_test, scalarV, gradMag);
				maxgradmag = MAX(gradMag, maxgradmag);
				fGradMag[ind] = gradMag;
				fSecDeriv[ind] = secDeriv;


				maxSecDeriv = MAX(secDeriv, maxSecDeriv);
				minSecDeriv = MIN(secDeriv, minSecDeriv);
			}
		}
	}
	cout << "Max gradient magnitude = " << maxgradmag << endl;


	cout << "Max 2nd derivative = " << maxSecDeriv << " Min 2nd derivative = " << minSecDeriv << endl;
	// Do normalization and quantization
	void* tempGM = 0;
	if (maxgradmag != 0.0f)
	{
		// No normalization
		float oneOverMaxMag = 1.0f;
		// With normalization
		//if ( normalize )
		oneOverMaxMag = 1.0f / maxgradmag;
		tempGM = malloc(_dim.volume() * UINT64(_voxelByteSize));
		memset(tempGM, 0, _dim.volume() * UINT64(_voxelByteSize));
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			float normGradMag = fGradMag[ind] * oneOverMaxMag;
			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempGM)[ind] = char(normGradMag * 128.0f);
				else
					((unsigned char*)tempGM)[ind] = unsigned char(normGradMag * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempGM)[ind] = short(normGradMag * numeric_limits<short>::max());
				else
					((unsigned short*)tempGM)[ind] = unsigned short(normGradMag * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempGM)[ind] = float(normGradMag);
				else
				{
					if (_isSigned)
						((int*)tempGM)[ind] = int(normGradMag * numeric_limits<int>::max());
					else
						((unsigned int*)tempGM)[ind] = int(normGradMag * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volGradMag = new VolumeData(_dim, _voxelByteSize, tempGM);
	}
	// For second derivative
	void* tempSD = 0;
	tempSD = malloc(_dim.volume() * UINT64(_voxelByteSize));
	memset(tempSD, 0, _dim.volume() * UINT64(_voxelByteSize));
	//if (minSecDeriv != maxSecDeriv)
	{

		float oneOverRange = 1.0f;
		//if (normalize)
		{
			if (minSecDeriv != maxSecDeriv)
				oneOverRange /= (maxSecDeriv - minSecDeriv);
		}
		float normSecDeriv = 0.5f;
		for (UINT64 ind = 0; ind < _dim.volume(); ind++)
		{
			if (minSecDeriv != maxSecDeriv) // do normalization only when the min max are not the same
				normSecDeriv = (fSecDeriv[ind] - minSecDeriv) * oneOverRange;

			if (_voxelByteSize == 1)
			{
				if (_isSigned)
					((char*)tempSD)[ind] = char(normSecDeriv * 128.0f);
				else
					((unsigned char*)tempSD)[ind] = unsigned char(normSecDeriv * 255.0f);
			}
			else if (_voxelByteSize == 2)
			{
				if (_isSigned)
					((short*)tempSD)[ind] = short(normSecDeriv * numeric_limits<short>::max());
				else
					((unsigned short*)tempSD)[ind] = unsigned short(normSecDeriv * numeric_limits<unsigned short>::max());
			}
			else if (_voxelByteSize == 4)
			{
				if (_isFloat)
					((float*)tempSD)[ind] = float(normSecDeriv);
				else
				{
					if (_isSigned)
						((int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<int>::max());
					else
						((unsigned int*)tempSD)[ind] = int(normSecDeriv * numeric_limits<unsigned int>::max());
				}
			}
		}
		*volDirSecDeriv = new VolumeData(_dim, _voxelByteSize, tempSD);
	}
	free(tempGM);
	free(tempSD);
	SAFE_DELETE_ARRAY(fGrad);
	SAFE_DELETE_ARRAY(fGradMag);
	SAFE_DELETE_ARRAY(fSecDeriv);

}


VolumeData* VolumeData::downSampleOneLevel()
{
	VolumeData* vol = 0;

	UINT64 ddimX = _dim.x >> 1; 
	UINT64 ddimY = _dim.y >> 1;
	UINT64 ddimZ = _dim.z >> 1;
	float* tempf = (float*)malloc(ddimX * ddimY * ddimZ * sizeof(float)); //for more precision
	memset(tempf, 0, ddimX * ddimY * ddimZ * sizeof(float));
	void* temp = malloc(ddimX * ddimY * ddimZ * UINT64(_voxelByteSize));
	memset(temp, 0, ddimX * ddimY * ddimZ * UINT64(_voxelByteSize));

	if(_voxelByteSize == 1)
	{
		unsigned char* dataChar = (unsigned char*)_data;
		for(UINT64 i=0; i<_dim.x; i++)
			for(UINT64 j=0; j<_dim.y; j++)
				for(UINT64 k=0; k<_dim.z; k++)
				{								
					tempf[k/2*ddimX*ddimY+ j/2*ddimY+ i/2]
					+= dataChar[k*_dim.x*_dim.y + j*_dim.y + i];
				}	
		for(UINT64 i=0; i<ddimX*ddimY*ddimZ; i++)
		{
				((unsigned char*)(temp))[i] = (unsigned char)(tempf[i]/8.0);
		}
	}
	else if(_voxelByteSize == 2)
	{
		unsigned short* dataShort = (unsigned short*)_data;
		for(UINT64 i=0; i<_dim.x; i++)
			for(UINT64 j=0; j<_dim.y; j++)
				for(UINT64 k=0; k<_dim.z; k++)
				{								
					tempf[k/2*ddimX*ddimY+ j/2*ddimY+ i/2]
					+= dataShort[k*_dim.x*_dim.y + j*_dim.y + i];
				}	
		for(UINT64 i=0; i<ddimX*ddimY*ddimZ; i++)
		{
				((unsigned short*)(temp))[i] = (unsigned short)(tempf[i]/8.0);
		}
	}

	free(tempf);

	vol = new VolumeData(UINT64VECTOR3(ddimX, ddimY, ddimZ), _voxelByteSize, temp);

	return vol;
}

void VolumeData::convert2Float(float** ppFloatVol, bool normalize) const
{

	// Convert to float vol
	SAFE_DELETE_ARRAY(*ppFloatVol);
	if( _isSigned || _isFloat )
		return; // don't even bother for signed or float data!
	*ppFloatVol = new float[_dim.volume()];
	float UCharNormalizer = 1.0f / 255.0f;
	float UShortNormalizer = 1.0f / float(numeric_limits<unsigned short>::max());
	if(_voxelByteSize == 1)
	{
		for(UINT64 i = 0; i < _dim.volume(); i++ ){
			(*ppFloatVol)[i] = (float)reinterpret_cast<unsigned char*>(_data)[i];
			if( normalize )
				(*ppFloatVol)[i] *= UCharNormalizer;
		}
	}
	else if(_voxelByteSize == 2)
	{
		for(UINT64 i = 0; i < _dim.volume(); i++ ){
			(*ppFloatVol)[i] = (float)reinterpret_cast<unsigned short*>(_data)[i];
			if( normalize )
				(*ppFloatVol)[i] *= UShortNormalizer;
		}
	}
}

int VolumeData::findBackgrndVal()
{
	if(_histogram1D.empty())
	{
		computeHistogram1D();
	}
	UINT64 maxCnt = 0;
	int    backgroundVal = -1;
	for(size_t i = 0; i < _histogram1D.size(); i++)
	{
		if( maxCnt < _histogram1D[i] ){
			maxCnt = _histogram1D[i];
			backgroundVal = int(i);
		}
	}
	// We just take background values at two ends of value range
	if( backgroundVal >= 5 && backgroundVal <= 250 )
		backgroundVal = -1;

	_backgrndVal = backgroundVal;
	return backgroundVal;
}

float VolumeData::computeDataSampleOnPlane(const PlaneCoeffs& P, const float& x, const float& y)
{
	float z = (P.c0 - P.c1 * x - P.c2 * y) / P.c3;
	return z;
}

void VolumeData::rotateDataSampleToPlane(const PlaneCoeffs& P, float& x, float& y, float& z)
{
	 // Assume originally, x, y, z sit on a plane with normal pointing to z+

	FLOATVECTOR3  newnorm(P.c1, P.c2, P.c3);
	newnorm.normalize();
	FLOATVECTOR3  oldnorm(0, 0, 1);
	FLOATVECTOR3 v = newnorm % oldnorm; // v = a x b cross product
	float s = v.length();
	float c = newnorm ^ oldnorm; // v = a.b dot product
	FLOATMATRIX3 R; // R = I + [v] + v^2(1-c)/(s^2)
	FLOATMATRIX3 Vx(
		0.0f, -v.z, v.y,
		v.z, 0.0f, -v.x,
		-v.y, v.x, 0.0f);
	FLOATMATRIX3 R3 = Vx * Vx * (1.0f - c) / (s * s);
	R = R + Vx + R3;
	// Transform point
	FLOATVECTOR3 oldPt(x, y, z);
	FLOATVECTOR3 newPt = R * oldPt;
	x = newPt.x;
	y = newPt.y;
	z = newPt.z;
}


bool VolumeData::loadFromFile(const char* filename)
{

	std::string fname = std::string(filename);
	std::string extName = SysTools::ToLowerCase(SysTools::GetExt(fname));
	if (extName != "raw")
	{
		if (extName == "nrrd" || extName == "nhdr")
		{
			cout << "Loading nrrd file..." << endl;
			return loadFromNrrdFile(filename, true, g_params.IsNormalizeData());
		}
		else {
			cout << "File format is not supported!" << endl;
			return false;
		}
	}
	// parse file name
	size_t lXyz = fname.find("-xyz-");
	size_t lDimX = lXyz + string("-xyz-").length();
	size_t lDimY = fname.find("X");
	size_t lDimZ = fname.find("X", lDimY + 1);
	size_t lEnd = fname.find(".", lDimZ + 1);

	FILE* pFile;
	pFile = fopen(filename, "rb");
	if (NULL == pFile) {
		cout << "Can't open file: " << filename << "!" << endl;
		return false;
	}

	_dim.x = atoi(fname.substr(lDimX, lDimY - lDimX).c_str());
	_dim.y = atoi(fname.substr(lDimY + 1, lDimZ - lDimY).c_str());
	_dim.z = atoi(fname.substr(lDimZ + 1, lEnd - lDimZ).c_str());
	_voxelByteSize = sizeof(float); // as float by default
	size_t volSize = fread(_data, UINT64(_voxelByteSize), _dim.volume(), pFile);
	fclose(pFile);
	if (volSize != _dim.volume())
		return false;
	cout << "Loaded: " << fname.c_str() << endl;
	return true;
}


bool VolumeData::loadFromNrrdFile(const char* filename, bool convertToFloat, bool doNormalize)
{
	int xres, yres, zres;
	Nrrd* nin;
	NrrdIoState* nio;
	char* err;
	//create a nrrd; at this point this is just an empty container 
	nin = nrrdNew();
	nio = nrrdIoStateNew();

	if (nrrdIoStateSet(nio, nrrdIoStateSkipData, AIR_TRUE)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "%s: trouble reading \"%s\":\n%s", __FUNCTION__, filename, err);
		free(err);
		return false;
	}

	//res = nrrdLoad(nin, filename, nio);
	//if (res)
	//{
	//	fprintf(stderr, "cannot open %s.\n", filename );
	//	return false;
	//}

	/* read in the nrrd from file */
	if (nrrdLoad(nin, filename, nio)) {
		err = biffGetDone(NRRD);
		fprintf(stderr, "%s: trouble reading \"%s\":\n%s", __FUNCTION__, filename, err);
		free(err);
		return false;
	}


	xres = nin->axis[0].size;
	yres = nin->axis[1].size;
	if (nin->axis[2].size == 0)
		zres = 1;
	else
		zres = nin->axis[2].size;

	fprintf(stderr, "(xres, yres, zres)=(%d, %d, %d)\n", xres, yres, zres);

	int byteSize = 1;
	unsigned char* tmp_array = NULL;
	bool isSigned = false;
	bool isFloat = false;
	switch (nin->type)
	{
	case nrrdTypeUChar:
		byteSize = 1;
		isSigned = false;
		isFloat = false;
		break;
	case nrrdTypeShort:
		byteSize = 2;
		isSigned = true;
		isFloat = false;
		break;
	case nrrdTypeUShort:
		byteSize = 2;
		isSigned = false;
		isFloat = false;
		break;
	case nrrdTypeFloat:
		byteSize = 4;
		isSigned = true;
		isFloat = true;
		break;
	case nrrdTypeInt:
		byteSize = 4;
		isSigned = true;
		isFloat = true;
		break;
	case nrrdTypeDouble:
		byteSize = 8;
		isSigned = true;
		isFloat = true;
		break;
	default:
		fprintf(stderr, "%s data type unsupported.\n", filename);
		return false;

	}
	_dim = UINT64VECTOR3(xres, yres, zres);
	_voxelByteSize = byteSize;
	_isFloat = isFloat;
	_isSigned = isSigned;

	UINT64 dataSize = UINT64(_voxelByteSize) * _dim.volume();
	if (_data)
		SAFE_DELETE_ARRAY(_data);

	_data = new unsigned char[dataSize];

	nin->data = _data;
	// read again to get actual data
	if (nrrdLoad(nin, filename, NULL))
	{
		fprintf(stderr, "%s appears to be NRRD, but corrupted", filename);
		delete[] tmp_array;
		nrrdNix(nin);
		return false;
	}
	//// Cast to float
	if (convertToFloat)
	{
		if (nin->type != nrrdTypeFloat)
		{
			float* fData = NULL;
			castDataToFloat(&fData);
			SAFE_DELETE_ARRAY(_data);
			_data = fData;
			_voxelByteSize = 4;
		}
		_isFloat = true;
	}

	// Need to normalize if it's float
	float minVal = 0.0f;
	float maxVal = 0.0f;
	//if (nin->type != nrrdTypeUChar && nin->type != nrrdTypeChar)
	{
		if (doNormalize)
		{
			normalize(minVal, maxVal);
			_minval = 0.0f;
			_maxval = 1.0f;
		}
		else
			findMinMax(_minval, _maxval);
	}
	//else
	//{
	//	if (doNormalize) {
	//		normalize(minVal, maxVal);
	//		_minval = 0.0f;
	//		_maxval = 1.0f;
	//	}
	//}

	return true;
}


void VolumeData::findMinMax(float& minVal, float& maxVal)
{
	minVal = numeric_limits<float>::max();
	maxVal = -numeric_limits<float>::max();
	// First pass find min,max
	for (UINT64 z = 0; z < _dim.z; z++)
	{
		for (UINT64 y = 0; y < _dim.y; y++)
		{
			for (UINT64 x = 0; x < _dim.x; x++)
			{
				float voxVal = getVoxel(UINT64VECTOR3(x, y, z));
				if (isnan(voxVal)) { // Empty value!
				//	setVoxel(UINT64VECTOR3(x, y, z), -10.0f);
					continue;
				}
				minVal = MIN(minVal, voxVal);
				maxVal = MAX(maxVal, voxVal);
			}
		}
	}
	cout << "Vol val min = " << minVal << ", max = " << maxVal << endl;
	this->_minval = minVal;
	this->_maxval = maxVal;
}


void VolumeData::normalizeToRange(const float minVal, const float maxVal)
{
	// Second pass, normalization
	float oneOverRange = (minVal == maxVal) ? 0.0f : 1.0f / (maxVal - minVal);
	for (UINT64 z = 0; z < _dim.z; z++)
	{
		for (UINT64 y = 0; y < _dim.y; y++)
		{
			for (UINT64 x = 0; x < _dim.x; x++)
			{
				float voxVal = getVoxel(UINT64VECTOR3(x, y, z));
				float normVal = (voxVal - minVal) * oneOverRange;
				setVoxel(UINT64VECTOR3(x, y, z), normVal);
			}
		}
	}
}