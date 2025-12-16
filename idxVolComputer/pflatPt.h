#ifndef PFLATPT_H
#define PFLATPT_H
#include "StdDefines.h"
#include <iostream>
// representation of p-flat points
// special case for p == 0
class zeroFlatPt
{
public:
	zeroFlatPt():
		_clusterId(0)
	{}
	~zeroFlatPt(){ _pos.clear(); }

	zeroFlatPt(const std::vector<float>& pos, int clusterId = 0)
	{
		_pos = pos;
		_clusterId = clusterId;
	}

	void set(const std::vector<float>& pos, int clusterId = 0)
	{
		_pos = pos;
		_clusterId = clusterId;
	}
	// Must implement coord() and oneMinusPercentile
	const float* coord() const { return (const float*)&(_pos[0]); }
	float oneMinusPercentile() const { return 1.0f; }// Always pass the test in KD tree searching
	const std::vector<float>& val() const { return _pos; }

	const int clusterId() const { return _clusterId; }
private:
	std::vector<float>  _pos; // ND coordinates of the p-flat : IN Inselberg's PCP coordinate!!!
	int                 _clusterId; // id of which the cluster this point belongs
};
// General case for p >= 1
class pFlatPt{
public:
	pFlatPt()
	{
		_pos[0] = _pos[1] = 0.0f;
		_rawPtId = 0;
		_subDimId = 0;
		_strength = 0.0f;
		_oneMinusPercentile = 1.0f; // full percentile
		_isrotated = false;
	}

	pFlatPt(float x, float y, int subDimId, UINT64 rawPtId, float strength = 0.0f, float oneMinusPercentile = 0.0f, bool isRotated = false)
	{
		_pos[0] = x;
		_pos[1] = y;
		_subDimId = subDimId;
		_rawPtId = rawPtId;
		_strength = strength;
		_oneMinusPercentile = oneMinusPercentile;
		_isrotated = isRotated;
	}
	void set(float x, float y, int subDimId, UINT64 rawPtId, float strength = 0.0f, float oneMinusPercentile = 0.0f, bool isRotated = false)
	{
		_pos[0] = x;
		_pos[1] = y;
		_subDimId = subDimId;
		_rawPtId = rawPtId;
		_strength = strength;
		_oneMinusPercentile = oneMinusPercentile;
		_isrotated = isRotated;
	}
	// accessors: must implement coord() for the KD tree 
	const float* coord() const { return _pos; }
	float x() const { return _pos[0]; }
	float y() const { return _pos[1]; }
	
	UINT64 rawId() const { return _rawPtId; }
	void   setRawId(UINT64 id) { _rawPtId = id; }
	
	int    subDimId() const { return _subDimId; }
	void   setSubDimId(int subDimId) { _subDimId = subDimId; }
	
	float strength() const { return _strength; }
	float oneMinusPercentile() const { return _oneMinusPercentile; }
	
	void  setIsRotated(bool rotate) { _isrotated = rotate; }
	bool  isRotated() const { return _isrotated; }

	friend std::ostream& operator<<(std::ostream& out, const pFlatPt& pt)
	{
		return out << pt.coord()[0] << " " << pt.coord()[1] << " " << pt.subDimId() <<" "<< pt.rawId();
	}

	friend std::ostream& operator<<(std::ostream& out, pFlatPt* ppt)
	{
		return out << ppt->coord()[0] << " " << ppt->coord()[1] << " " << ppt->subDimId() << " " << ppt->rawId();
	}
private:
	float  _pos[2]; // 2D coordinates of the p-flat : IN Inselberg's PCP coordinate!!!
	UINT64 _rawPtId; // id of the original point that generates the p-flat
	int    _subDimId; // which sub-dimensions this p-flat sits in 
	float  _strength; // the length of the norm resulting in this p-flat!
	float  _oneMinusPercentile; // 100% - oneMinusPercentile of this p-flat idx pt inside its subspace
	bool   _isrotated; // ONLY for Nguyen-Rosen tvcg case! Indicates whether the point is rotated! 
};


#endif