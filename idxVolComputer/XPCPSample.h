#ifndef XPCPSAMPLE_H_
#define XPCPSAMPLE_H_
#include "MyVectors.h"
#include "pFlatPt.h"
#include <vector>
#include <iostream>

// !!!! NOTE
// The 1-flat and 2-flat use Inselberg's PCP coordinates, not normalized!!!
// 
// Structure for extended PCP sample with 1-flat and 2-flat

//-----------------------------------------------------------
// PlainData for XPCP
//|0 -- N-1 | N -- 2 * (N-1) + N | 2N-1 -- 2N-1 + 2 * (N-2) |
//|ND value | 1-flat             | 2-flat                   |
//

enum XPCP_TYPE
{
	XPCP_REG = 0,	// regular certain xpcp
	XPCP_UNCERTAIN  // uncertain xpcp
};

class XPCPSample // for N-D data
{
public:
	XPCPSample():
		data_type(XPCP_REG),
		num_pts(1)
	{};

	XPCPSample(int rawDataDim):
		data_type(XPCP_REG),
		num_pts(1)
	{
		x = std::vector<float>(rawDataDim, 0);
		pcp_1flat = std::vector<FLOATVECTOR2>(rawDataDim - 1);
		pcp_2flat = std::vector<FLOATVECTOR2>(rawDataDim - 2);
		strength_1flat = std::vector<FLOATVECTOR2>(rawDataDim - 1); // strength of 1-flats
		strength_2flat = std::vector<FLOATVECTOR2>(rawDataDim - 2); // strength of 2-flats
		isRotated_1flat = std::vector<bool>(rawDataDim - 1, false);
	};

	XPCPSample(const XPCPSample& other):
		data_type(XPCP_REG)
	{
		x = other.x;
		pcp_1flat = other.pcp_1flat;
		pcp_2flat = other.pcp_2flat;
		strength_1flat = other.strength_1flat;
		strength_2flat = other.strength_2flat;
		num_pts = other.num_pts;
		isRotated_1flat = other.isRotated_1flat;
	};

	inline void toStdVector(std::vector<float>& plainData)
	{
		// Since 2-flat could be empty in some cases, we still have to make ALL dataitem the same size, therefore, we calculate the data length using pcp_1flat sizes
		plainData.resize(x.size() + pcp_1flat.size() * 2 + (pcp_1flat.size() - 1) * 2, 0.0f);
		size_t cnt = 0;
		std::copy(x.begin(), x.end(), plainData.begin());
		cnt += x.size(); 
		for (size_t i = 0; i < pcp_1flat.size(); i++){
			plainData[(cnt + 2 * i) ] = pcp_1flat[i].x;
			plainData[(cnt + 2 * i) + 1] = pcp_1flat[i].y;
		}
		cnt += pcp_1flat.size() * 2; 
		for (size_t i = 0; i < pcp_2flat.size(); i++){
			plainData[(cnt + 2 * i)] = pcp_2flat[i].x;
			plainData[(cnt + 2 * i) + 1] = pcp_2flat[i].y;
		}
	}

	inline void to4DVec(std::vector<FLOATVECTOR4>& v4Data, int pFlat, UINT64 ptId)
	{
		// Record:
		// |  ptId | 2D coords in PCP (Inselberg's def) | id: sub2D space |
		// the 2D coords of the flats and the ID of which sub 2D space they come from
		size_t N = x.size();
		v4Data.clear();
		switch (pFlat) // p-flat can only be 1,2,3...
		{
		case 1:
			v4Data.resize(N - 1);
			for (size_t i = 0; i < pcp_1flat.size(); i++)
			{
				v4Data[i] = FLOATVECTOR4(float(ptId), pcp_1flat[i].x, pcp_1flat[i].y, float(i));
			}
			break;
		case 2:
			v4Data.resize(N - 2);
			for (size_t i = 0; i < pcp_2flat.size(); i++)
			{
				v4Data[i] = FLOATVECTOR4(float(ptId), pcp_2flat[i].x, pcp_2flat[i].y, float(i));
			}
			break;
		default:
			break; // not supported yet
		}
	}

	
	inline void toStdVector(std::vector<float>& plainData, int component) // Get the #component from the tuple: component = 0: x, 1: pcp_1flat, 2: pcp_2flat
	{
		size_t N = x.size();
		plainData.clear();
		switch (component)
		{
		case 0:
			plainData.resize(N, 0.0f);
			std::copy(x.begin(), x.end(), plainData.begin());
			break;
		case 1:
			plainData.resize(2 * (N - 1), 0.0f);
			for (size_t i = 0; i < pcp_1flat.size(); i++){
				plainData[2 * i] = pcp_1flat[i].x;
				plainData[2 * i + 1] = pcp_1flat[i].y;
			}
			break;
		case 2:
			plainData.resize(2 * (N - 2), 0.0f);
			for (size_t i = 0; i < pcp_2flat.size(); i++){
				plainData[2 * i] = pcp_2flat[i].x;
				plainData[2 * i + 1] = pcp_2flat[i].y;
			}
			break;
		default:
			plainData.resize(N, 0.0f);
			std::copy(x.begin(), x.end(), plainData.begin());
			break;
		}
	}

	// Convert to a list of pFlatPt
	virtual void to_pFlatRec(std::vector<pFlatPt*>& pFlatRec, int pFlat, UINT64 ptId, bool repeat_pcp = false)
	{
		pFlatRec.clear();
		switch (pFlat)
		{
		case 1:
			for (size_t i = 0; i < pcp_1flat.size(); i++)
			{
				int j = int(i / num_pts);
				int subdimId = (num_pts == 1)? i : int(i) % num_pts;
		
				pFlatRec.push_back(new pFlatPt(pcp_1flat[i].x, pcp_1flat[i].y, subdimId,
					ptId, strength_1flat[j].x, strength_1flat[j].y, isRotated_1flat[i]));
			}
			break;
		case 2:
			for (size_t i = 0; i < pcp_2flat.size(); i++)
			{
				int j = i / num_pts;
				int subdimId = (num_pts == 1) ? i : int(i) % num_pts;
				if (repeat_pcp){
					subdimId /= 4;
					j /= 4;
				}

				pFlatRec.push_back(new pFlatPt(pcp_2flat[i].x, pcp_2flat[i].y, subdimId, 
					ptId, strength_2flat[j].x, strength_2flat[j].y));
			}
			break;
		default:
			break;

		}
	}

	virtual XPCP_TYPE getXPCPtype() const
	{
		return data_type;
	}
public:
	std::vector<float>     x; // N data values
	std::vector<FLOATVECTOR2>     pcp_1flat; // N-1 1-flat data points in PCP
	std::vector<FLOATVECTOR2>     pcp_2flat; // (N-2) 2-flat data points in PCP. NOTE: we omit the second index point for a 2-flat to save some space!!! The actual requirement fo r
	//==========================================================
	// Also records the strength of each p-flat
	// for each p-flat indexed point, record (1) some measurement compared to its full space, (2) the relative ranking inside its subspace
	std::vector<FLOATVECTOR2> strength_1flat;
	std::vector<FLOATVECTOR2> strength_2flat;
	XPCP_TYPE				  data_type;
	int                       num_pts;
	std::vector<bool>         isRotated_1flat; // Used only in Nguyen-Rosen's case
};

#endif