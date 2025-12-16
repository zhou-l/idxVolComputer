#ifndef UNCERTAIN_XPCP_SAMPLE_H
#define UNCERTAIN_XPCP_SAMPLE_H

#include "XPCPSample.h"


class UncertainXPCPSample : public XPCPSample
{
public:
	UncertainXPCPSample():
		XPCPSample()
	{
		data_type = XPCP_UNCERTAIN;
	}

	UncertainXPCPSample(int rawDataDim, int n_pts)
	{
		data_type = XPCP_UNCERTAIN;
		x = std::vector<float>(rawDataDim, 0);
		pcp_1flat = std::vector<FLOATVECTOR2>(n_pts* (rawDataDim - 1));
		pcp_2flat = std::vector<FLOATVECTOR2>(n_pts* (rawDataDim - 2));
		strength_1flat = std::vector<FLOATVECTOR2>(rawDataDim - 1);
		strength_2flat = std::vector<FLOATVECTOR2>(rawDataDim - 2);
		num_pts = n_pts;
	}

	UncertainXPCPSample(const UncertainXPCPSample& other)
	{
		data_type = XPCP_UNCERTAIN;
		x = other.x;
		pcp_1flat = other.pcp_1flat;
		pcp_2flat = other.pcp_2flat;
		strength_1flat = other.strength_1flat;
		strength_2flat = other.strength_2flat;
		num_pts = other.num_pts;

	}

	virtual ~UncertainXPCPSample()
	{}

	// Convert to a list of pFlatPt
  virtual void to_pFlatRec(std::vector<pFlatPt*>& pFlatRec, int pFlat, UINT64 ptId, bool repeat_pcp = false)
	{
	  pFlatRec.clear();
	  int dim1D = int(pcp_1flat.size()) / num_pts;
	  int dim2D = int(pcp_2flat.size()) / num_pts;

		switch (pFlat)
		{
		case 1:
			for (size_t i = 0; i < pcp_1flat.size(); i++)
			{
				int k = i % dim1D;
				pFlatRec.push_back(new pFlatPt(pcp_1flat[i].x, pcp_1flat[i].y, k, 
					ptId, strength_1flat[k].x, strength_1flat[k].y));
			}
			break;
		case 2:
			for (size_t i = 0; i < pcp_2flat.size(); i++)
			{
				int k = i;
				if (repeat_pcp)
					k /= 4;
				k = k % dim2D;
				pFlatRec.push_back(new pFlatPt(pcp_2flat[i].x, pcp_2flat[i].y, int(k), ptId, strength_2flat[k].x, strength_2flat[k].y));
			}
			break;
		default:
			break;

		}
	}

};

#endif