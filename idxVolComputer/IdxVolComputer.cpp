#include "IdxVolComputer.h"
#include "global.h"
#include "PointMonteCarloSampler.h"
#include "PointDataHandler.h"
#include "XPCPdataAnalyzer.h"

IdxVolComputer::IdxVolComputer()
{
}

IdxVolComputer::~IdxVolComputer()
{
}

void IdxVolComputer::processSpatialDomainContinuousIdxPt(const std::vector<VolumeData*>& volData)
{
    cout << "Processing continuous indexed points in the spatial domain" << endl;
    if (volData.empty()) {
        cout << "VolData List is empty! Skip spatial domain processing!" << endl;
        return;
    }
    //////////////////////////////////////
    //Output a slice of each volume for importance sampling test in Matlab
    for (size_t i = 0; i < volData.size(); i++)
    {
        string filename = "VolData" + number2String(i);
        volData[i]->writeLayerToPPMFile(filename.c_str(), volData[i]->getDim().z / 2);
    }
    //////////////////////////////////////
    // get spatial samples first
    vector<vector<float>> sampleLocs;
    //float sampleRatio = 0.002f;
    // Dimension of the indexed points volume
    UINT64VECTOR3 idxVolDim = volData[0]->getDim();
    // Sample over a grid or do monte carlo sampling?
    bool sampleOnGrid = true;
    PointMonteCarloSampler* sampler = dynamic_cast<PointMonteCarloSampler*>(_dataHandler);
    if (sampleOnGrid)
    { // sample on a regular grid

        float gridSampleRate = g_params.gridSampleRate();
        idxVolDim = UINT64VECTOR3(UINT64(gridSampleRate * float(idxVolDim.x)),
            UINT64(gridSampleRate * float(idxVolDim.y)),
            UINT64(gridSampleRate * float(idxVolDim.z)));
        // get sample locations and sample values
        sampleLocs.resize(idxVolDim.volume());
        g_sampled_raw.resize(idxVolDim.volume());
        cout << "g_sampled_raw size = " << g_sampled_raw.size() << endl;
        g_samplePortion = gridSampleRate * gridSampleRate * gridSampleRate;
        cout << "g_sampledPortion = " << g_samplePortion << endl;
        vector<float> pos(3, 0.0f);
        for (UINT64 z = 0; z < idxVolDim.z; z++)
        {
            for (UINT64 y = 0; y < idxVolDim.y; y++)
            {
                for (UINT64 x = 0; x < idxVolDim.x; x++)
                {
                    UINT64 idx = (z * idxVolDim.y + y) * idxVolDim.x + x;

                    pos[0] = float(x) / float(idxVolDim.x - 1);
                    pos[1] = float(y) / float(idxVolDim.y - 1);
                    pos[2] = float(z) / float(idxVolDim.z - 1);

                    sampleLocs[idx] = pos;

                    g_sampled_raw[idx] = sampler->sampleVoxelWithInterpNormCoords(FLOATVECTOR3(pos[0], pos[1], pos[2]));
                }
            }
        }
    }
    else
        //Monte carlo sampling
        sampleLocs = sampler->sampleSpatialDomain(g_samplePortion, g_sampled_raw/*_sampledData*/);
    // This should be made a global variable
    int num_neighbors = g_params.SpatialNeighborNum();
    // Should be set in the configuration dialog
    FLOATVECTOR3 R = g_params.SpatialNeighborSize();//g_spatialSampleNeighborSize; // FLOATVECTOR3(5.0f, 5.0f, 5.0f);
    cout << "R = " << R << endl;
    bool writeFile = true;

    // NOTE: this should be made a global variable!!!
    //bool doCompareValDomain = true;
    if (g_processFullSpace)
    {
        g_params.setDataVols(volData);
        // To compare with the value domain computation, we need another pass
        if (g_params.IsCompareValueDomain()) {
            _xpcpAnalyzer->processRawData(g_sampled_raw, g_numNN_valDomain);
            // And, record the value domain results in the variables for comparison
            g_xpcpDataS2 = g_xpcpData;
            g_1flat_listS2 = g_1flat_list;

            // Set volume data set 2
            vector<VolumeData*> idxPtVolsS2;
            saveIdxPtDataToVolumes(idxPtVolsS2, idxVolDim, int(volData.size()), sampleLocs, g_xpcpDataS2);

            g_params.setIdxPtVolsS2(idxPtVolsS2);
        }
        // Compute indexed points in the spatial domain
        _xpcpAnalyzer->analyzeSpatialRawDataContIdxPts(volData, sampleLocs, idxVolDim, num_neighbors, R, writeFile);

    }
    else
        _xpcpAnalyzer->processRawDataInSubSpaces(g_sampled_raw, g_numNN_valDomain);


    vector<int> labels;
    if (g_params.doClustering())
    {
        XPCPdataAnalyzer::cluster(g_sampled_raw, labels, g_params.numClusters());
        g_params.clusterLabelData(labels);
    }
}

bool IdxVolComputer::generateIdxPtVols(const QStringList&, bool negative)
{
    return false;
}

bool IdxVolComputer::loadVolumeFile(const QString& fileName, VolumeData** vol, float& minVal, float& maxVal, bool doNormalize)
{
	return false;
}

void IdxVolComputer::saveIdxPtDataToVolumes(std::vector<VolumeData*>& idxPtVols, UINT64VECTOR3 idxVolDim, int nVols, const std::vector<md_val>& sampleLocs, const std::vector<XPCPSample>& xpcpData)
{
}
