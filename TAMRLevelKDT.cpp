#include "TAMRLevelKDT.h"

namespace ospray {
  namespace tamr {

    TAMRLevelKDT::TAMRLevelKDT(const TAMRData &input, int level) 
    {
      std::unordered_map<int,TAMRLevel>::const_iterator itr = input.voxelsInLevel.find(level);

      if(itr == input.voxelsInLevel.end()){
        throw std::runtime_error("An wrong AMR level is specified");
      }

      const TAMRLevel&  levelInput = itr->second;

      this->level.cellWidthInModel = levelInput.cellWidthInModel;
      this->level.cellWidth = levelInput.cellWidth;
      this->level.halfCellWidth = levelInput.halfCellWidth;
      this->level.rcpCellWidth = levelInput.rcpCellWidth;

      this->worldBounds = levelInput.bounds;

      PRINT(levelInput.voxels.size());

      node.resize(1);
      buildRec(0,worldBounds, levelInput.voxels);

    }

    TAMRLevelKDT::~TAMRLevelKDT(){
      for(auto &l : leaf){
        l.voxels.clear();
      }

      leaf.clear();
      node.clear();
    }

    void TAMRLevelKDT::makeLeaf(index_t nodeID,
                    const box3f &bounds,
                    const std::vector<TAMRVoxel> &voxels)
    {
      node[nodeID].dim = 3; 
      node[nodeID].ofs = this->leaf.size();
      node[nodeID].numItems = voxels.size();

      TAMRLevelKDT::Leaf newLeaf; 
      newLeaf.bounds = bounds;
      newLeaf.voxels = voxels;

      this->leaf.push_back(newLeaf);

    }

    void TAMRLevelKDT::makeInner(index_t nodeID, int dim, float pos, int childID)
    {
      node[nodeID].dim = dim;
      node[nodeID].pos = pos;
      node[nodeID].ofs = childID;  
    }

    void TAMRLevelKDT::buildRec(int nodeID, 
                                const box3f &bounds, 
                                std::vector<TAMRVoxel> voxels)
    {

      if(voxels.empty())
        return;

      int bestDim = -1;

      vec3f bs = bounds.size() + vec3f(1);
      bestDim = (bs.x >= bs.y) ? ((bs.x >= bs.z) ? 0 : 2) : ((bs.y >= bs.z) ? 1 : 2); 

      size_t numSlot = (size_t)bs.product();
      if(numSlot == voxels.size())
      {
        makeLeaf(nodeID, bounds, voxels);
      }else{
#if 0
        float bestPos = bounds.lower[bestDim] + 0.5 * bs[bestDim] - 1;
#else 
       float bestPos = getBestPos(bounds,voxels,bestDim);
#endif 



        std::vector<TAMRVoxel> l, r;
        box3f lBounds, rBounds;
        for(auto & voxel : voxels){
          if(voxel.lower[bestDim] <= bestPos)
          {
            l.push_back(voxel);
            lBounds.extend(voxel.lower);
          }
          else{
            r.push_back(voxel);
            rBounds.extend(voxel.lower);
          }
        }

        // PRINT(bounds);
        // PRINT(bestDim);
        // PRINT(bestPos);
        // PRINT(lBounds);
        // PRINT(rBounds);
        // PRINT("====================");

        int newNodeID = node.size();
        makeInner(nodeID, bestDim, bestPos, newNodeID);

        node.push_back(TAMRLevelKDT::Node());
        node.push_back(TAMRLevelKDT::Node());

        buildRec(newNodeID+0,lBounds,l);
        buildRec(newNodeID+1,rBounds,r);
      }
    }


    float TAMRLevelKDT::getBestPos(const box3f& bounds, std::vector<TAMRVoxel> voxels, int dim)
    {

      // count the point number in each grid point on specific dimension.
      // e.g. bounds = [[0,0,0],[3,3,2]], dim = 0
      std::unordered_map<float, int> pNumInDim; 
      for(auto& voxel: voxels){
        pNumInDim[voxel.lower[dim]]++;
      }

      vec3f bs = bounds.size() + vec3f(1);
      int filledNum(1);
      for(int i = 0 ; i < 3; i++){
        if(i != dim)
          filledNum *= bs[i];
      }



      std::vector<float> potentialPos;
      float start = bounds.lower[dim];
      float end = bounds.upper[dim];
      while(start < end){
        if(pNumInDim[start] != pNumInDim[start +1]){
          potentialPos.push_back(start);
        }
        start++;
      }
 
      assert(potentialPos.size()!=0);
 
      // if(filledNum == 713){
      //   PRINT(filledNum);
      //   PRINT(bounds);
      //   PRINT(dim);
      //   PRINT(potentialPos.size());
      //   PRINT("====================");

      // }

      float bestPos = std::numeric_limits<float>::infinity();
      float mid = bounds.center()[dim];

      if(potentialPos.size() == 0)
        bestPos = mid;

      for (const auto &split : potentialPos) {
          if (fabsf(split - mid) < fabsf(bestPos-mid))
            bestPos = split;
      }

      return bestPos;
    } 

  }  // namespace tamr
}  // namespace ospray
