#ifndef TAMRDATA_H_
#define TAMRDATA_H_

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <unordered_map>
#include "common/sg/common/Common.h"
#include "common/sg/geometry/Spheres.h"
#include "common/sg/importer/Importer.h"
#include "common/sg/transferFunction/TransferFunction.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"

using namespace ospray;
using namespace ospray::sg;

namespace ospray {
  namespace tamr {

    struct TAMRVoxel
    {
      vec3f lower;
      int level;
      size_t indexInBuffer;
    };

    struct TAMRLevel
    {
      float cellWidthInModel;
      float cellWidth;
      float rcpCellWidth;
      float halfCellWidth;
      int level;

      box3f bounds;
      std::vector<TAMRVoxel> voxels;

      inline void push_voxel(TAMRVoxel &voxel)
      {
        voxels.push_back(voxel);
        bounds.extend(voxel.lower);
      }
    };

    struct TAMRData
    {
      vec3f amrOrigin;
      float cellScale;
      std::unordered_map<int, TAMRLevel> voxelsInLevel;
    };

  }  // namespace tamr
}  // namespace ospray


#endif