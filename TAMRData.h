#include <iostream>                                                             
#include <algorithm>                                                            
#include <limits>                                                               
#include <array>                                                                
#include <unordered_map>                                                        
#include "common/sg/importer/Importer.h"
#include "common/sg/transferFunction/TransferFunction.h"
#include "common/sg/common/Common.h"
#include "common/sg/geometry/Spheres.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"



using namespace ospray;
using namespace ospray::sg;

#define MAXLEVELWIDTH 64


namespace ospray{
  namespace tamr{

  struct TAMRVoxel{
    vec3f lower;                                                              
    int level;    

    size_t indexInBuffer; 

    inline float getVoxelWidth() const{
      float width = (float)(1 << level); 
      return width / MAXLEVELWIDTH; 
    }

    inline vec3f getCenter() const{
      return lower + vec3f(0.5 * getVoxelWidth());
    }
  };


  struct TAMRLevel{
    float cellWidth;
    float rcpCellWidth;
    float halfCellWidth;
    int level;

    box3f bounds; 
  };

  
  struct TAMRData{
    vec3f amrOrigin; 

    float scale; 

    std::vector<TAMRLevel> levels; 

    std::vector<TAMRVoxel> voxels;
  };


  }// ::ospray::tamr
}// ::ospray
