#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "common/sg/common/Common.h"
#include "common/sg/geometry/Spheres.h"
#include "common/sg/importer/Importer.h"
#include "common/sg/transferFunction/TransferFunction.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"

#include "ospcommon/AffineSpace.h"
#include "ospcommon/box.h"
#include "ospcommon/math.h"
#include "ospcommon/memory/malloc.h"
#include "ospcommon/range.h"
#include "ospcommon/xml/XML.h"
#include "ospray/ospray.h"

#include "TAMRData.h"
#include "TAMRLevelKDT.h"

using namespace ospcommon;
using namespace ospray;
using namespace ospray::sg;

using namespace ospray::tamr;

static std::vector<float> tfn_colors = {
      1.0f, 0.28f, 0.28f, 0.86f, 
      1.0f, 0.f,   0.f,   0.36f, 
      1.0f, 0.f,   1.f,   1.f, 
      1.0f, 0.f,   0.5f,  0.f, 
      1.0f, 1.f,   1.f,   0.f, 
      1.0f, 1.f,   0.38f, 0.f, 
      1.0f, 0.42f, 0.f,   0.f, 
      1.0f, 0.88f, 0.3f,  0.3f 
};

inline vec3f findColorForValue(std::vector<float> colors, int i, int N, bool isRandom)
{
  vec3f c; 
  int colorNum = colors.size() / 4;
  if(N <= colorNum){
    c = vec3f(colors[i * 4 + 1], colors[i * 4 + 2],colors[i * 4 + 3]);
  }else{
    float a = (1.0f/(N -1)) * i  * (colorNum - 1);
    int lo = (int)std::floor(a);
    int up = lo + 1;
    float reminder = a - lo;
    if(isRandom){
      // int rd(0);
      // while(rd == lo)
      //   rd = random() % colors.size();

      // c = (1.0f - reminder) * vec3f(colors[rd * 4 + 1], colors[rd * 4 + 2],colors[rd * 4 + 3])
      //       + reminder * vec3f(colors[(rd +1) * 4 + 1], colors[(rd +1) * 4 + 2],colors[(rd +1) * 4 + 3]);

      int ci = i % colorNum;
      c= vec3f(colors[ci * 4 + 1], colors[ci * 4 + 2],colors[ci * 4 + 3]);
    }
    else{
      c = (1.0f - reminder) * vec3f(colors[lo * 4 + 1], colors[lo * 4 + 2],colors[lo * 4 + 3])
            + reminder * vec3f(colors[up * 4 + 1], colors[up * 4 + 2],colors[up * 4 + 3]);
    }

  }
  return c;
}


struct Hexahedron
{
  vec3i lower;
  int level;
};

void importExaJet(const std::shared_ptr<Node> world, const FileName fileName)
{
  int fd               = open(fileName.c_str(), O_RDONLY);
  struct stat stat_buf = {0};
  fstat(fd, &stat_buf);
  const size_t num_hexes = stat_buf.st_size / sizeof(Hexahedron);
  std::cout << "File " << fileName.c_str() << "\n"
            << "size: " << stat_buf.st_size << "\n"
            << "#hexes: " << num_hexes << "\n";
  void *mapping = mmap(NULL, stat_buf.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (mapping == MAP_FAILED) {
    std::cout << "Failed to map file\n";
    perror("mapping file");
    return;
  }

  ospcommon::containers::AlignedVector<vec4f> points;
  ospcommon::containers::AlignedVector<vec4uc> colors;

  points.reserve(num_hexes);
  colors.reserve(num_hexes);

  Hexahedron *hexes = static_cast<Hexahedron *>(mapping);

  ospray::tamr::TAMRData data;
  data.amrOrigin = hexes[0].lower;
  int maxLevel   = 0;

  size_t showVoxelNumber = num_hexes;//* 0.001;

  for (size_t i = 0; i < showVoxelNumber; ++i) {
    Hexahedron &h = hexes[i];

    ospray::tamr::TAMRVoxel voxel;
    voxel.level         = h.level;
    maxLevel            = max(maxLevel, h.level);
    float model2world   = 1.0 / (1 << voxel.level);
    voxel.lower         = (h.lower - data.amrOrigin) * model2world;
    voxel.indexInBuffer = i;

    data.voxelsInLevel[h.level].push_voxel(voxel);
  }
  // Cell width in model space. Scale to 1 in world space
  data.cellScale = (float)(1 << maxLevel);

  for (auto &lv : data.voxelsInLevel) {
    lv.second.cellWidthInModel = (float)(1 << lv.first);
    lv.second.cellWidth        = lv.second.cellWidthInModel / data.cellScale;
    lv.second.halfCellWidth    = 0.5 * lv.second.cellWidth;
    lv.second.rcpCellWidth     = 1.f / lv.second.cellWidth;

    std::cout << "Level " << lv.first << " Num: " << lv.second.voxels.size()
              << " bounds" << lv.second.bounds << "\n";
  }

  munmap(mapping, stat_buf.st_size);
  close(fd);

  TAMRLevelKDT accel(data, 6);
  PRINT(accel.leaf.size());
  int leafIndex(0);
  for(auto & l : accel.leaf)
  {
    //std::cout<<"leaf bounds:" << l.bounds<<std::endl;
    
    vec3f c = findColorForValue(tfn_colors, leafIndex, (int)accel.leaf.size(),true);
    for(auto &v: l.voxels){
        float radii  = 0.5 * accel.level.cellWidth;
        vec3f center = (v.lower + vec3f(0.5)) * accel.level.cellWidth;
        points.push_back(vec4f(center, radii));
        colors.push_back(vec4uc(c.x * 255.0, c.y * 255.0, c.z * 255.0, 255));
    }

    leafIndex++;
  }



  // for (auto &lv : data.voxelsInLevel) {
  //     for (auto &v : lv.second.voxels) {
  //       float radii  = 0.5 * lv.second.cellWidth;
  //       vec3f center = (v.lower + vec3f(0.5)) * lv.second.cellWidth;
  //       points.push_back(vec4f(center, radii));

  //       vec3f c;
  //       switch (lv.first) {
  //       case 3:
  //         c = vec3f(1.0, 0.0, 0.0);
  //         break;
  //       case 4:
  //         c = vec3f(0.0, 1.0, 0.0);
  //         break;
  //       case 5:
  //         c = vec3f(0.0, 0.0, 1.0);
  //         break;
  //       default:
  //         c = vec3f(1.0, 1.0, 0.0);
  //       }
  //       colors.push_back(vec4uc(c.x * 255.0, c.y * 255.0, c.z * 255.0, 255));
  //     }
  // }

  auto exajetGeom = createNode(fileName, "Spheres")->nodeAs<Spheres>();
  exajetGeom->createChild("bytes_per_sphere", "int", int(sizeof(vec4f)));
  exajetGeom->createChild("offset_center", "int", int(0));
  exajetGeom->createChild("radius", "float", 0.5f);
  exajetGeom->createChild("offset_radius", "int", int(sizeof(vec3f)));

  auto materials = exajetGeom->child("materialList").nodeAs<MaterialList>();
  materials->item(0)["Ks"] = vec3f(0.f);

  auto spheres = std::make_shared<DataVectorT<vec4f, OSP_RAW>>();
  spheres->setName("spheres");
  spheres->v = std::move(points);

  auto colorData = std::make_shared<DataVectorT<vec4uc, OSP_UCHAR4>>();
  colorData->setName("color");
  colorData->v = std::move(colors);

  exajetGeom->add(spheres);
  exajetGeom->add(colorData);

  world->add(exajetGeom);
}

extern "C" OSPRAY_DLLEXPORT void ospray_init_module_exajet_import()
{
  std::cout << "Loading NASA exajet importer module\n";
}

OSPSG_REGISTER_IMPORT_FUNCTION(importExaJet, bin);
