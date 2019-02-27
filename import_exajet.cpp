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
// #include "importer/Importer.h"
// #include "transferFunction/TransferFunction.h"
// #include "common/Common.h"
// #include "geometry/Spheres.h"
#include "common/sg/volume/Volume.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"

#include "ospcommon/AffineSpace.h"
#include "ospcommon/box.h"
// #include "ospcommon/ospmath.h"
#include "ospcommon/memory/malloc.h"
#include "ospcommon/range.h"
#include "ospcommon/xml/XML.h"
#include "ospray/ospray.h"

#include "sparsepp/spp.h"
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


//NATHAN: Looks like this maps to the Exajet file format. It specifies:
//* The lower left corner of the hex.
//* The AMR level.
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
  //         c = vVolume
  //         breakVolume
  //       case 4:Volume
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


struct HexVert {
  size_t x, y, z;
  HexVert() : x(-1), y(-1), z(-1) {}
  HexVert(size_t x, size_t y, size_t z) : x(x), y(y), z(z) {}
  HexVert(const vec3i &v) : x(v.x), y(v.y), z(v.z) {}
};

bool operator==(const HexVert &a, const HexVert &b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

struct HashHexVert {
  // Hard to really decide how to hash each vert to a unique index,
  // due to the AMR layout..
  size_t operator()(const HexVert &a) const {
    return a.x + 1232128 * (a.y + a.z * 1259072);
  }
};

void importUnstructured(const std::shared_ptr<Node> world, const FileName fileName){
  // Open the hexahedron data file
  int hexFd = open(fileName.c_str(), O_RDONLY);
  struct stat statBuf = {0};
  fstat(hexFd, &statBuf);
  const size_t numHexes = statBuf.st_size / sizeof(Hexahedron);
  std::cout << "File " << fileName.c_str() << "\n"
    << "size: " << statBuf.st_size << "\n"
    << "#hexes: " << numHexes << "\n";
  void *hexMapping = mmap(NULL, statBuf.st_size, PROT_READ, MAP_PRIVATE, hexFd, 0);
  if (hexMapping == MAP_FAILED) {
    std::cout << "Failed to map hexes file\n";
    perror("hex_mapping file");
    return;
  }

  // Open the field data file
  const std::string cellFieldName = "y_vorticity.bin";
  const FileName fieldFile = fileName.path() + cellFieldName;
  std::cout << "Loading field file: " << fieldFile << "\n";
  int fieldFd = open(fieldFile.c_str(), O_RDONLY);
  struct stat fieldStatBuf = {0};
  fstat(fieldFd, &fieldStatBuf);
  std::cout << "File " << fieldFile.c_str() << "\n"
    << "size: " << fieldStatBuf.st_size << "\n";
  void *fieldMapping = mmap(NULL, fieldStatBuf.st_size, PROT_READ,
                            MAP_PRIVATE, fieldFd, 0);
  if (fieldMapping == MAP_FAILED) {
    std::cout << "Failed to map field file\n";
    perror("field_mapping file");
    return;
  }

  ospcommon::containers::AlignedVector<vec3f> verts;
  ospcommon::containers::AlignedVector<vec4i> indices;
  ospcommon::containers::AlignedVector<float> cellVals;
                                                                                   
  const Hexahedron *hexes = static_cast<const Hexahedron*>(hexMapping);
  const float *cellField = static_cast<const float*>(fieldMapping);

  const int desiredLevel = -1;
  const size_t memLimit = 0;//size_t(5)*size_t(1024)*size_t(1024)*size_t(1024);

  // TODO: Using this vertex index re-mapping will help us save
  // a ton of memory and indices by re-using existing vertices, but will
  // really hurt load performance.
#define REMAP_INDICES
#ifdef REMAP_INDICES
  spp::sparse_hash_map<HexVert, int32_t, HashHexVert> vertsMap;
#endif

  const vec3i gridMin = vec3i(1232128, 1259072, 1238336);
  const float voxelScale = 0.0005;
  const vec3f worldMin = vec3f(-1.73575, -9.44, -3.73281);

  for (size_t i = 1; i < numHexes; ++i) {
    const Hexahedron &h = hexes[i];
    if (verts.size() + 8 >= std::numeric_limits<int32_t>::max()) {
      std::cout << "Index size limit reached, terminating mesh load\n";
      break;
    }
    if (desiredLevel == -1 || h.level == desiredLevel) {
      // Pretty inefficient, no vertex re-use, but rendering the octree AMR
      // as an unstructured mesh is a bad route anyways that we won't do beyond
      // quick testing/previewing
      const vec3i hexSize = vec3i(1 << h.level);

      // Verts ordering for a hex cell:
      // four bottom verts counter-clockwise
      // four top verts counter-clockwise
      for (int k = 0; k < 2; ++k) {
        vec4i idx;
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < 2; ++i) {
            // We want to go x_low -> x_hi if y_low, and x_hi -> x_low if y_hi
            const int x = (i + j) % 2;
            const vec3i p = h.lower + hexSize * vec3i(x, j, k);
            const vec3f worldPos = vec3f(p - gridMin) * voxelScale + worldMin;

#ifdef REMAP_INDICES
            HexVert hexVert(p);
            auto fnd = vertsMap.find(hexVert);
            if (fnd == vertsMap.end()) {
              vertsMap[hexVert] = verts.size();
              idx[j * 2 + i] = verts.size();
              verts.push_back(worldPos);
            } else {
              idx[j * 2 + i] = fnd->second;
            }
#else
            idx[j * 2 + i] = verts.size();
            verts.push_back(worldPos);
#endif
          }
        }
        indices.push_back(idx);
      }
      cellVals.push_back(cellField[i]);
    }
    const size_t memSize = verts.size() * sizeof(vec3f)
                           + indices.size() * sizeof(vec4i)
                           + cellVals.size() * sizeof(float);
    if (memLimit != 0 && memSize >= memLimit) {
      break;
    }
  }

  std::cout << "Imported " << cellVals.size() << " hexahedrons\n";

  munmap(hexMapping, statBuf.st_size);
  munmap(fieldMapping, fieldStatBuf.st_size);
  close(hexFd);
  close(fieldFd);

  //NATHAN: Here is where we create the unstructured volume. This code uses
  //OSPRay's scene graph functionality. We should probably avoid using OSPRay's
  //scene graph for now because we are starting out by rendering just one volume.
  //Furthermore, Will ays that the scene graph is poorly documented (I think
  //that's what he said?)  
 
  //NATHAN: "Jet" appears to refer to the jet airplane in the data set, not the
  //"jet" colormap.
  auto jet = createNode(fileName, "UnstructuredVolume")->nodeAs<Volume>();
  jet->createChild("cellFieldName", "string", cellFieldName);

  auto vertsData = std::make_shared<DataVectorT<vec3f, OSP_FLOAT3>>();
  vertsData->setName("vertices");
  vertsData->v = std::move(verts);

  auto indicesData = std::make_shared<DataVectorT<vec4i, OSP_INT4>>();
  indicesData->setName("indices");
  indicesData->v = std::move(indices);

  auto cellFieldData = std::make_shared<DataVector1f>();
  cellFieldData->setName("0");
  cellFieldData->v = std::move(cellVals);

  //NATHAN: Since we are using cellField, I believe that we are using cell-centered data here.
  auto fieldList = std::make_shared<NodeList<DataVector1f>>();
  fieldList->setName("cellFields");
  fieldList->push_back(cellFieldData);

  std::vector<sg::Any> cellFieldNames = {cellFieldName};
  jet->createChild("cellFieldName", "string", cellFieldName).setWhiteList(cellFieldNames);

  jet->add(vertsData);
  jet->add(indicesData);
  jet->add(fieldList);

  world->add(jet);
}

extern "C" OSPRAY_DLLEXPORT void ospray_init_module_exajet_import() {
  std::cout << "Loading NASA exajet importer module\n";
}

OSPSG_REGISTER_IMPORT_FUNCTION(importExaJet, bin);
OSPSG_REGISTER_IMPORT_FUNCTION(importUnstructured, jetunstr);

