#include <sys/types.h>                                                          
#include <sys/stat.h>                                                           
#include <sys/mman.h>                                                           
#include <fcntl.h>                                                              
#include <unistd.h>   

#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "importer/Importer.h"
#include "transferFunction/TransferFunction.h"
#include "common/Common.h"
#include "geometry/Spheres.h"
#include "volume/Volume.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"

#include "ospcommon/AffineSpace.h"
#include "ospcommon/box.h"
#include "ospcommon/ospmath.h"
#include "ospcommon/memory/malloc.h"
#include "ospcommon/range.h"
#include "ospcommon/xml/XML.h"
#include "ospray/ospray.h"


#include "TAMRData.h"


using namespace ospcommon;
using namespace ospray;
using namespace ospray::sg;

using namespace ospray::tamr;


struct Hexahedron {  
  vec3i lower;                                                              
  int level;                                                              
};


void importExaJet(const std::shared_ptr<Node> world, const FileName fileName){
  int fd = open(fileName.c_str(), O_RDONLY);                                                
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
                                                                                   
  Hexahedron *hexes = static_cast<Hexahedron*>(mapping);                           
  std::unordered_map<int, box3i> level_bounds;     

  TAMRData data; 
  data.amrOrigin = hexes[0].lower;
  data.voxels.reserve(num_hexes);
  for (size_t i = 0; i < num_hexes; ++i) {                                      
    Hexahedron &h = hexes[i];                                             

    TAMRVoxel voxel;
    voxel.level         = h.level;
    float model2world   = 1.0 /(1 << voxel.level);
    voxel.lower         = (h.lower - data.amrOrigin) * model2world * voxel.getVoxelWidth() ;
    voxel.indexInBuffer = i;
    data.voxels.push_back(voxel);

    level_bounds[h.level].extend(voxel.lower);   

    // if(h.level ==6)
    // {
    points.push_back(vec4f(voxel.getCenter(),0.5 * voxel.getVoxelWidth()));

    vec3f c;
    switch (h.level) {
    case 3:
      c = vec3f(1.0, 0.0, 0.0);
      break;
    case 4:
      c = vec3f(0.0, 1.0, 0.0);
      break;
    case 5:
      c = vec3f(0.0, 0.0, 1.0);
      break;
    default:
      c = vec3f(1.0, 1.0, 0.0);
    }
    colors.push_back(vec4uc(c.x * 255.0, c.y * 255.0, c.z * 255.0, 255));
    // }
  }

  for (auto &l : level_bounds) {     
    float rcpWidth = 1.0 / (1 << l.first);
    box3f box = box3f(l.second.lower * rcpWidth, l.second.upper * rcpWidth);                                       
    std::cout << "Level " << l.first << " bounds " << l.second <<"\n";         
  }                                                                                                                                                          
  munmap(mapping, stat_buf.st_size);                                            
  close(fd); 

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
  const std::string cellFieldName = "x_vorticity.bin";
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

  for (size_t i = 0; i < numHexes; ++i) {
    const Hexahedron &h = hexes[i];
    if (desiredLevel != -1 && h.level == desiredLevel) {
      // Pretty inefficient, no vertex re-use, but rendering the octree AMR
      // as an unstructured mesh is a bad route anyways that we won't do beyond
      // quick testing/previewing
      const vec3i hexSize = vec3i(1 << h.level);
      const vec3i upper = h.lower + hexSize;
      // Verts ordering for a hex cell:
      // four bottom verts counter-clockwise
      vec4i idx;
      idx.x = verts.size();
      verts.push_back(vec3f(h.lower));

      idx.y = verts.size();
      verts.push_back(vec3f(upper.x, h.lower.y, h.lower.z));

      idx.z = verts.size();
      verts.push_back(vec3f(upper.x, upper.y, h.lower.z));

      idx.w = verts.size();
      verts.push_back(vec3f(h.lower.x, upper.y, h.lower.z));
      indices.push_back(idx);


      // four top verts counter-clockwise
      idx.x = verts.size();
      verts.push_back(vec3f(h.lower.x, h.lower.y, upper.z));

      idx.y = verts.size();
      verts.push_back(vec3f(upper.x, h.lower.y, upper.z));

      idx.z = verts.size();
      verts.push_back(vec3f(upper.x, upper.y, upper.z));

      idx.w = verts.size();
      verts.push_back(vec3f(h.lower.x, upper.y, upper.z));
      indices.push_back(idx);

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

