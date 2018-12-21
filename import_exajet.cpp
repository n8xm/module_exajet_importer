#include <iostream>                                                             
#include <algorithm>                                                            
#include <limits>                                                               
#include <array>                                                                
#include <unordered_map>                                                        
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


#include "common/sg/importer/Importer.h"
#include "common/sg/transferFunction/TransferFunction.h"
#include "common/sg/common/Common.h"
#include "common/sg/geometry/Spheres.h"
#include "ospcommon/containers/AlignedVector.h"
#include "ospcommon/vec.h"

using namespace ospray;
using namespace ospray::sg;

struct HexSphere{
  vec3f center;
  float radius;
};

struct Hexahedron {  
  vec3i lower;                                                              
  int level;                                                              
  
  HexSphere getHexSphere(){
    HexSphere hs; 
    int size = 1 << level;
    hs.radius = size * 0.5;
    hs.center = lower + vec3f(hs.radius);
    return hs;
  }

};



struct Box {                                                                       
  std::array<int, 3> lower, upper;                                                 
                                                                                   
  Box() {                                                                          
    lower.fill(std::numeric_limits<int>::max());                                   
    upper.fill(std::numeric_limits<int>::min());                                   
  }                                                                                
  void extend(int x, int y, int z) {                                               
    lower[0] = std::min(lower[0], x);                                              
    lower[1] = std::min(lower[1], y);                                              
    lower[2] = std::min(lower[2], z);                                              
                                                                                   
    upper[0] = std::max(upper[0], x);                                              
    upper[1] = std::max(upper[1], y);                                              
    upper[2] = std::max(upper[2], z);                                              
  }                                                                                
};  

std::ostream& operator<<(std::ostream &os, const Box &b) {                         
  os << "[" << b.lower[0] << ", " << b.lower[1] << ", " << b.lower[2]              
    << "], [" << b.upper[0] << ", " << b.upper[1] << ", " << b.upper[2] << "]"; 
  return os;                                                                       
} 



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
  std::unordered_map<int, Box> level_bounds;                                    
  for (size_t i = 0; i < num_hexes; ++i) {                                      
    Hexahedron &h = hexes[i];                                             
    level_bounds[h.level].extend(h.lower.x, h.lower.y, h.lower.z);      

    HexSphere hs = h.getHexSphere();

    points.push_back(vec4f(hs.center,hs.radius));
    //points.push_back(vec3f(h.lower.x,h.lower.y,h.lower.z));          

    vec3f c;
    switch(h.level){
      case 3:
        c = vec3f(1.0,0.0,0.0);
        break;
      case 4:
        c = vec3f(0.0,1.0,0.0);
        break;
      case 5:
        c = vec3f(0.0,0.0,1.0);
        break;
      default:
        c = vec3f(1.0,1.0,0.0);
    }
    colors.push_back(vec4uc(c.x * 255.0, c.y * 255.0, c.z * 255.0, 255));                
  }                                                                             
                                                                                
  for (auto &l : level_bounds) {                                                
    std::cout << "Level " << l.first << " bounds " << l.second << "\n";         
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

extern "C" OSPRAY_DLLEXPORT void ospray_init_module_exajet_import() {
  std::cout << "Loading NASA exajet importer module\n";
}

OSPSG_REGISTER_IMPORT_FUNCTION(importExaJet, bin);

