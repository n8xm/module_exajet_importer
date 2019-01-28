#include "TAMRData.h"

using namespace ospray;
using namespace ospray::sg;

namespace ospray {
  namespace tamr {

    struct TAMRLevelKDT
    {
      TAMRLevelKDT(const TAMRData &input, int level);
      ~TAMRLevelKDT();

      /*! precomputed values per level, so we can easily compute
      logicla coordinates, find any level's cell width, etc */
      struct Level
      {
        float cellWidthInModel;
        float cellWidth;
        float rcpCellWidth;
        float halfCellWidth;
        int level;
      };

      struct Leaf
      {
        Leaf() {}
        Leaf(const Leaf &l) : voxels(l.voxels), bounds(l.bounds) {}

        std::vector<TAMRVoxel> voxels;
        box3f bounds;
      };

      /*! each node in the tree refers to either a pair ofo child
      nodes (using split plane pos, split plane dim, and offset in
      node[] array), or to a list of 'numItems' bricks (using
      dim=3, and offset pointing to the start of the list of
      numItems brick pointers in 'brickArray[]'. In case of a
      leaf, the 'num' bricks stored in this leaf are (in theory),
      exactly one brick per level, in sorted order */
      struct Node
      {
        inline bool isLeaf() const
        {
          return dim == 3;
        }
        // first dword
        uint32 ofs : 30;  // offset in node[] array (if inner), or brick ID (if
                          // leaf)
        uint32 dim : 2;   // upper two bits: split dimension. '3' means 'leaf
        // second dword
        union
        {
          float pos;
          uint32 numItems;
        };
      };

      //! list of levels
      Level level;
      //! list of inner nodes
      std::vector<Node> node;
      //! list of leaf nodes
      std::vector<Leaf> leaf;
      //! world bounds of domain
      box3f worldBounds;

     private:
      void makeLeaf(index_t nodeID,
                    const box3f &bounds,
                    const std::vector<TAMRVoxel> &voxels);
      void makeInner(index_t nodeID, int dim, float pos, int childID);
      void buildRec(int nodeID,
                    const box3f &bounds,
                    std::vector<TAMRVoxel> voxels);

      float getBestPos(const box3f& bounds, std::vector<TAMRVoxel> voxels, int dim);
    };

  }  // namespace tamr
}  // namespace ospray
