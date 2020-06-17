#ifndef MESH_H
#define MESH_H

#include <vector>

#include "Vector2f.h"
#include "Vector3f.h"
#include "object3d.hpp"
#include "triangle.hpp"

class KDTree;

class Mesh : public Object3D {
   public:
    KDTree *kdtree;

    Mesh(const char *filename, Material *m);

    struct TriangleIndex {
        TriangleIndex() {
            x[0] = 0;
            x[1] = 0;
            x[2] = 0;

            texID[0] = -1;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front
        // face
        int x[3]{};
        int texID[3]{};
    };

    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    std::vector<Vector2f> texCoord;
    bool intersect(const Ray &r, Hit &h, double tmin) override;

    // void getAABB(std::vector<AABB> &vec) {
    //     for (int triId = 0; triId < (int)t.size(); ++triId) {
    //         TriangleIndex &triIndex = t[triId];
    //         Triangle *triangle =
    //             new Triangle(v[triIndex[0]], v[triIndex[1]], v[triIndex[2]],
    //                          n[triId], material);
    //         vec.push_back(AABB(triangle));
    //     }
    // }

   private:
    // Normal can be used for light estimation
    void computeNormal();
};

#endif
