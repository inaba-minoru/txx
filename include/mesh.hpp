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

    std::vector<Triangle *> triangles;
    std::vector<double> area_presum;
    void sampleLight(unsigned short *Xi, Vector3f &p, Vector3f &light) {
        double temp = erand48(Xi) * area;
        int l = -1, r = area_presum.size() - 1;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (area_presum[m] < temp)
                l = m;
            else
                r = m;
        }
        triangles[r]->sampleLight(Xi, p, light);
    }

    void transform(const Matrix4f &mat) {
        area = 0;
        for (int i = 0; i < triangles.size(); i++) {
            triangles[i]->transform(mat);
            area += triangles[i]->area;
            area_presum[i] = area;
        }
    }

   private:
    // Normal can be used for light estimation
    void computeNormal();
};

#endif
