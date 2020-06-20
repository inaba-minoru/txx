#include "mesh.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <kdtree.hpp>
#include <sstream>
#include <utility>
#include <vector>

bool Mesh::intersect(const Ray &r, Hit &h, double tmin) {
    return kdtree->intersect(r, h, tmin);

    // Optional: Change this brute force method into a faster one.
    // bool result = false;
    // for (int triId = 0; triId < (int)t.size(); ++triId) {
    //     // TriangleIndex &triIndex = t[triId];
    //     // Triangle triangle(v[triIndex[0]], v[triIndex[1]], v[triIndex[2]],
    //     //   material);
    //     // triangle.normal = n[triId];
    //     // result |= triangle.intersect(r, h, tmin);

    //     result |= triangles[triId]->intersect(r, h, tmin);
    // }
    // return result;
}

Mesh::Mesh(const char *filename, Material *material) : Object3D(material) {
    // Optional: Use tiny obj loader to replace this simple one.
    std::ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        std::cout << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::string vTok("v");
    std::string fTok("f");
    std::string texTok("vt");
    char bslash = '/', space = ' ';
    std::string tok;
    // int texID;
    bool flag = 0;

    while (true) {
        std::getline(f, line);
        if (f.eof()) {
            break;
        }
        if (line.size() < 3) {
            continue;
        }
        if (line.at(0) == '#') {
            continue;
        }
        std::stringstream ss(line);
        ss >> tok;
        if (tok == vTok) {
            Vector3f vec;
            ss >> vec[0] >> vec[1] >> vec[2];
            v.push_back(vec);
            n.push_back(0);
        } else if (tok == fTok) {
            if (line.find(bslash) != std::string::npos) {
                std::replace(line.begin(), line.end(), bslash, space);
                std::stringstream facess(line);
                TriangleIndex trig;
                facess >> tok;
                for (int ii = 0; ii < 3; ii++) {
                    facess >> trig[ii] >> trig.texID[ii];
                    trig[ii]--;
                    trig.texID[ii]--;
                }
                t.push_back(trig);
                double t1, t2;
                if (facess >> trig[1] >> trig.texID[1]) {
                    trig[1]--;
                    trig.texID[1]--;
                    std::swap(trig[0], trig[2]);
                    std::swap(trig.texID[0], trig.texID[2]);
                    t.push_back(trig);
                }

            } else {
                TriangleIndex trig;
                for (int ii = 0; ii < 3; ii++) {
                    ss >> trig[ii];
                    trig[ii]--;
                }
                t.push_back(trig);
                if (ss >> trig[1]) {
                    trig[1]--;
                    std::swap(trig[0], trig[2]);
                    t.push_back(trig);
                }
            }
        } else if (tok == texTok) {
            Vector2f texcoord;
            ss >> texcoord[0];
            ss >> texcoord[1];
            texCoord.push_back(texcoord);
        } else if (tok == "flag")
            flag = 1;
    }
    computeNormal();

    f.close();

    // build kdtree
    vector<AABB> vec;
    for (unsigned int triId = 0; triId < t.size(); ++triId) {
        TriangleIndex &triIndex = t[triId];
        Triangle *triangle;
        if (flag)
            triangle = new Triangle(v[triIndex[0]], v[triIndex[1]],
                                    v[triIndex[2]], material);
        else
            triangle = new Triangle(v[triIndex[0]], v[triIndex[1]],
                                    v[triIndex[2]], n[triIndex[0]],
                                    n[triIndex[1]], n[triIndex[2]], material);
        if (triIndex.texID[0] != -1)
            triangle->setTexCoords(texCoord[triIndex.texID[0]],
                                   texCoord[triIndex.texID[1]],
                                   texCoord[triIndex.texID[2]]);
        vec.push_back(AABB(triangle));

        triangles.push_back(triangle);
        area += triangle->area;
        area_presum.push_back(area);
    }
    kdtree = new KDTree(vec.begin(), vec.size());
    // end build kdtree
}

void Mesh::computeNormal() {
    n.resize(t.size());
    for (int triId = 0; triId < (int)t.size(); ++triId) {
        TriangleIndex &triIndex = t[triId];
        Vector3f a = v[triIndex[1]] - v[triIndex[0]];
        Vector3f b = v[triIndex[2]] - v[triIndex[0]];
        b = Vector3f::cross(a, b);
        for (size_t i = 0; i < 3; i++) n[triIndex[i]] += b;

        // n[triId] = b / b.length();
    }

    for (size_t i = 0; i < v.size(); i++) n[i].normalize();
}
