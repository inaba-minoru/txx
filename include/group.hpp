#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <vector>

#include "hit.hpp"
#include "object3d.hpp"
#include "ray.hpp"

// TODO: Implement Group - add data structure to store a list of Object* ------
// OK
class Group : public Object3D {
   private:
    std::vector<Object3D *> objects;

   public:
    Group() = default;
    explicit Group(int num_objects) { objects.resize(num_objects); }

    ~Group() override {
        for (size_t i = 0; i < objects.size(); i++) {
            delete objects[i];
        }
    }

    bool intersect(const Ray &r, Hit &h, double tmin) {
        bool isIntersect = false;
        for (size_t i = 0; i < objects.size(); i++) {
            isIntersect |= objects[i]->intersect(r, h, tmin);
        }
        return isIntersect;
    }

    void addObject(int index, Object3D *obj) { objects[index] = obj; }

    int getGroupSize() { return objects.size(); }

    void init() {
        for (size_t i = 0; i < objects.size(); i++) {
            if (objects[i]->material->emission != Vector3f::ZERO) {
                area += objects[i]->area;
                illum_obj.push_back(objects[i]);
                area_presum.push_back(area);
            }
        }
    }

    std::vector<Object3D *> illum_obj;
    std::vector<double> area_presum;
    void sampleLight(unsigned short *Xi, Vector3f &p, Vector3f &light,
                     double &_pdf) {
        double temp = erand48(Xi) * area;
        int l = -1, r = area_presum.size() - 1;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (area_presum[m] < temp)
                l = m;
            else
                r = m;
        }
        illum_obj[r]->sampleLight(Xi, p, light);
        _pdf = area;
    }
};

#endif
