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
};

#endif
