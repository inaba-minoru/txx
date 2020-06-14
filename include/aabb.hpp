#ifndef AABB_H
#define AABB_H

#include <vecmath.h>

#include <object3d.hpp>
#include <ray.hpp>
#include <sphere.hpp>
#include <triangle.hpp>

class AABB {
   private:
   public:
    Vector3f pmin, pmax, pmid;  // 左下角与右上角
    Object3D *obj = nullptr;

    AABB() {}
    AABB(const Vector3f &pmin, const Vector3f &pmax)
        : pmin(pmin), pmax(pmax), pmid((pmin + pmax) / 2) {}
    AABB(Triangle *obj) : obj(obj) {
        pmin = min(min(obj->vertices[0], obj->vertices[1]), obj->vertices[2]);
        pmax = max(max(obj->vertices[0], obj->vertices[1]), obj->vertices[2]);
        pmid = (pmin + pmax) / 2;
    }
    AABB(Sphere *obj) : obj(obj) {
        pmin = Vector3f(obj->centre.x() - obj->radius,
                        obj->centre.y() - obj->radius,
                        obj->centre.z() - obj->radius);
        pmax = Vector3f(obj->centre.x() + obj->radius,
                        obj->centre.y() + obj->radius,
                        obj->centre.z() + obj->radius);
        pmid = obj->centre;
    }
    ~AABB() {}

    Vector3f min(const Vector3f &x, const Vector3f &y) {
        return Vector3f(std::min(x.x(), y.x()), std::min(x.y(), y.y()),
                        std::min(x.z(), y.z()));
    }
    Vector3f max(const Vector3f &x, const Vector3f &y) {
        return Vector3f(std::max(x.x(), y.x()), std::max(x.y(), y.y()),
                        std::max(x.z(), y.z()));
    }

    bool intersect(const Ray &r) const {
        Vector3f orig = r.getOrigin(), dir = r.getDirection();
        float tin = __FLT_MIN__, tout = __FLT_MAX__, t1, t2;

        for (int i = 0; i < 3; i++) {
            if (dir[i] == 0)                                   // 平行
                if (orig[i] >= pmin[i] && orig[i] <= pmax[i])  // 内部
                    continue;                                  //
                else                                           // 外部
                    return 0;                                  //

            t1 = (pmin[i] - orig[i]) / dir[i];
            t2 = (pmax[i] - orig[i]) / dir[i];
            if (t1 > t2) std::swap(t1, t2);  // t1 in , t2 out
            tin = std::max(tin, t1);
            tout = std::min(tout, t2);
        }

        return tin <= tout;
    }
};

#endif  // AABB_H