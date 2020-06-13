#ifndef AABB_H
#define AABB_H

#include <vecmath.h>

#include <ray.hpp>

class AABB {
   private:
    Vector3f pmin, pmax;  // 左下角与右上角

   public:
    AABB(const Vector3f &pmin, const Vector3f &pmax) : pmin(pmin), pmax(pmax) {}
    ~AABB() {}

    bool intersect(const Ray &r) {
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