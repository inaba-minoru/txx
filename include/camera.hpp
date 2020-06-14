#ifndef CAMERA_H
#define CAMERA_H

#include <float.h>
#include <vecmath.h>

#include <cmath>

#include "ray.hpp"

class Camera {
   public:
    Camera(const Vector3f &center, const Vector3f &direction,
           const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up);
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

   protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// TODO: Implement Perspective camera ------ OK
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {
   private:
    float f;

   public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle)
        : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        f = imgH / tan(angle / 2) / 2;
    }

    Ray generateRay(const Vector2f &point) override {
        Vector3f dc = Vector3f(point.x() - width / 2, point.y() - height / 2, f)
                          .normalized();
        return Ray(center, Matrix3f(horizontal, up, direction) * dc);
    }
};

#endif  // CAMERA_H
