#ifndef PT_H
#define PT_H

#include <hit.hpp>
#include <kdtree.hpp>
#include <light.hpp>
#include <ray.hpp>
#include <scene_parser.hpp>

Vector3f calcDiffuseOutDir(const Vector3f &in_dir, const Vector3f &norm) {
    return (cbrt(drand48()) * (2 * Vector3f(drand48(), drand48(), drand48()) -
                               Vector3f(1, 1, 1))
                                  .normalized() -
            norm)
        .normalized();
}
Vector3f calcReflectOutDir(const Vector3f &in_dir, const Vector3f &norm) {
    return in_dir + 2 * Vector3f::dot(-in_dir, norm) * norm;
}
bool calcRefractOutDir(const Vector3f &in_dir, const Vector3f &norm,
                       const float &eta, Vector3f &out_dir) {
    // Vector3f vertical = (Vector3f::dot(-in_dir, norm) * norm + in_dir) / eta;
    // out_dir = cos(asin(vertical.length())) * -norm + vertical;
    // return 1;
    float cos1 = Vector3f::dot(-in_dir, norm);
    float cos2 = sqrt(1 - eta * eta * (1 - cos1 * cos1));

    if (cos2 > 0) {
        out_dir = eta * in_dir + norm * (eta * cos1 - cos2);
        return 1;
    }

    return 0;
}

Vector3f radiance(const SceneParser &scene_parser, const KDTree &kdtree,
                  const Ray &ray, const float &tmin, const int &depth,
                  const float &n) {
    if (depth > 8) return Vector3f::ZERO;

    // 没有来源
    Hit hit;
    if (!kdtree.intersect(ray, hit, tmin)) return Vector3f::ZERO;
    //

    Vector3f hit_point = ray.pointAtParameter(hit.getT());
    Material *material = hit.getMaterial();

    // 光源提供来源（直接）
    Vector3f light_color = Vector3f::ZERO;
    for (int i = 0; i < scene_parser.getNumLights(); ++i) {
        Vector3f dir, col;
        Hit temphit;
        scene_parser.getLight(i)->getIllumination(hit_point, dir, col);
        if (!kdtree.intersect(Ray(hit_point, dir), temphit, tmin))  // 无遮挡
            light_color += hit.getMaterial()->Shade(ray, hit, dir, col);
    }
    //

    // 漫反射、反射、折射提供来源（间接）

    Vector3f diffuse_color = Vector3f::ZERO;
    Vector3f reflect_color = Vector3f::ZERO;
    Vector3f refract_color = Vector3f::ZERO;
    if (drand48() < 0.5) {
        Ray diffuse_ray = Ray(
            hit_point, calcDiffuseOutDir(ray.getDirection(), hit.getNormal()));
        diffuse_color = 0.5 * radiance(scene_parser, kdtree, diffuse_ray, tmin,
                                       depth + 1, n);
    } else if (drand48() < 0.0) {
        Ray reflect_ray = Ray(
            hit_point, calcReflectOutDir(ray.getDirection(), hit.getNormal()));
        reflect_color = 0.5 * radiance(scene_parser, kdtree, reflect_ray, tmin,
                                       depth + 1, n);
    } else {
        Vector3f out_dir;
        if (calcRefractOutDir(ray.getDirection(), hit.getNormal(),
                              n < 1 ? 1 / 0.8 : 0.8, out_dir)) {
            Ray refract_ray = Ray(hit_point, out_dir);
            refract_color = radiance(scene_parser, kdtree, refract_ray, tmin,
                                     depth + 1, n < 1 ? 1 : 0.8);
        } else {
            Ray reflect_ray =
                Ray(hit_point,
                    calcReflectOutDir(ray.getDirection(), hit.getNormal()));
            reflect_color = 0.5 * radiance(scene_parser, kdtree, reflect_ray,
                                           tmin, depth + 1, n);
        }
    }

    //
    return light_color + diffuse_color + reflect_color + refract_color;
}

#endif  // PT_H