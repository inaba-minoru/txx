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
    return -in_dir + 2 * (Vector3f::dot(-in_dir, norm) * norm + in_dir);
}
Vector3f calcRefractOutDir(const Vector3f &in_dir, const Vector3f &norm) {
    return in_dir;
}

Vector3f radiance(const SceneParser &scene_parser, const KDTree &kdtree,
                  const Ray &ray, const int &tmin) {
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

    // Ray diffuse_ray =
    //     Ray(hit_point, calcDiffuseOutDir(ray.getDirection(),
    //     hit.getNormal()));
    // Vector3f diffuse_color = radiance(scene_parser, kdtree, diffuse_ray,
    // tmin);

    Ray reflect_ray =
        Ray(hit_point, calcReflectOutDir(ray.getDirection(), hit.getNormal()));
    Vector3f reflect_color = radiance(scene_parser, kdtree, reflect_ray, tmin);

    // Ray refract_ray =
    //     Ray(hit_point, calcRefractOutDir(ray.getDirection(),
    //     hit.getNormal()));
    // Vector3f refract_color = radiance(scene_parser, kdtree, refract_ray,
    // tmin);

    //

    // return diffuse_color;
    return light_color + reflect_color;
}

#endif  // PT_H