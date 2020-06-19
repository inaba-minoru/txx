#ifndef PT_H
#define PT_H

#include <assert.h>

#include <ggx.hpp>
#include <group.hpp>
#include <hit.hpp>
// #include <kdtree.hpp>
#include <vecmath.h>

#include <light.hpp>
#include <ray.hpp>
#include <scene_parser.hpp>

Vector3f min(const Vector3f &x, const Vector3f &y) {
    return Vector3f(std::min(x.x(), y.x()), std::min(x.y(), y.y()),
                    std::min(x.z(), y.z()));
}

Vector3f calcDiffuseOutDir(const Vector3f &in_dir, const Vector3f &norm) {
    return (cbrt(drand48()) * (2 * Vector3f(drand48(), drand48(), drand48()) -
                               Vector3f(1, 1, 1))
                                  .normalized() +
            norm)
        .normalized();
}
Vector3f calcReflectOutDir(const Vector3f &in_dir, const Vector3f &norm) {
    return 2 * Vector3f::dot(in_dir, norm) * norm - in_dir;
}
bool calcRefractOutDir(const Vector3f &in_dir, const Vector3f &norm,
                       const double &eta, Vector3f &out_dir) {
    double cos1 = Vector3f::dot(-in_dir, norm);
    double cos2 = sqrt(1 - eta * eta * (1 - cos1 * cos1));

    if (cos2 > 0) {
        out_dir = eta * in_dir + norm * (eta * cos1 - cos2);
        return 1;
    }

    return 0;
}

// Vector3f radiance(unsigned short *Xi, const SceneParser &scene_parser,
//                   const KDTree &kdtree, const Ray &ray, const double &tmin,
//                   const int &depth, double n) {
Vector3f radiance(unsigned short *Xi, const SceneParser &scene_parser,
                  const Ray &ray, const double &tmin, const int &depth,
                  double n, bool accept_light) {
    if (depth > 8) return Vector3f::ZERO;

    // 没有来源
    Hit hit;
    if (!scene_parser.getGroup()->intersect(ray, hit, tmin))
        return Vector3f::ZERO;

    Vector3f hit_point = ray.pointAtParameter(hit.getT());
    Material *material = hit.getMaterial();

    if (material->emission != Vector3f::ZERO) {
        if (accept_light)
            return material->color;
        else
            return Vector3f::ZERO;
    }

    // 采样光源
    Vector3f light_p, light_color;
    double light__pdf;
    scene_parser.getGroup()->sampleLight(Xi, light_p, light_color, light__pdf);
    // //  light_p.print();
    // puts("here");
    Vector3f temp = light_p - hit_point;
    Vector3f light_dir = temp.normalized();
    Ray temp_ray(hit_point, light_dir);
    Hit temp_hit;
    scene_parser.getGroup()->intersect(temp_ray, temp_hit, tmin);
    Vector3f dir_light;
    if ((temp_ray.pointAtParameter(temp_hit.getT()) - light_p).length() <
        1e-4) {
        dir_light = light_color *
                    material->f(light_dir, hit.getNormal(), -ray.getDirection(),
                                n > 1 ? 1 : material->eta, n) *
                    Vector3f::dot(light_dir, hit.getNormal()) *
                    Vector3f::dot(-light_dir, temp_hit.getNormal()) /
                    temp.squaredLength() * light__pdf;

        if (depth) dir_light = min(dir_light, light_color);
    }

    Vector3f obj_color;
    if (hit.hasTex && material->texture.valid())
        obj_color = material->texture(hit.texCoord.x(), hit.texCoord.y());
    else
        obj_color = material->color;

    // double P_RR = std::min(
    //     std::max(obj_color.x(), std::max(obj_color.y(), obj_color.z())),
    //     0.99);
    // if (depth > 4)
    //     if (erand48(Xi) > P_RR)
    //         return Vector3f::ZERO;
    //     else
    //         obj_color *= P_RR;

    // 漫反射、反射、折射提供来源（间接）

    Vector3f o;
    double weight, new_eta = n;
    // sample(Xi, -ray.getDirection(), hit.getNormal(), o, weight, new_eta,
    //        n > 1 ? 1 : material->eta, n, material->alpha_g,
    //        material->alpha_g2);
    material->sample(Xi, -ray.getDirection(), hit.getNormal(), n,
                     n > 1 ? 1 : material->eta, o, weight, new_eta,
                     accept_light);

    Vector3f color = radiance(Xi, scene_parser, Ray(hit_point, o), tmin,
                              depth + 1, new_eta, accept_light) *
                     weight * obj_color;

    return color + material->emission + dir_light * obj_color;
}

#endif  // PT_H
