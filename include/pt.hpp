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

Vector3f radiance(unsigned short *Xi, const SceneParser &scene_parser, Ray ray,
                  const double &tmin, int depth, double n, bool accept_light) {
    Vector3f ret, k = 1;

    while (1) {
        if (depth > 8) break;

        // 没有来源
        Hit hit;
        if (!scene_parser.getGroup()->intersect(ray, hit, tmin)) break;
        //
        

        Vector3f hit_point = ray.pointAtParameter(hit.getT());
        Material *material = hit.getMaterial();

        if (material->emission != Vector3f::ZERO) {
            if (accept_light) ret += k * material->emission;
            break;
        }

        // 采样光源
        Vector3f dir_light;
        if (scene_parser.getGroup()->illum_obj.size()) {
            Vector3f light_p, light_color;
            double light__pdf;
            scene_parser.getGroup()->sampleLight(Xi, light_p, light_color,
                                                 light__pdf);

            Vector3f temp = light_p - hit_point;
            Vector3f light_dir = temp.normalized();
            Ray temp_ray(hit_point, light_dir);
            Hit temp_hit;
            scene_parser.getGroup()->intersect(temp_ray, temp_hit, tmin);
            if ((temp_ray.pointAtParameter(temp_hit.getT()) - light_p)
                    .length() < 1e-8) {
                dir_light =
                    light_color *
                    material->f(light_dir, hit.getNormal(), -ray.getDirection(),
                                n > 1 ? 1 : material->eta, n) *
                    Vector3f::dot(light_dir, hit.getNormal()) *
                    Vector3f::dot(-light_dir, temp_hit.getNormal()) /
                    temp.squaredLength() * light__pdf;

                if (depth) dir_light = min(dir_light, light_color);
            }
        }

        Vector3f obj_color;
        if (hit.hasTex && material->texture.valid())
            obj_color = material->texture(hit.texCoord.x(), hit.texCoord.y());
        else
            obj_color = material->color;

        double P_RR =
            std::max(obj_color.x(), std::max(obj_color.y(), obj_color.z()));
        if (depth > 4)
            if (erand48(Xi) > P_RR)
                break;
            else
                obj_color *= P_RR;

        // 漫反射、反射、折射提供来源（间接）

        Vector3f o;
        double weight, new_eta = n;
        material->sample(Xi, -ray.getDirection(), hit.getNormal(), n,
                         n > 1 ? 1 : material->eta, o, weight, new_eta,
                         accept_light);

        ray = Ray(hit_point, o);
        depth++;
        n = new_eta;

        ret += k * dir_light * obj_color;
        //  + obj_color*k;
        k = k * weight * obj_color;

        if (k.x() < 1e-4 && k.y() < 1e-4 && k.z() < 1e-4) break;
    }

    return ret;
}

#endif  // PT_H
