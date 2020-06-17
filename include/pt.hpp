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
                  double n) {
    // if (erand48(Xi) < 0.01) return 0;
    if (depth > 8) return Vector3f::ZERO;

    // 没有来源
    Hit hit;
    // if (!kdtree.intersect(ray, hit, tmin)) return Vector3f::ZERO;
    if (!scene_parser.getGroup()->intersect(ray, hit, tmin))
        return Vector3f::ZERO;
    //

    Vector3f hit_point = ray.pointAtParameter(hit.getT());
    Material *material = hit.getMaterial();
    // hit.getNormal().normalized().print();
    assert(Vector3f::dot(-ray.getDirection(), hit.getNormal()) >= 0);

    // if (hit_point.z() < 100 && hit_point.z() >70) hit_point.print();
    // Ray rr(Vector3f(50,52,295.6), Vector3f(0.506013,0.385619,0.771524));
    // Hit hh;
    // scene_parser.getGroup()->intersect(rr, hh, tmin);
    // rr.pointAtParameter(hh.getT()).print();
    // rr.getDirection().print();
    // rr.pointAtParameter(1).print();

    // if (material->emission.length() != 0)

    // 光源提供来源（直接）
    Vector3f light_color = Vector3f::ZERO;
    for (int i = 0; i < scene_parser.getNumLights(); ++i) {
        Vector3f dir, col;
        Hit temphit;
        scene_parser.getLight(i)->getIllumination(hit_point, dir, col);
        if (!scene_parser.getGroup()->intersect(Ray(hit_point, dir), temphit,
                                                tmin))  //无遮挡
            light_color += hit.getMaterial()->Shade(ray, hit, dir, col);
        // light_color += hit.getMaterial()->CookTorrance(ray, hit, dir, col);
        // light_color += fs(dir, -ray.getDirection(), hit.getNormal(),
        //                   n > 1 ? 1 : 30, n, 0.0001) *
        //                col * material->diffuseColor;
        // *M_PI;
    }
    //

    // 漫反射、反射、折射提供来源（间接）

    // Vector3f o = calcReflectOutDir(-ray.getDirection(), hit.getNormal());
    Vector3f o;
    //  = calcDiffuseOutDir(-ray.getDirection(), hit.getNormal());
    double weight, new_eta;
    sample(Xi, -ray.getDirection(), hit.getNormal(), o, weight, new_eta,
           n > 1 ? 1 : material->eta, n, material->alpha_g, material->alpha_g2);
    Vector3f myColor;
    if (hit.hasTex && material->texture.valid())
        myColor = material->texture(hit.texCoord.x(), hit.texCoord.y());
    else
        myColor = material->diffuseColor;

    Vector3f color = radiance(Xi, scene_parser, Ray(hit_point, o), tmin,
                              depth + 1, new_eta) *
                     //  material->diffuseColor *
                     //   fs(o, -ray.getDirection(), hit.getNormal(),
                     //  n > 1 ? 1 : material->eta, n, material->alpha_g2) *
                     weight * myColor;

    // if (hit.getNormal().y()<-0.9)
    // {
    //     #pragma omp critical
    //     {
    //     printf("y %f\n", weight);
    //     color.print();
    //     }
    // }
    // if (hit.getNormal().z()>0.9)
    // {
    //     #pragma omp critical
    //     {
    //     printf("z %f\n", weight);
    //     color.print();
    //     }
    // }

    //  *
    //  Vector3f::dot(-ray.getDirection(), hit.getNormal());

    // double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
    // Vector3f w = hit.getNormal(),
    //          u = Vector3f::cross(
    //                  fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0,
    //                  0), w) .normalized(),
    //          v = Vector3f::cross(w, u).normalized();
    // Vector3f d =
    //     (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 -
    //     r2)).normalized();
    // return radiance(Xi, scene_parser, kdtree, Ray(hit_point, d), tmin,
    //                 depth + 1, n) *
    //            material->diffuseColor +
    //        material->emission;
    // ;

    //  material->emission.print();
    // color.print();

    // if (!depth) {
    // ray.getDirection().print();
    // hit_point.print();
    // material->diffuseColor.print();
    // color.print();
    // }

    //
    return color + material->emission + light_color;
}

#endif  // PT_H
