#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <kdtree.hpp>
#include <pt.hpp>
#include <string>

#include "camera.hpp"
#include "group.hpp"
#include "image.hpp"
#include "light.hpp"
#include "scene_parser.hpp"

// using namespace std;

// #include <ggx.hpp>

int main(int argc, char* argv[]) {
    // Vector3f i = Vector3f(0, -.5, 1).normalized();
    // Vector3f n = Vector3f(0, 0, 1).normalized();
    // Vector3f m = Vector3f(0, 0, 1).normalized();
    // double im = Vector3f::dot(i, m);
    // double in = Vector3f::dot(i, n);
    // double eta_i = 30;
    // double eta_t = 1;
    // double eta = eta_i / eta_t;
    // // Vector3f o =
    // //     (eta * im - (in >= 0 ? 1 : -1) * sqrt(1 + eta * eta * (im * im
    // -
    // // 1)))
    // //     *
    // //         m -
    // //     eta * i;
    // Vector3f o = 2 * im * m - i;
    // i.print();
    // o.normalized().print();
    // double on = Vector3f::dot(o, n);
    // Vector3f ht = -(eta_i * i + eta_t * o).normalized();
    // ht.print();
    // // (2 * im * m - i).print();
    // // printf("%f\n", im * G(i, o, m, n, 0.1) / (in * Vector3f::dot(m,
    // n)));
    // printf("%f\n", F(i, m, eta_t, eta_i));
    // // printf("%f\n", D(ht, n, 0.1));
    // printf("%f\n", ft(i, o, n, eta_t, eta_i, 0.1, in, on));
    // printf("%f\n", fr(i, o, n, eta_t, eta_i, 0.1, in, on));

    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum]
                  << std::endl;
    }

    if (argc != 3) {
        std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>"
                  << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];  // only bmp is allowed.

    // TODO: Main RayCasting Logic ------ OK
    // First, parse the scene using SceneParser.
    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.

    SceneParser sceneParser(inputFile.c_str());
    std::vector<AABB>::iterator it = sceneParser.a_aabb.begin();
    KDTree kdtree(sceneParser.a_aabb.begin(), sceneParser.a_aabb.size());
    Camera* camera = sceneParser.getCamera();

    Image image(camera->getWidth(), camera->getHeight());

    double dir[8][2] = {{-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5},
                       {-0.5, 0},    {0, -0.5},   {0.5, 0},    {0, 0.5}};

    // #pragma omp parallel
    //     {
    //         unsigned short Xi[3];

    // #pragma omp for schedule(guided)
    //         for (int pixel = 0; pixel < camera->getWidth() *
    //         camera->getHeight();
    //              ++pixel) {
    //             int x = pixel / camera->getHeight();
    //             int y = pixel % camera->getHeight();
    //             // Xi[2] = pixel * pixel * pixel;

    //             Vector3f color_sum = Vector3f::ZERO;

    //             Ray camRay = camera->generateRay(Vector2f(x, y));

    //             for (int iter = 0; iter < 100; ++iter) {
    //                 color_sum += radiance(Xi, sceneParser, kdtree, camRay, 1,
    //                 0, 1);
    //             }

    // #pragma omp critical
    //             image.SetPixel(x, y, color_sum / 100);
    //         }
    //     }

    int h = camera->getHeight();
    int w = camera->getWidth();
    int samps = 100 / 4;

    Vector3f r;

#pragma omp parallel for schedule(dynamic, 1) private(r)  // OpenMP
    for (int y = 0; y < camera->getHeight(); y++) {  // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4,
                100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y};
             x < camera->getWidth(); x++, r = Vector3f()) {  // Loop cols
            for (int sy = 0; sy < 2; sy++)        // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++) {  // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        double r1 = 2 * erand48(Xi),
                               dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi),
                               dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Ray camRay = camera->generateRay(
                            Vector2f(sx + dx + x, sy + dy + y));
                        r = r +
                            radiance(Xi, sceneParser, kdtree, camRay, 1e-4, 0, 1) *
                                (1. / samps);
                    }  // Camera rays are pushed ^^^^^ forward to start in
                       // interior
                    // c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))
                    // * .25;
                }
#pragma omp critical
            image.SetPixel(x, y, r / 4);
        }
    }

    image.SaveBMP(outputFile.c_str());

    std::cout << "Hello! Computer Graphics!" << endl;
    return 0;
}

// #include <math.h>    // smallpt, a Path Tracer by Kevin Beason, 2008
// #include <stdio.h>   //        Remove "-fopenmp" for g++ version < 4.2
// #include <stdlib.h>  // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt

// #include <pt.hpp>
// struct Vec {         // Usage: time ./smallpt 5000 && xv image.ppm
//     double x, y, z;  // position, also color (r,g,b)
//     Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
//         x = x_;
//         y = y_;
//         z = z_;
//     }
//     Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z +
//     b.z); } Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y,
//     z - b.z); } Vec operator*(double b) const { return Vec(x * b, y * b, z *
//     b); } Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z *
//     b.z); } Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z
//     * z)); } double dot(const Vec &b) const {
//         return x * b.x + y * b.y + z * b.z;
//     }  // cross:
//     Vec operator%(const Vec &b) {
//         return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
//     }
// };
// enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
// inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
// inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
// Vector3f toVector3f(const Vec &x) { return Vector3f(x.x, x.y, x.z); }
// Vec toVec(const Vector3f &x) { return Vec(x.x(), x.y(), x.z()); }
// int main(int argc, char *argv[]) {
//     for (int argNum = 1; argNum < argc; ++argNum) {
//         std::cout << "Argument " << argNum << " is: " << argv[argNum]
//                   << std::endl;
//     }

//     if (argc != 3) {
//         std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>"
//                   << endl;
//         return 1;
//     }
//     string inputFile = argv[1];
//     string outputFile = argv[2];  // only bmp is allowed.

//     // TODO: Main RayCasting Logic ------ OK
//     // First, parse the scene using SceneParser.
//     // Then loop over each pixel in the image, shooting a ray
//     // through that pixel and finding its intersection with
//     // the scene.  Write the color at the intersection to that
//     // pixel in your output image.

//     SceneParser sceneParser(inputFile.c_str());
//     std::vector<AABB>::iterator it = sceneParser.a_aabb.begin();
//     KDTree kdtree(sceneParser.a_aabb.begin(), sceneParser.a_aabb.size());
//     Camera *camera = sceneParser.getCamera();

//     Image image(camera->getWidth(), camera->getHeight());

//     int w = 1024, h = 768,
//         samps = argc == 2 ? atoi(argv[1]) / 4 : 1;              // # samples
//     Ray cam(Vector3f(50, 52, 295.6), Vector3f(0, -0.042612,
//     -1).normalized());  // cam pos, dir Vec cx = Vec(w * .5135 / h), cy = (cx
//     % toVec(cam.getDirection())).norm() * .5135, r,
//         *c = new Vec[w * h];
// #pragma omp parallel for schedule(dynamic, 1) private(r)  // OpenMP
//     for (int y = 0; y < h; y++) {  // Loop over image rows
//         fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4,
//                 100. * y / (h - 1));
//         for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w;
//              x++)  // Loop cols
//             for (int sy = 0, i = (h - y - 1) * w + x; sy < 2;
//                  sy++)  // 2x2 subpixel rows
//                 for (int sx = 0; sx < 2;
//                      sx++, r = Vec()) {  // 2x2 subpixel cols
//                     for (int s = 0; s < samps; s++) {
//                         double r1 = 2 * erand48(Xi),
//                                dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
//                         double r2 = 2 * erand48(Xi),
//                                dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
//                         Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
//                                 cy * (((sy + .5 + dy) / 2 + y) / h - .5) +
//                                 toVec(cam.getDirection());
//                         r = r + toVec(radiance(Xi, sceneParser, kdtree,
//                                                Ray(toVector3f(toVec(cam.getOrigin())
//                                                + d * 140),
//                                                    toVector3f(d.norm())),
//                                                1e-1, 0, 1) *
//                                       (1. / samps));
//                     }  // Camera rays are pushed ^^^^^ forward to start in
//                        // interior
//                     c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) *
//                     .25;
//                 }
//     }
//     FILE *f = fopen("image.ppm", "w");  // Write image to PPM file.
//     fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
//     for (int i = 0; i < w * h; i++)
//         fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
// }
