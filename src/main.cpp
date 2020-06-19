#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
// #include <kdtree.hpp>
#include <algorithm>
#include <pt.hpp>
#include <string>

#include "camera.hpp"
#include "group.hpp"
#include "image.hpp"
#include "light.hpp"
#include "scene_parser.hpp"

// using namespace std;

// #include <ggx.hpp>

double clamp(const double& x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
Vector3f clamp(const Vector3f& x) {
    return Vector3f(clamp(x.x()), clamp(x.y()), clamp(x.z()));
}
Vector3f gamma(const Vector3f& color) {
    return Vector3f(pow(clamp(color.x()), 1 / 2.2),
                    pow(clamp(color.y()), 1 / 2.2),
                    pow(clamp(color.z()), 1 / 2.2));
}

int main(int argc, char* argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum]
                  << std::endl;
    }

    if (argc != 3) {
        std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>"
                  << std::endl;
        return 1;
    }
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];  // only bmp is allowed.

    // TODO: Main RayCasting Logic ------ OK
    // First, parse the scene using SceneParser.
    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.

    SceneParser sceneParser(inputFile.c_str());
    Camera* camera = sceneParser.getCamera();

    Image image(camera->getWidth(), camera->getHeight());

    sceneParser.getGroup()->init();

    int N = 10;
    int N2 = N * N;
    int samps_per_grid = 1;
    int samps = N2 * samps_per_grid;

    // int N = 9;
    // int samps = 100 / N;

#pragma omp parallel
    {
        unsigned short Xi[3];

#pragma omp for schedule(dynamic, 1)
        for (int pixel = 0; pixel < camera->getWidth() * camera->getHeight();
             pixel++) {
            fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps,
                    100. * pixel / (camera->getHeight() * camera->getWidth()));

            int x = pixel / camera->getHeight();
            int y = pixel % camera->getHeight();

            Xi[2] = pixel * pixel * pixel;

            Vector3f sum;

            for (size_t xx = 0; xx < N; xx++) {
                for (size_t yy = 0; yy < N; yy++) {
                    double sx = x + (xx + erand48(Xi)) / N;
                    double sy = y + (yy + erand48(Xi)) / N;
                    // Vector3f temp;
                    Ray cam_ray = camera->generateRay(Xi, Vector2f(sx, sy));
                    for (size_t i = 0; i < samps_per_grid; i++) {
                        sum +=
                            radiance(Xi, sceneParser, cam_ray, 1e-4, 0, 1, 1) /
                            samps_per_grid / N2;
                    }
                    // sum += temp / N2;
                }
            }
#pragma omp critical
            image.SetPixel(x, y, gamma(sum));

            //             bool rooks[9]{};
            //             Vector3f sum(0);
            //             for (int i = 0, j, idx; i < N; i++) {
            //                 //
            //                 idx = erand48(Xi) * (N - i) + 1;
            //                 for (j = 0; j < N; j++) {
            //                     idx -= !rooks[j];
            //                     if (!idx) break;
            //                 }
            //                 rooks[j] = 1;
            //                 //
            //                 assert(j < N);
            //                 Ray camRay = camera->generateRay(
            //                     Vector2f(x - (N / 2. + i) / N, y - (N / 2. +
            //                     j) / N));
            //                 Vector3f temp(0);
            //                 for (int k = 0; k < samps; ++k)
            //                     temp +=
            //                         radiance(Xi, sceneParser, camRay, 1e-4,
            //                         0, 1) / samps;
            //                 sum += clamp(temp) / N;
            //             }
            // #pragma omp critical
            //             image.SetPixel(x, y, gamma(sum));

            // Vector3f color_sum = Vector3f::ZERO;

            // Ray camRay = camera->generateRay(Vector2f(x, y));

            // for (int s = 0; s < samps; s++)
            //     color_sum +=
            //         radiance(Xi, sceneParser, camRay, 1e-4, 0, 1) / samps;

            // #pragma omp critical
            //             image.SetPixel(x, y, gamma(color_sum));
        }
    }

    //     int h = camera->getHeight();
    //     int w = camera->getWidth();
    //     int samps = 100 / 4;

    //     Vector3f r, t;

    // #pragma omp parallel for schedule(dynamic, 1) private(r, t)  // OpenMP
    //     for (int y = 0; y < camera->getHeight(); y++) {  // Loop over image
    //     rows
    //         fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4,
    //                 100. * y / (h - 1));
    //         for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y};
    //              x < camera->getWidth(); x++, r = Vector3f()) {  // Loop cols
    //             for (int sy = 0; sy < 2; sy++)  // 2x2 subpixel rows
    //                 for (int sx = 0; sx < 2;
    //                      sx++, t = Vector3f()) {  // 2x2 subpixel cols
    //                     for (int s = 0; s < samps; s++) {
    //                         double r1 = 2 * erand48(Xi),
    //                                dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 -
    //                                r1);
    //                         double r2 = 2 * erand48(Xi),
    //                                dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 -
    //                                r2);
    //                         Ray camRay = camera->generateRay(
    //                             Vector2f(sx + dx + x, sy + dy + y));
    //                         t += radiance(Xi, sceneParser, camRay, 1e-4, 0,
    //                         1) *
    //                              (1. / samps);
    //                     }  // Camera rays are pushed ^^^^^ forward to start
    //                     in r += clamp(t) / 4;
    //                 }
    // #pragma omp critical
    //             image.SetPixel(x, y, gamma(r));
    //         }
    //     }

    // image.Gauss();
    image.SaveBMP(outputFile.c_str());

    std::cout << "Hello! Computer Graphics!" << std::endl;
    return 0;
}