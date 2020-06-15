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

int main(int argc, char *argv[]) {
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
    Camera *camera = sceneParser.getCamera();

    Image image(camera->getWidth(), camera->getHeight());

    float dir[8][2] = {{-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5},
                       {-0.5, 0},    {0, -0.5},   {0.5, 0},    {0, 0.5}};

#pragma omp parallel for schedule(guided)
    for (int pixel = 0; pixel < camera->getWidth() * camera->getHeight();
         ++pixel) {
        int x = pixel / camera->getHeight();
        int y = pixel % camera->getHeight();

        Vector3f color_sum = Vector3f::ZERO;

        Ray camRay = camera->generateRay(Vector2f(x, y));

        for (int iter = 0; iter < 1000; ++iter)
            color_sum += radiance(sceneParser, kdtree, camRay, 1e-3, 0, 1);

#pragma omp critical
        image.SetPixel(x, y, color_sum / 1000);
    }
    image.SaveBMP(outputFile.c_str());

    std::cout << "Hello! Computer Graphics!" << endl;
    return 0;
}
