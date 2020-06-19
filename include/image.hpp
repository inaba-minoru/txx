#ifndef IMAGE_H
#define IMAGE_H

#include <vecmath.h>

#include <cassert>

// Simple image class
class Image {
   public:
    Image(int w, int h) {
        width = w;
        height = h;
        data = new Vector3f[width * height];
    }

    ~Image() { delete[] data; }

    int Width() const { return width; }

    int Height() const { return height; }

    const Vector3f &GetPixel(int x, int y) const {
        assert(x >= 0 && x < width);
        assert(y >= 0 && y < height);
        return data[y * width + x];
    }

    void SetAllPixels(const Vector3f &color) {
        for (int i = 0; i < width * height; ++i) {
            data[i] = color;
        }
    }

    void SetPixel(int x, int y, const Vector3f &color) {
        assert(x >= 0 && x < width);
        assert(y >= 0 && y < height);
        data[y * width + x] = color;
    }

    static Image *LoadPPM(const char *filename);

    void SavePPM(const char *filename) const;

    static Image *LoadTGA(const char *filename);

    void SaveTGA(const char *filename) const;

    int SaveBMP(const char *filename);

    void SaveImage(const char *filename);

    void Gauss() {
        int N = 3;
        int Ndiv2 = N / 2;
        double sigma = 0.8;
        double sigma22 = 2 * sigma * sigma;
        double gauss[N][N], sum;
        for (int i = -Ndiv2; i <= Ndiv2; i++)
            for (int j = -Ndiv2; j <= Ndiv2; j++)
                gauss[i + Ndiv2][j + Ndiv2] = exp(-i * i - j * j / sigma22);

        for (size_t x = 0; x < width; x++) {
            for (size_t y = 0; y < height; y++) {
                Vector3f temp = 0;
                double sum = 0;
                for (int i = -Ndiv2; i <= Ndiv2; i++)
                    for (int j = -Ndiv2; j <= Ndiv2; j++) {
                        int tx = x + i;
                        int ty = y + j;
                        int ti = i + Ndiv2;
                        int tj = j + Ndiv2;
                        if (tx >= 0 && tx < width && ty >= 0 && ty < height) {
                            temp += data[ty * width + tx] * gauss[ti][tj];
                            sum += gauss[ti][tj];
                        }
                    }
                data[y * width + x] = temp / sum;
            }
        }
    }

   private:
    int width;
    int height;
    Vector3f *data;
};

#endif  // IMAGE_H
