#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "lodepng.h"

using namespace std;

static constexpr int H = 200;
static constexpr int W = 200;

void step_dotted_line(unsigned char* const image, const int& width, const int& height,
                      const int& x0, const int& y0, const int& x1, const int& y1, const int& count, const int& color) {
    double const step = 1.0/(double)count;
    for (double i = 0; i < 1; i += step) {
        int const x = round((1.0 - i)*x0 + i*x1);
        int const y = round((1.0 - i)*y0 + i*y1);
        if (x < width && y < height) image[x + width * y] = color;
    }
}

void dotted_line(unsigned char* const image, const int& width, const int& height,
                 const int& x0, const int& y0, const int& x1, const int& y1, const int& color) {
    double const count = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
    double const step = 1.0/count;
    for (double i = 0; i < 1; i += step) {
        int const x = round((1.0 - i)*x0 + i*x1);
        int const y = round((1.0 - i)*y0 + i*y1);
        if (x < width && y < height) image[x + width * y] = color;
    }
}

void x_loop_line(unsigned char* const image, const int& width, const int& height,
                 int x0, int y0, int x1, int y1, const int& color) {
    bool xchange = false;
    if (abs(x0 - x1) < abs(y0 - y1)) {
        swap(x0, y0);
        swap(x1, y1);
        xchange = true;
    }
    if (x0 > x1) {
        swap(x0, x1);
        swap(y0, y1);
    }

    int y = y0;
    int dx = x1 - x0;
    int dy = abs(y1 - y0);
    int error = 0;
    int y_step = (y1 > y0) ? 1 : -1;

    for (int x = x0; x <= x1; x++) {
        int pixel_index;
        if (xchange) {
            if (y >= 0 && y < width && x >= 0 && x < height)
                pixel_index = y + x * width;
            else
                pixel_index = -1;
        }
        else {
            if (x >= 0 && x < width && y >= 0 && y < height)
                pixel_index = x + y * width;
            else
                pixel_index = -1;
        }

        if (pixel_index >= 0 && pixel_index < width * height)
            image[pixel_index] = color;

        error += dy;
        if (2 * error >= dx) {
            y += y_step;
            error -= dx;
        }
    }
}


void star1(unsigned char* image, const int& width, const int& height, const int& color,
           const int& radius, const int& points) {
    int const xc = width / 2;
    int const yc = height / 2;

    for (int i = 0; i < points; i++) {
        double const angle = 2 * M_PI * (double)i/points;
        int const x1 = xc + radius * cos(angle);
        int const y1 = yc + radius * sin(angle);
        step_dotted_line(image, width, height, xc, yc, x1, y1, 500, color);
    }
}

void star2(unsigned char* image, const int& width, const int& height, const int& color,
           const int& radius, const int& points) {
    int const xc = width / 2;
    int const yc = height / 2;

    for (int i = 0; i < points; i++) {
        double const angle = 2 * M_PI * (double)i/points;
        int const x1 = xc + radius * cos(angle);
        int const y1 = yc + radius * sin(angle);
        dotted_line(image, width, height, xc, yc, x1, y1, color);
    }
}

void star3(unsigned char* image, const int& width, const int& height, const int& color,
           const int& radius, const int& points) {
    int const xc = width / 2;
    int const yc = height / 2;

    for (int i = 0; i < points; i++) {
        double const angle = 2 * M_PI * (double)i/points;
        int const x1 = xc + radius * cos(angle);
        int const y1 = yc + radius * sin(angle);
        x_loop_line(image, width, height, xc, yc, x1, y1, color);
    }
}

struct Vertex {
    double x, y, z;
};

struct Face {
    Vertex v1, v2, v3;
};

void readObj(const string& filename, vector<Vertex>& vertices, vector<Face>& faces) {
    ifstream file(filename);
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string type;
        iss >> type;

        if (type == "v") {
            Vertex v{};
            iss >> v.x >> v.y >> v.z;
            vertices.push_back(v);
        }
        if (type == "f") {
            Face f{};
            int i;
            iss >> i;
            f.v1 = vertices[i - 1];
            iss.ignore(INT_MAX, ' ');
            iss >> i;
            f.v2 = vertices[i - 1];
            iss.ignore(INT_MAX, ' ');
            iss >> i;
            f.v3 = vertices[i - 1];
            faces.push_back(f);
        }
    }
}

void drawObj(unsigned char* image, const int& width, const int& height,
             const vector<Vertex>& vertices, const vector<Face>& faces) {
    for (Vertex v : vertices) {
        int const x = 8000 * v.x + 500;
        int const y = 8000 * v.y + 100;
        if (x > 0 && y > 0 && x < width && y < height) image[y * width + x] = 255;
    }
    for (Face f : faces) {
        int const x1 = 8000 * f.v1.x + 500;
        int const y1 = 8000 * f.v1.y + 100;
        int const x2 = 8000 * f.v2.x + 500;
        int const y2 = 8000 * f.v2.y + 100;
        int const x3 = 8000 * f.v3.x + 500;
        int const y3 = 8000 * f.v3.y + 100;
        x_loop_line(image, width, height, x1, y1, x2, y2, 255);
        x_loop_line(image, width, height, x2, y2, x3, y3, 255);
        x_loop_line(image, width, height, x3, y3, x1, y1, 255);
    }
}

int main() {
    // Генерация одноканальной картинки
    unsigned char* greyimage = new unsigned char[W * H];

    // Чёрная картинка
    for (int i = 0; i < H * W; i++) greyimage[i] = 0;
    lodepng::encode("blackkartinka.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated black." << std::endl;

    // Белая картинка
    for (int i = 0; i < H * W; i++) greyimage[i] = 255;
    lodepng::encode("whitekartinka.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated white." << std::endl;

    // Произвольная картинка
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            greyimage[x + W * y] = (x + y) / 2 < 256 ? (x + y) / 2 : 255;

    lodepng::encode("gradkartinka.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated gradient." << std::endl;

    // Произвольная картинка #2
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            greyimage[x + x * y] = (x + y) % 256;

    lodepng::encode("kartinka?.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated something..." << std::endl;

    // Звезда 1
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star1(greyimage, W, H, 255, 95, 13);
    lodepng::encode("star1.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated star 1." << std::endl;

    // Звезда 2
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star2(greyimage, W, H, 255, 95, 13);
    lodepng::encode("star2.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated star 2." << std::endl;

    // Звезда 3
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star3(greyimage, W, H, 255, 95, 13);
    lodepng::encode("star3.png", greyimage, W, H, LCT_GREY, 8);
    std::cout << "Generated star 2." << std::endl;

    delete greyimage;

    // Генерация трёхканальной картинки
    unsigned char* rgbimage = new unsigned char[W * H * 3];

    // Красная картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 255; rgbimage[i + 1] = 0; rgbimage[i + 2] = 0;
    }
    lodepng::encode("redkartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated red." << std::endl;

    // Зелёная картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 0; rgbimage[i + 1] = 255; rgbimage[i + 2] = 0;
    }
    lodepng::encode("greenkartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated green." << std::endl;

    // Синяя картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 0; rgbimage[i + 1] = 0; rgbimage[i + 2] = 255;
    }
    lodepng::encode("bluekartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated blue." << std::endl;

    delete rgbimage;

    // 3D-Модель
    unsigned char* modelimage = new unsigned char[1000 * 1000];
    vector<Vertex> vertices;
    vector<Face> faces;
    string filename = "model_1.obj";

    readObj(filename, vertices, faces);
    drawObj(modelimage, 1000, 1000, vertices, faces);
    lodepng::encode("model_1.png", modelimage, 1000, 1000, LCT_GREY, 8);

    delete modelimage;

    return 0;
}
