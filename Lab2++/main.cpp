#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "lodepng.h"
#include "giflib/gif_lib.h"

using namespace std;

static constexpr int H = 200;
static constexpr int W = 200;

static constexpr int scale = 1800;
static constexpr int offsetX = 100;
static constexpr int offsetY = 10;

struct GIF {
    vector<unsigned char*> frames;
    vector<int> delays;
    int width, height;

    GIF(const int& width, const int& height) : width(width), height(height) {}

    void add_frame(unsigned char* frame, const int& width, const int& height) {
        int delay_ms = 100;
        if (width == this->width && height == this->height) {
            frames.push_back(frame);
            delays.push_back(delay_ms / 10);
        }
    }

    void clear() {
        for (unsigned char* frame : this->frames)
            delete[] frame;
        frames.clear();
        delays.clear();
    }

    int get_frame_count() const {
        return frames.size();
    }

    bool save(const string& filename, const bool& LZW) {
        int error = 0;
        GifFileType* gif = EGifOpenFileName(filename.c_str(), false, &error);
        if (!gif) {
            cerr << "Failed to open file for writing: " << GifErrorString(error) << endl;
            return false;
        }

        // Создаем цветовую карту
        int color_size = LZW ? 16 : 256;
        ColorMapObject* colorMap = GifMakeMapObject(color_size, NULL);
        if (!colorMap) {
            cerr << "Failed to create color map" << endl;
            EGifCloseFile(gif, &error);
            return false;
        }

        // Заполняем цветовую карту градациями серого
        for (int i = 0; i < color_size; i++) {
            int value = LZW ? (i * 255 / 15) : i;
            colorMap->Colors[i].Red = value;
            colorMap->Colors[i].Green = value;
            colorMap->Colors[i].Blue = value;
        }

        // Устанавливаем глобальную цветовую карту
        if (EGifPutScreenDesc(gif, width, height, LZW ? 4 : 8, 0, colorMap) != GIF_OK) {
            cerr << "Failed to set screen description" << endl;
            GifFreeMapObject(colorMap);
            EGifCloseFile(gif, &error);
            return false;
        }

        GifFreeMapObject(colorMap);

        // Netscape Application Extension
        unsigned char netscapeExt[] = "NETSCAPE2.0";
        unsigned char loopData[3] = { 0x01, 0x00, 0x00 }; // loop forever

        // Application Extension block
        if (EGifPutExtensionLeader(gif, APPLICATION_EXT_FUNC_CODE) != GIF_OK) {
            cerr << "Failed to write application extension leader" << endl;
        }

        // Application Identifier and Authentication Code (11 bytes)
        if (EGifPutExtensionBlock(gif, 11, netscapeExt) != GIF_OK) {
            cerr << "Failed to write Netscape extension block" << endl;
        }

        // Loop data (3 bytes)
        if (EGifPutExtensionBlock(gif, 3, loopData) != GIF_OK) {
            cerr << "Failed to write loop data" << endl;
        }

        // Terminate extension
        if (EGifPutExtensionTrailer(gif) != GIF_OK) {
            cerr << "Failed to write extension trailer" << endl;
        }

        // Записываем кадры
        for (size_t i = 0; i < frames.size(); i++) {
            // Graphic Control Extension
            unsigned char gce[4] = {
                0x04, // disposal method
                static_cast<unsigned char>(delays[i] & 0xFF),
                static_cast<unsigned char>((delays[i] >> 8) & 0xFF),
                0x00  // no transparent color
            };

            if (EGifPutExtension(gif, 0xF9, 4, gce) != GIF_OK) {
                cerr << "Failed to write GCE for frame " << i << endl;
                continue;
            }

            // Image descriptor - используем глобальную цветовую карту
            if (EGifPutImageDesc(gif, 0, 0, width, height, false, NULL) != GIF_OK) {
                cerr << "Failed to set image description for frame " << i << endl;
                continue;
            }

            // Записываем данные построчно
            const unsigned char* frame_data = frames[i];
            for (int y = 0; y < height; y++) {
                if (EGifPutLine(gif, const_cast<GifByteType*>(&frame_data[y * width]), width) != GIF_OK) {
                    cerr << "Failed to write line " << y << " of frame " << i << endl;
                    break;
                }
            }
        }

        if (EGifCloseFile(gif, &error) != GIF_OK) {
            cerr << "Failed to close GIF file. Error: " << GifErrorString(error) << endl;
            return false;
        }

        cout << "GIF successfully saved: " << filename << " (" << frames.size() << " frames)" << endl;
        return true;
    }

    ~GIF() {
        clear();
    }
};

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

void x_loop_line(unsigned char* const image, const int& width, const int& height,
                     int x0, int y0, int x1, int y1,
                     const int& r, const int& g, const int& b) {
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
                pixel_index = (y + x * width) * 3;
            else
                pixel_index = -1;
        }
        else {
            if (x >= 0 && x < width && y >= 0 && y < height)
                pixel_index = (x + y * width) * 3;
            else
                pixel_index = -1;
        }

        if (pixel_index >= 0 && pixel_index < width * height * 3) {
            image[pixel_index] = r;
            image[pixel_index + 1] = g;
            image[pixel_index + 2] = b;
        }

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

struct Point {
    double x, y;
};

struct Vertex {
    double x, y, z;
};

struct Face {
    int v1, v2, v3;
};

struct Vector3D {
    double x, y, z;

    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    Vector3D normalize() const {
        double const len = length();
        if (len > 0) return {x/len, y/len, z/len};
        return {0, 0, 0};
    }
};

struct Matrix4 {
    double m[4][4] = {0};

    Matrix4() {
        // Единичная матрица по умолчанию
        for (int i = 0; i < 4; i++) m[i][i] = 1.0;
    }
};

Matrix4 rotationX(double angle) {
    Matrix4 mat;
    double c = cos(angle);
    double s = sin(angle);
    mat.m[1][1] = c;
    mat.m[1][2] = s;
    mat.m[2][1] = -s;
    mat.m[2][2] = c;
    return mat;
}

Matrix4 rotationY(double angle) {
    Matrix4 mat;
    double c = cos(angle);
    double s = sin(angle);
    mat.m[0][0] = c;
    mat.m[0][2] = -s;
    mat.m[2][0] = s;
    mat.m[2][2] = c;
    return mat;
}

Matrix4 rotationZ(double angle) {
    Matrix4 mat;
    double c = cos(angle);
    double s = sin(angle);
    mat.m[0][0] = c;
    mat.m[0][1] = s;
    mat.m[1][0] = -s;
    mat.m[1][1] = c;
    return mat;
}

static Vector3D LDIR = Vector3D(0, 0, 1);

void read_obj(const string& filename, vector<Vertex>& vertices, vector<Face>& faces) {
    ifstream file(filename);
    string line;
    vertices.clear();
    faces.clear();
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
            f.v1 = i - 1;
            iss.ignore(INT_MAX, ' ');
            iss >> i;
            f.v2 = i - 1;
            iss.ignore(INT_MAX, ' ');
            iss >> i;
            f.v3 = i - 1;
            faces.push_back(f);
        }
    }
}

array<double, 3> get_barycentric_coordinates(int x, int y, double x0, double y0, double x1, double y1,
    double x2, double y2) {
    double lambda0 = ((x - x2) * (y1 - y2) - (x1 - x2) * (y - y2)) / ((x0 - x2) * (y1 - y2) - (x1 - x2) * (y0 - y2));
    double lambda1 = ((x0 - x2) * (y - y2) - (x - x2) * (y0 - y2)) / ((x0 - x2) * (y1 - y2) - (x1 - x2) * (y0 - y2));
    double lambda2 = 1.0 - lambda0 - lambda1;

    return {lambda0, lambda1, lambda2};
}

Vector3D calculateTriangleNormal(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
    Vector3D edge1 = Vector3D(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
    Vector3D edge2 = Vector3D(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);

    Vector3D normal = Vector3D(
        edge1.y * edge2.z - edge1.z * edge2.y,
        edge1.z * edge2.x - edge1.x * edge2.z,
        edge1.x * edge2.y - edge1.y * edge2.x
    );

    return normal.normalize();
}

void shift_vertex(Vertex& vertex, const Vector3D& shift) {
    vertex.x += shift.x;
    vertex.y += shift.y;
    vertex.z += shift.z;
}

void shift_vertex(Vector3D& vector, const Vector3D& shift) {
    vector.x += shift.x;
    vector.y += shift.y;
    vector.z += shift.z;
}

Matrix4 multiply_matrices(const Matrix4& a, const Matrix4& b) {
    Matrix4 result;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.m[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                result.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return result;
}

void multiply_matrix_vertex(Vertex& vertex, const Matrix4& mat) {
    double x = vertex.x;
    double y = vertex.y;
    double z = vertex.z;

    vertex.x = x * mat.m[0][0] + y * mat.m[1][0] + z * mat.m[2][0] + mat.m[3][0];
    vertex.y = x * mat.m[0][1] + y * mat.m[1][1] + z * mat.m[2][1] + mat.m[3][1];
    vertex.z = x * mat.m[0][2] + y * mat.m[1][2] + z * mat.m[2][2] + mat.m[3][2];
    double w = x * mat.m[0][3] + y * mat.m[1][3] + z * mat.m[2][3] + mat.m[3][3];

    if (w != 0.0) {
        vertex.x /= w; vertex.y /= w; vertex.z /= w;
    }
}

void multiply_matrix_vertex( Vector3D& vector, const Matrix4& mat) {
    double x = vector.x;
    double y = vector.y;
    double z = vector.z;

    vector.x = x * mat.m[0][0] + y * mat.m[1][0] + z * mat.m[2][0] + mat.m[3][0];
    vector.y = x * mat.m[0][1] + y * mat.m[1][1] + z * mat.m[2][1] + mat.m[3][1];
    vector.z = x * mat.m[0][2] + y * mat.m[1][2] + z * mat.m[2][2] + mat.m[3][2];
    double w = x * mat.m[0][3] + y * mat.m[1][3] + z * mat.m[2][3] + mat.m[3][3];

    if (w != 0.0) {
        vector.x /= w; vector.y /= w; vector.z /= w;
    }
}

void transform(vector<Vertex>& vertices, const Vector3D& shift, const Matrix4& rotation) {
    for (Vertex& v : vertices) {
        multiply_matrix_vertex(v, rotation);
        shift_vertex(v, shift);
    }
}

void draw_triangle(unsigned char* image, const int& width, const int& height, const vector<Vertex>& vertices,
    const Face& face) {
    Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);

    double const x0 = scale * vertices[face.v1].x + offsetX;
    double const y0 = scale * vertices[face.v1].y + offsetY;
    double const x1 = scale * vertices[face.v2].x + offsetX;
    double const y1 = scale * vertices[face.v2].y + offsetY;
    double const x2 = scale * vertices[face.v3].x + offsetX;
    double const y2 = scale * vertices[face.v3].y + offsetY;

    int x_min = max(0, (int)floor(min(min(x0, x1), x2)));
    int x_max = min(width-1, (int)ceil(max(max(x0, x1), x2)));
    int y_min = max(0, (int)floor(min(min(y0, y1), y2)));
    int y_max = min(height-1, (int)ceil(max(max(y0, y1), y2)));

    int red = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length());
    int green = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length());
    int blue = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length());

    for (int i = y_min; i <= y_max; i++) {
        for (int j = x_min; j <= x_max; j++) {
            array<double, 3> barycentric = get_barycentric_coordinates(j, i, x0, y0, x1, y1, x2, y2);
            if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                if (i >= 0 && i < height && j >= 0 && j < width) {
                    image[(j + i * width) * 3] = red;
                    image[(j + i * width) * 3 + 1] = green;
                    image[(j + i * width) * 3 + 2] = blue;
                }
            }
        }
    }
}

void draw_triangle(unsigned char* image, const int& width, const int& height, const vector<Vertex>& vertices,
    const Face& face, unsigned char r, unsigned char g, unsigned char b) {

    double const x0 = scale * vertices[face.v1].x + offsetX;
    double const y0 = scale * vertices[face.v1].y + offsetY;
    double const x1 = scale * vertices[face.v2].x + offsetX;
    double const y1 = scale * vertices[face.v2].y + offsetY;
    double const x2 = scale * vertices[face.v3].x + offsetX;
    double const y2 = scale * vertices[face.v3].y + offsetY;

    int x_min = max(0, (int)floor(min(min(x0, x1), x2)));
    int x_max = min(width-1, (int)ceil(max(max(x0, x1), x2)));
    int y_min = max(0, (int)floor(min(min(y0, y1), y2)));
    int y_max = min(height-1, (int)ceil(max(max(y0, y1), y2)));

    for (int i = y_min; i <= y_max; i++) {
        for (int j = x_min; j <= x_max; j++) {
            array<double, 3> barycentric = get_barycentric_coordinates(j, i, x0, y0, x1, y1, x2, y2);
            if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                if (i >= 0 && i < height && j >= 0 && j < width) {
                    image[(j + i * width) * 3] = r;
                    image[(j + i * width) * 3 + 1] = g;
                    image[(j + i * width) * 3 + 2] = b;
                }
            }
        }
    }
}

void draw_triangle_zbuffer(unsigned char* image, float* zbuffer, const int& width, const int& height,
                               const vector<Vertex>& vertices, const Face& face) {
    Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);

    double const x0 = scale * vertices[face.v1].x + offsetX;
    double const y0 = scale * vertices[face.v1].y + offsetY;
    double const z0 = vertices[face.v1].z;  // Z-координаты вершин

    double const x1 = scale * vertices[face.v2].x + offsetX;
    double const y1 = scale * vertices[face.v2].y + offsetY;
    double const z1 = vertices[face.v2].z;

    double const x2 = scale * vertices[face.v3].x + offsetX;
    double const y2 = scale * vertices[face.v3].y + offsetY;
    double const z2 = vertices[face.v3].z;

    int x_min = max(0, (int)floor(min(min(x0, x1), x2)));
    int x_max = min(width-1, (int)ceil(max(max(x0, x1), x2)));
    int y_min = max(0, (int)floor(min(min(y0, y1), y2)));
    int y_max = min(height-1, (int)ceil(max(max(y0, y1), y2)));

    // Освещение
    double dot_light = - (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z);
    dot_light = max(0.0, dot_light);
    int intensity = (int)(250 * dot_light);
    intensity = max(0, min(255, intensity));

    for (int i = y_min; i <= y_max; i++) {
        for (int j = x_min; j <= x_max; j++) {
            array<double, 3> barycentric = get_barycentric_coordinates(j, i, x0, y0, x1, y1, x2, y2);

            if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                // Вычисляем z-координату через барицентрические координаты
                double z_interpolated = barycentric[0] * z0 + barycentric[1] * z1 + barycentric[2] * z2;

                int pixel_index = j + i * width;

                // Z-буфер: рисуем только если пиксель ближе к камере
                if (z_interpolated < zbuffer[pixel_index]) {
                    zbuffer[pixel_index] = z_interpolated; // Обновляем z-буфер

                    int img_index = pixel_index * 3;
                    image[img_index] = intensity;
                    image[img_index + 1] = intensity;
                    image[img_index + 2] = intensity;
                }
            }
        }
    }
}

void draw_triangle_zbuffer(unsigned char* image, float* zbuffer, const int& width, const int& height,
                               const vector<Vertex>& vertices, const Face& face,
                               const unsigned char& r, const unsigned char& g, const unsigned char& b) {
    Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);

    double const x0 = scale * vertices[face.v1].x + offsetX;
    double const y0 = scale * vertices[face.v1].y + offsetY;
    double const z0 = vertices[face.v1].z;  // Z-координаты вершин

    double const x1 = scale * vertices[face.v2].x + offsetX;
    double const y1 = scale * vertices[face.v2].y + offsetY;
    double const z1 = vertices[face.v2].z;

    double const x2 = scale * vertices[face.v3].x + offsetX;
    double const y2 = scale * vertices[face.v3].y + offsetY;
    double const z2 = vertices[face.v3].z;

    int x_min = max(0, (int)floor(min(min(x0, x1), x2)));
    int x_max = min(width-1, (int)ceil(max(max(x0, x1), x2)));
    int y_min = max(0, (int)floor(min(min(y0, y1), y2)));
    int y_max = min(height-1, (int)ceil(max(max(y0, y1), y2)));

    // Освещение
    double dot_light = - (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z);
    dot_light = max(0.0, dot_light);
    int intensity = (int)(250 * dot_light);
    intensity = max(0, min(255, intensity));

    for (int i = y_min; i <= y_max; i++) {
        for (int j = x_min; j <= x_max; j++) {
            array<double, 3> barycentric = get_barycentric_coordinates(j, i, x0, y0, x1, y1, x2, y2);

            if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                // Вычисляем z-координату через барицентрические координаты
                double z_interpolated = barycentric[0] * z0 + barycentric[1] * z1 + barycentric[2] * z2;

                int pixel_index = j + i * width;

                // Z-буфер: рисуем только если пиксель ближе к камере
                if (z_interpolated < zbuffer[pixel_index]) {
                    zbuffer[pixel_index] = z_interpolated; // Обновляем z-буфер

                    int img_index = pixel_index * 3;
                    image[img_index] = r * intensity;
                    image[img_index + 1] = g * intensity;
                    image[img_index + 2] = b * intensity;
                }
            }
        }
    }
}

void draw_wireframe_obj(unsigned char* image, const int& width, const int& height,
             const vector<Vertex>& vertices, const vector<Face>& faces) {
    for (Vertex v : vertices) {
        int const x = scale * v.x + 500;
        int const y = scale * v.y + offsetY;
        if (x > 0 && y > 0 && x < width && y < height) image[y * width + x] = 255;
    }
    for (const auto& face : faces) {
        int const x1 = scale * vertices[face.v1].x + offsetX;
        int const y1 = scale * vertices[face.v1].y + offsetY;
        int const x2 = scale * vertices[face.v2].x + offsetX;
        int const y2 = scale * vertices[face.v2].y + offsetY;
        int const x3 = scale * vertices[face.v3].x + offsetX;
        int const y3 = scale * vertices[face.v3].y + offsetY;
        x_loop_line(image, width, height, x1, y1, x2, y2, 255);
        x_loop_line(image, width, height, x2, y2, x3, y3, 255);
        x_loop_line(image, width, height, x3, y3, x1, y1, 255);
    }
}

void draw_solid_obj(unsigned char* image, const int& width, const int& height, const vector<Vertex>& vertices,
    const vector<Face>& faces) {
    for (const auto& face : faces) {
        Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);
        if ((normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length()) < 0)
            draw_triangle(image, width, height, vertices, face);

        // int const x1 = 8000 * f.v1.x + 500;
        // int const y1 = 8000 * f.v1.y + 100;
        // int const x2 = 8000 * f.v2.x + 500;
        // int const y2 = 8000 * f.v2.y + 100;
        // int const x3 = 8000 * f.v3.x + 500;
        // int const y3 = 8000 * f.v3.y + 100;
        // x_loop_line(image, width, height, x1, y1, x2, y2, 255);
        // x_loop_line(image, width, height, x2, y2, x3, y3, 255);
        // x_loop_line(image, width, height, x3, y3, x1, y1, 255);
    }
}

void draw_solid_obj_zbuffer(unsigned char* image, const int& width, const int& height, const vector<Vertex>& vertices,
    const vector<Face>& faces) {

    float* zbuffer = new float[width * height];
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = 1000.0f; // Большое значение - всё далеко
    }

    for (const auto& face : faces) {
        Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);
        if ((normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length()) < 0)
            draw_triangle_zbuffer(image, zbuffer, width, height, vertices, face);

        // int const x1 = 8000 * f.v1.x + 500;
        // int const y1 = 8000 * f.v1.y + 100;
        // int const x2 = 8000 * f.v2.x + 500;
        // int const y2 = 8000 * f.v2.y + 100;
        // int const x3 = 8000 * f.v3.x + 500;
        // int const y3 = 8000 * f.v3.y + 100;
        // x_loop_line(image, width, height, x1, y1, x2, y2, 255);
        // x_loop_line(image, width, height, x2, y2, x3, y3, 255);
        // x_loop_line(image, width, height, x3, y3, x1, y1, 255);
    }

    delete[] zbuffer;
}

void draw_solid_obj_zbuffer(unsigned char* image, const int& width, const int& height, const vector<Vertex>& vertices,
    const vector<Face>& faces, const unsigned int& r, const unsigned int& g, const unsigned int& b) {

    float* zbuffer = new float[width * height];
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = 1000.0f; // Большое значение - всё далеко
    }

    for (const auto& face : faces) {
        Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);
        if ((normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z) / (normal.length() * LDIR.length()) < 0)
            draw_triangle_zbuffer(image, zbuffer, width, height, vertices, face, r, g, b);

        // int const x1 = 8000 * f.v1.x + 500;
        // int const y1 = 8000 * f.v1.y + 100;
        // int const x2 = 8000 * f.v2.x + 500;
        // int const y2 = 8000 * f.v2.y + 100;
        // int const x3 = 8000 * f.v3.x + 500;
        // int const y3 = 8000 * f.v3.y + 100;
        // x_loop_line(image, width, height, x1, y1, x2, y2, 255);
        // x_loop_line(image, width, height, x2, y2, x3, y3, 255);
        // x_loop_line(image, width, height, x3, y3, x1, y1, 255);
    }

    delete[] zbuffer;
}

void draw_normals(unsigned char* image, const int& width, const int& height,
                  const vector<Vertex>& vertices, const vector<Face>& faces) {

    for (const auto& face : faces) {
        Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);

        // Центр треугольника
        double cx = (vertices[face.v1].x + vertices[face.v2].x + vertices[face.v3].x) / 3;
        double cy = (vertices[face.v1].y + vertices[face.v2].y + vertices[face.v3].y) / 3;

        // Отрисовка нормали в 3D пространстве
        int x1 = scale * cx + offsetX;
        int y1 = scale * cy + offsetY;
        int x2 = scale * (cx + normal.x * 0.001) + offsetX;
        int y2 = scale * (cy + normal.y * 0.001) + offsetY;

        // Цвет нормали в зависимости от направления
        // int color = 128;
        // if (normal.z > 0.5) color = 255; // Передние грани - белые
        // else if (normal.z < -0.5) color = 64; // Задние грани - тёмные

        int red = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
                    / (normal.length() * LDIR.length());
        int green = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
                    / (normal.length() * LDIR.length());
        int blue = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
                    / (normal.length() * LDIR.length());

        x_loop_line(image, width, height, x1, y1, x2, y2, red, green, blue);
    }
}

void draw_normals(unsigned char* image, const int& width, const int& height,
                  const vector<Vertex>& vertices, const vector<Face>& faces,
                  const unsigned int& r, const unsigned int& g, const unsigned int& b) {

    float* zbuffer = new float[width * height];
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = 1000.0f;
    }

    for (const auto& face : faces) {
        Vector3D normal = calculateTriangleNormal(vertices[face.v1], vertices[face.v2], vertices[face.v3]);

        // Центр треугольника
        double cx = (vertices[face.v1].x + vertices[face.v2].x + vertices[face.v3].x) / 3;
        double cy = (vertices[face.v1].y + vertices[face.v2].y + vertices[face.v3].y) / 3;

        // Отрисовка нормали в 3D пространстве
        int x1 = scale * cx + offsetX;
        int y1 = scale * cy + offsetY;
        int x2 = scale * (cx + normal.x * 0.001) + offsetX;
        int y2 = scale * (cy + normal.y * 0.001) + offsetY;

        // Цвет нормали в зависимости от направления
        // int color = 128;
        // if (normal.z > 0.5) color = 255; // Передние грани - белые
        // else if (normal.z < -0.5) color = 64; // Задние грани - тёмные

        // int red = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
        //          / (normal.length() * LDIR.length());
        // int green = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
        //          / (normal.length() * LDIR.length());
        // int blue = -250 * (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z)
        //          / (normal.length() * LDIR.length());

        double dot_light = - (normal.x * LDIR.x + normal.y * LDIR.y + normal.z * LDIR.z);
        dot_light = max(0.0, dot_light);
        int intensity = (int)(250 * dot_light);
        intensity = max(0, min(255, intensity));

        for (int i = y1; i <= y2; i++) {
            for (int j = x1; j <= x2; j++) {
                array<double, 3> barycentric = get_barycentric_coordinates(j, i,
                    vertices[face.v1].x, vertices[face.v1].y,
                    vertices[face.v2].x, vertices[face.v2].y,
                    vertices[face.v3].x, vertices[face.v3].y);

                if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                    // Вычисляем z-координату через барицентрические координаты
                    double z_interpolated = barycentric[0] * vertices[face.v1].z + barycentric[1]
                                            * vertices[face.v2].z + barycentric[2] * vertices[face.v3].z;

                    int pixel_index = j + i * width;

                    // Z-буфер: рисуем только если пиксель ближе к камере
                    if (z_interpolated < zbuffer[pixel_index]) {
                        zbuffer[pixel_index] = z_interpolated; // Обновляем z-буфер

                        int img_index = pixel_index * 3;
                        image[img_index] = pixel_index + r * intensity;
                        image[img_index + 1] = pixel_index + g * intensity;
                        image[img_index + 2] = pixel_index + b * intensity;
                    }
                }
            }
        }
    }

    delete[] zbuffer;
}

void create_animation(const string& obj_filename, const string& output_gif = "animation.gif",
                      const int& frames = 36) {
    vector<Vertex> vertices;
    vector<Face> faces;
    read_obj(obj_filename, vertices, faces);
    double angle = M_PI;
    Matrix4 rot_matrix = rotationX(angle);
    Vector3D shift = {0, 0.1, 0};
    transform(vertices, shift, rot_matrix);

    GIF bnuuy_gif(200, 200);

    cout << "Creating " << frames << " frames..." << endl;
    system("mkdir frames");

    // Генерируем кадры
    for (int frame = 0; frame < frames; frame++) {
        vector<Vertex> temp_vertices = vertices;

        angle = - 2.0 * M_PI * frame / frames;
        rot_matrix = rotationY(angle);
        shift = {0, 0, 0};

        transform(temp_vertices, shift, rot_matrix);

        // Создаем RGB изображение
        unsigned char* rgb_image = new unsigned char[200 * 200 * 3];
        memset(rgb_image, 0, 200 * 200 * 3);

        draw_solid_obj_zbuffer(rgb_image, 200, 200, temp_vertices, faces, 0, 255, 255);

        // Конвертируем RGB в 8-bit grayscale для GIF
        unsigned char* gif_frame = new unsigned char[200 * 200];
        for (int i = 0; i < 200 * 200; i++) {
            // Простая конвертация RGB в grayscale
            int r = rgb_image[i * 3];
            int g = rgb_image[i * 3 + 1];
            int b = rgb_image[i * 3 + 2];
            gif_frame[i] = static_cast<unsigned char>((r + g + b) / 3);
        }

        bnuuy_gif.add_frame(gif_frame, 200, 200); // 200ms delay

        // Сохраняем PNG для отладки
        string frame_name = "frames/frame_" + to_string(frame) + ".png";
        lodepng::encode(frame_name, rgb_image, 200, 200, LCT_RGB, 8);

        delete[] rgb_image;

        // Прогресс
        if (frame % 1 == 0) {
            cout << "Progress: " << frame << "/" << frames << endl;
        }
    }

    if (bnuuy_gif.save(output_gif, false)) {
        cout << "Animation saved as: " << output_gif << endl;
    } else {
        cout << "Failed to save animation!" << endl;
    }
}

int main() {
    system("mkdir colors");
    system("mkdir shapes");
    system("mkdir model");
    // Генерация одноканальной картинки
    unsigned char* greyimage = new unsigned char[W * H];

    // Чёрная картинка
    for (int i = 0; i < H * W; i++) greyimage[i] = 0;
    lodepng::encode("colors/blackkartinka.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated black." << endl;

    // Белая картинка
    for (int i = 0; i < H * W; i++) greyimage[i] = 255;
    lodepng::encode("colors/whitekartinka.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated white." << endl;

    // Произвольная картинка
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            greyimage[x + W * y] = (x + y) / 2 < 256 ? (x + y) / 2 : 255;

    lodepng::encode("colors/gradkartinka.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated gradient." << endl;

    // Произвольная картинка #2
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            greyimage[x + x * y] = (x + y) % 256;

    lodepng::encode("color/kartinka?.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated something..." << endl;

    // Звезда 1
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star1(greyimage, W, H, 255, 95, 13);
    lodepng::encode("shapes/star1.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated star 1." << endl;

    // Звезда 2
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star2(greyimage, W, H, 255, 95, 13);
    lodepng::encode("shapes/star2.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated star 2." << endl;

    // Звезда 3
    for (int i = 0; i < W * H; i++) greyimage[i] = 0;
    star3(greyimage, W, H, 255, 95, 13);
    lodepng::encode("shapes/star3.png", greyimage, W, H, LCT_GREY, 8);
    cout << "Generated star 3." << endl;

    delete greyimage;

    // Генерация трёхканальной картинки
    unsigned char* rgbimage = new unsigned char[W * H * 3];

    // Красная картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 255; rgbimage[i + 1] = 0; rgbimage[i + 2] = 0;
    }
    lodepng::encode("colors/redkartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated red." << endl;

    // Зелёная картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 0; rgbimage[i + 1] = 255; rgbimage[i + 2] = 0;
    }
    lodepng::encode("colors/greenkartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated green." << endl;

    // Синяя картинка
    for (int i = 0; i < H * W * 3; i += 3) {
        rgbimage[i] = 0; rgbimage[i + 1] = 0; rgbimage[i + 2] = 255;
    }
    lodepng::encode("colors/bluekartinka.png", rgbimage, W, H, LCT_RGB, 8);
    cout << "Generated blue." << endl;

    delete rgbimage;

    // 3D-Модель
    unsigned char* modelimage = new unsigned char[1000 * 1000];
    vector<Vertex> vertices;
    vector<Face> faces;
    Matrix4 rot_matrix;
    Vector3D shift;
    string filename = "model_1.obj";

    read_obj(filename, vertices, faces);
    draw_wireframe_obj(modelimage, W, H, vertices, faces);
    lodepng::encode("model/wireframe.png", modelimage, 1000, 1000, LCT_GREY, 8);
    cout << "Generated wireframe model." << endl;

    for (int i = 0; i < 1000 * 1000; i++) modelimage[i] = 0;

    // Повёрнутый

    read_obj(filename, vertices, faces);
    rot_matrix = multiply_matrices(rotationY(-90*M_PI/180), rotationX(0*M_PI/180));
    rot_matrix = multiply_matrices(rot_matrix, rotationZ(180*M_PI/180));
    shift = {0,0.1,0};
    transform(vertices, shift, rot_matrix);

    draw_wireframe_obj(modelimage, 1000, 1000, vertices, faces);
    lodepng::encode("model/wireframe_rot.png", modelimage, 1000, 1000, LCT_GREY, 8);
    cout << "Generated wireframe model rotated." << endl;

    delete modelimage;

    // Разными цветами
    unsigned char* rgbmodelimage = new unsigned char[1000 * 1000 * 3];
    read_obj(filename, vertices, faces);
    // rot_matrix = rotationZ(180*M_PI/180);
    // shift = {0, -100, 0};
    // transform(vertices, shift, rot_matrix);
    draw_solid_obj_zbuffer(rgbmodelimage, 1000, 1000, vertices, faces);
    lodepng::encode("model/solid.png", rgbmodelimage, 1000, 1000, LCT_RGB, 8);
    cout << "Generated solid model." << endl;

    for (int i = 0; i < 1000 * 1000 * 3; i++) rgbmodelimage[i] = 0;

    // Повёрнутый
    read_obj(filename, vertices, faces);
    rot_matrix = multiply_matrices(rotationY(-90*M_PI/180), rotationX(0*M_PI/180));
    rot_matrix = multiply_matrices(rot_matrix, rotationZ(180*M_PI/180));
    shift = {0,0.1,0};
    transform(vertices, shift, rot_matrix);

    // draw_solid_obj_zbuffer(rgbmodelimage, 1000, 1000, vertices, faces, 255, 0, 0);
    draw_normals(rgbmodelimage, 1000, 1000, vertices, faces, 255, 0, 0);
    lodepng::encode("model/solid_rot1.png", rgbmodelimage, 1000, 1000, LCT_RGB, 8);
    cout << "Generated solid model rotated." << endl;

    for (int i = 0; i < 1000 * 1000 * 3; i++) rgbmodelimage[i] = 0;

    // Повёрнутый
    read_obj(filename, vertices, faces);
    rot_matrix = multiply_matrices(rotationY(-90*M_PI/180), rotationX(90*M_PI/180));
    rot_matrix = multiply_matrices(rot_matrix, rotationZ(180*M_PI/180));
    shift = {0,0.05,0};
    transform(vertices, shift, rot_matrix);

    draw_solid_obj_zbuffer(rgbmodelimage, 1000, 1000, vertices, faces);
    // debug_normals(rgbmodelimage, 1000, 1000, vertices, faces);
    lodepng::encode("model/solid_rot2.png", rgbmodelimage, 1000, 1000, LCT_RGB, 8);
    cout << "Generated solid model rotated." << endl;

    delete rgbmodelimage;

    // Отрисовка треугольников
    unsigned char* triangleimage = new unsigned char[1000 * 1000 * 3];
    vector<Vertex> ver = {{0.031709, 0.025774, -0.020613},
        {0.018666, 0.026160, -0.028820},
        {0.034755, 0.028198, -0.016804},
        {0.1, 0.2, 0.2},
        {0.1, -0.9, 0.01},
        {0.1, -0.4, 0.3},
        {0.2, -0.1, 0.2},
        {0.15, 0.01, -0.2},
        {0.1, 0.05, -0.11}};

    Face f1 = {0,1,2};
    Face f2 = {3,4,5};
    Face f3 = {6,7,8};
    draw_triangle(triangleimage, 1000, 1000, ver, {f1});
    draw_triangle(triangleimage, 1000, 1000, ver, {f2});
    draw_triangle(triangleimage, 1000, 1000, ver, {f3});

    lodepng::encode("shapes/triangle1.png", triangleimage, 1000, 1000, LCT_RGB, 8);
    cout << "Generated triangle." << std::endl;

    delete triangleimage;

    // GIF-анимация

    create_animation("model_1.obj");

    return 0;
}
