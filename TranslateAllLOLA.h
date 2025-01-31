
#include <string>
#include <cstdint> 

#define SCALING_FACTOR        0.5
#define OFFSET                1737400

typedef struct
{
    uint8_t header[80];
    uint32_t num_triangles;
}
stl_header;

typedef struct
{
    float x;
    float y;
    float z;
}
vertex;

typedef struct
{
    vertex normal;
    vertex v1;
    vertex v2;
    vertex v3;
}
triangle;

uint16_t attribute_byte_count;

typedef struct {
    int width;
    int height;
    uint8_t *data;
    size_t size;
}
ppm_image;

void translate_file(char* full_path, int latitude, int longitude);
void translate_file_to_stl(int latitude, int longitude);
void translate_file_to_uchar(char* full_path, int latitude, int longitude);
void translate_file_to_ushort(char* full_path, int latitude, int longitude);
void translate_file_to_bmp(char* full_path, int latitude, int longitude);
void translate_file_to_pgm(char* full_path, int latitude, int longitude);
float radius(float height);
double deg2rad (double degrees);