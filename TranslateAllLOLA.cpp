// A simple program that translates LARGE LOLA IMG to 1024 by 1024 sections

#include "TranslateAllLOLA.h"

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <chrono>
#include <memory>
#include <iostream>
#include <cmath>
#include <png.h>
#include "libbmp.h"
#include <sys/stat.h>
#include <errno.h>

#define LINES			15360
#define LINE_SAMPLES	30720	// LSB_INTEGER
#define SAMPLE_BITS		16
#define TILE_SIZE		1024

// Function to ensure required directories exist
void ensure_output_directories()
{
    const char* dirs[] = 
    {
        "./one_deg_grids",
        "./one_deg_grids_uchar",
        "./one_deg_grids_ushort",
        "./lola_to_pgm"
    };

    for (const char* dir : dirs)
    {
        struct stat st;
        if (stat(dir, &st) != 0)
        {
            // Directory doesn't exist, create it
            #ifdef _WIN32
                if (mkdir(dir) != 0)
            #else
                if (mkdir(dir, 0755) != 0)
            #endif
                {
                    if (errno != EEXIST)
                    { 
                        // Check if another process didn't create it meanwhile
                        std::cerr << "Failed to create directory: " << dir << std::endl;
                        throw std::runtime_error("Failed to create required directory");
                    }
                }
        }
        else if (!S_ISDIR(st.st_mode))
        {
            // Path exists but is not a directory
            std::cerr << "Path exists but is not a directory: " << dir << std::endl;
            throw std::runtime_error("Path exists but is not a directory");
        }
    }
}

int main(int argc, char* argv[])
{
    try
    {
        // Ensure all required directories exist
        ensure_output_directories();

        char file_path[512];
        char file_name[512];
        char full_path[2048];  // Increased buffer size
        
        printf("main\n");
        
        // sprintf(file_path, "/home/jdneff/Projects/godot_projects/download_db_files/db_files/");
     	sprintf(file_path, "/home/jdneff/.var/app/org.godotengine.Godot/data/godot/app_userdata/Moon-13/lola_data/");
        
        for (int latitude=-90; latitude<90; latitude+=15)
        {
            for (int longitude=0; longitude<360; longitude+=30)
            {
                char full_path[2048];

                if (latitude >= 0)
                {
                    sprintf(file_name, "ldem_1024_%02dn_%02dn_%03d_%03d.img", latitude, latitude+15, longitude, longitude+30);
                }
                else
                {
                    sprintf(file_name, "ldem_1024_%02ds_%02ds_%03d_%03d.img", -latitude, -1*(latitude+15), longitude, longitude+30);
                }
                
                if (snprintf(full_path, sizeof(full_path), "%s%s", file_path, file_name) >= sizeof(full_path))
                {
                    std::cerr << "Path too long\n";
                    return 1;
                }
                printf("%s\n", full_path);

                translate_file(full_path, latitude, longitude);
                // translate_file_to_ushort(full_path, latitude, longitude);
                // translate_file_to_pgm(full_path, latitude, longitude);
            }
        }

        // translate_file_to_stl(23, 48);
        // translate_file_to_stl(24, 48);
        // translate_file_to_stl(25, 48);

        // translate_file_to_stl(23, 49);
        // translate_file_to_stl(24, 49);
        // translate_file_to_stl(25, 49);

        // translate_file_to_stl(23, 50);
        // translate_file_to_stl(24, 50);
        // translate_file_to_stl(25, 50);

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

void translate_file(char* full_path, int latitude, int longitude)
{
    std::ifstream image_file(full_path, std::ios::binary);
    if (!image_file)
    {
        throw std::runtime_error("Could not open input file");
    }

    const int x_tiles = LINE_SAMPLES/TILE_SIZE;
    const int y_tiles = LINES/TILE_SIZE;
    const long samples = LINES*LINE_SAMPLES;

    // Use vector instead of malloc for automatic memory management
    std::vector<unsigned short> file_array(samples);
    
    // Read the entire file at once
    image_file.read(reinterpret_cast<char*>(file_array.data()), samples * sizeof(unsigned short));
    image_file.close();

    for (int y_tile = 0; y_tile < y_tiles; y_tile++)
    {
        for (int x_tile = 0; x_tile < x_tiles; x_tile++)
        {
            char file_name[512];
            sprintf(file_name, "./one_deg_grids/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
            
            // Use RAII for file handling
            std::ofstream tile_file(file_name, std::ios::binary);
            if (!tile_file)
            {
                throw std::runtime_error("Could not open output file");
            }

            for (int line = 0; line < TILE_SIZE; line++)
            {
                long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

                for (int sample = 0; sample < TILE_SIZE; sample++)
                {
                    long sample_index = line_index + sample;
                    float height = static_cast<float>(file_array[sample_index]);
                    printf("%f\n", height);
                    tile_file.write(reinterpret_cast<const char*>(&height), sizeof(float));
                }
            }

            tile_file.close();

            printf("TILE DONE --------------------------------------------------\n");

            std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(1));

            // Read back for verification
            std::ifstream verify_file(file_name, std::ios::binary);
            if (!verify_file)
            {
                throw std::runtime_error("Could not open file for verification");
            }

            for (int y = 0; y < TILE_SIZE; y++)
            {
                for (int x = 0; x < TILE_SIZE; x++)
                {
                    float height;
                    verify_file.read(reinterpret_cast<char*>(&height), sizeof(float));
                    printf("%f\n", height);
                }
            }
        }
    }
    // No need for free() - vector handles cleanup automatically
}

// 1024 pix/deg by 1024 pix/deg
// one row for each 0.0009765625 degrees of latitude
// values are relative to a radius of 1737.4 km
// Map Resolution	1024 <PIX/DEG>
// Map Scale	29.6126469 <METERS/PIXEL>

// UNIT                  = METER
// SCALING_FACTOR        = 0.5
// OFFSET                = 1737400.

// NOTE:
// Conversion from Digital Number to HEIGHT, i.e. elevation in meters, is:
// HEIGHT = (DN * SCALING_FACTOR).
// The conversion from Digital Number to PLANETARY_RADIUS in meters is:
// PLANETARY_RADIUS = (DN * SCALING_FACTOR) + OFFSET,
// where OFFSET is the radius of a reference sphere.
// The planetopotential TOPOGRAPHY is PLANETARY_RADIUS - GEOID_RADIUS,
// where GEOID_RADIUS is derived from a gravitational equipotential model.
// By convention, the average GEOID_RADIUS at the equator is OFFSET.

void translate_file_to_stl(int latitude, int longitude)
{
	char file_path[512];
	char tile_file_name[512];
	char stl_file_name[512];
	char full_path[2048];
	float* file_array;
	
	FILE* stl_file;
	FILE* tile_file;
	stl_header header;
	vertex t_vertex;
	triangle t_triangle;

	sprintf(tile_file_name, "%03d_%03d_Tile.img", longitude, latitude);
	sprintf(stl_file_name, "%03d_%03d_tile.stl", longitude, latitude);
	sprintf(file_path, "/home/jdneff/Projects/godot_projects/translate_all_lola/one_deg_grids/");
	
	if (snprintf(full_path, sizeof(full_path), "%s%s", file_path, tile_file_name) >= sizeof(full_path))
    {
        std::cerr << "Path too long\n";
        return;
    }
	printf("full path: %s\n", full_path);
	tile_file = fopen(full_path, "rb");
	if (tile_file == NULL) printf("tile_file == NULL\n");

	file_array = new float[1024*1024*sizeof(float)];
	fread(file_array, sizeof(float), 1024*1024, tile_file);
	fclose(tile_file);

	stl_file = fopen(stl_file_name, "wb");
	
	header.num_triangles = 2*1023*1023;
	fwrite(&header, sizeof(header), 1, stl_file);

	for (int sub_longi = 0; sub_longi < 1023; sub_longi ++)
	{
		for (int sub_lati = 0; sub_lati < 1023; sub_lati ++)
		{
			float lat_long_frag = (float)1.0/1024;

			unsigned long i00 = sub_lati*1024 + sub_longi;
			unsigned long i10 = sub_lati*1024 + sub_longi + 1;
			unsigned long i11 = (sub_lati + 1)*1024 + sub_longi + 1;
			unsigned long i01 = (sub_lati + 1)*1024 + sub_longi;

			float h00 = file_array[i00];
			float h10 = file_array[i10];
			float h11 = file_array[i11];
			float h01 = file_array[i01];

			float radius00 = radius(h00);
			float radius10 = radius(h10);
			float radius11 = radius(h11);
			float radius01 = radius(h01);

			float phy_deg = latitude + sub_lati*lat_long_frag;
			float phy1_deg = latitude + (sub_lati + 1)*lat_long_frag;
			
			float theta_deg = longitude + sub_longi*lat_long_frag;
			float theta1_deg = longitude + (sub_longi + 1)*lat_long_frag;

			float phy = deg2rad(phy_deg);
			float phy1 = deg2rad(phy1_deg);
			
			float theta = deg2rad(theta_deg);
			float theta1 = deg2rad(theta1_deg);

			float x00 = radius00*cos(theta)*sin(phy);
			float y00 = radius00*cos(phy);
			float z00 = radius00*sin(theta)*sin(phy);

			float x10 = radius10*cos(theta1)*sin(phy);
			float y10 = radius10*cos(phy);
			float z10 = radius10*sin(theta1)*sin(phy);

			float x11 = radius11*cos(theta1)*sin(phy1);
			float y11 = radius11*cos(phy1);
			float z11 = radius11*sin(theta1)*sin(phy1);

			float x01 = radius01*cos(theta)*sin(phy1);
			float y01 = radius01*cos(phy1);
			float z01 = radius01*sin(theta)*sin(phy1);

			t_triangle.normal = {0.000000, 1.000000, 0.000000};
			t_triangle.v1 = {x01, y01, z01};
			t_triangle.v2 = {x11, y11, z11};
			t_triangle.v3 = {x00, y00, z00};

			fwrite(&t_triangle, sizeof(t_triangle), 1, stl_file);
			fwrite(&attribute_byte_count, sizeof(attribute_byte_count), 1, stl_file);

			t_triangle.normal = {0.000000, 1.000000, 0.000000};
			t_triangle.v1 = {x11, y11, z11};
			t_triangle.v2 = {x10, y10, z10};
			t_triangle.v3 = {x00, y00, z00};

			fwrite(&t_triangle, sizeof(t_triangle), 1, stl_file);
			fwrite(&attribute_byte_count, sizeof(attribute_byte_count), 1, stl_file);
		}
	}

	fclose(stl_file);

	delete file_array;
}

void translate_file_to_uchar(char* full_path, int latitude, int longitude)
{
    std::ifstream image_file(full_path, std::ios::binary);
    if (!image_file)
    {
        throw std::runtime_error("Could not open input file");
    }

    const int x_tiles = LINE_SAMPLES/TILE_SIZE;
    const int y_tiles = LINES/TILE_SIZE;
    const long samples = LINES*LINE_SAMPLES;

    // Use vector instead of malloc
    std::vector<unsigned short> file_array(samples);

    image_file.read(reinterpret_cast<char*>(file_array.data()), samples * sizeof(unsigned short));
    image_file.close();

    for (int y_tile = 0; y_tile < y_tiles; y_tile++)
    {
        for (int x_tile = 0; x_tile < x_tiles; x_tile++)
        {
            char file_name[512];
            sprintf(file_name, "./one_deg_grids_uchar/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
            
            std::ofstream tile_file(file_name, std::ios::binary);
            if (!tile_file)
            {
                throw std::runtime_error("Could not open output file");
            }

            for (int line = 0; line < TILE_SIZE; line++)
            {
                long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

                for (int sample = 0; sample < TILE_SIZE; sample++)
                {
                    long sample_index = line_index + sample;
                    unsigned short height = file_array[sample_index];
                    unsigned char uc_height = static_cast<unsigned char>(
                        (static_cast<float>(height)/32767.0f)*255.0f
                    );
                    tile_file.write(reinterpret_cast<const char*>(&uc_height), sizeof(unsigned char));
                }
            }

            tile_file.close();
            printf("x%i tile of %i, x%i tile of %i done\n", x_tile, x_tiles, y_tile, y_tiles);
        }
    }
    // Vector will automatically clean up
}

void translate_file_to_ushort(char* full_path, int latitude, int longitude)
{
    std::ifstream image_file(full_path, std::ios::binary);
    if (!image_file)
    {
        throw std::runtime_error("Could not open input file");
    }

    const int x_tiles = LINE_SAMPLES/TILE_SIZE;
    const int y_tiles = LINES/TILE_SIZE;
    const long samples = LINES*LINE_SAMPLES;

    // Use vector instead of malloc
    std::vector<unsigned short> file_array(samples);

    image_file.read(reinterpret_cast<char*>(file_array.data()), samples * sizeof(unsigned short));
    image_file.close();

    for (int y_tile = 0; y_tile < y_tiles; y_tile++)
    {
        for (int x_tile = 0; x_tile < x_tiles; x_tile++)
        {
            char file_name[512];
            sprintf(file_name, "./one_deg_grids_ushort/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
            
            std::ofstream tile_file(file_name, std::ios::binary);
            if (!tile_file)
            {
                throw std::runtime_error("Could not open output file");
            }

            for (int line = 0; line < TILE_SIZE; line++)
            {
                long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

                for (int sample = 0; sample < TILE_SIZE; sample++)
                {
                    long sample_index = line_index + sample;
                    unsigned short height = file_array[sample_index];
                    tile_file.write(reinterpret_cast<const char*>(&height), sizeof(unsigned short));
                }
            }

            tile_file.close();
            printf("x%i tile of %i, x%i tile of %i done\n", x_tile, x_tiles, y_tile, y_tiles);
        }
    }
    // Vector will automatically clean up
}

void translate_file_to_bmp(char* full_path, int latitude, int longitude)
{
	/*
	png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_bytepp row_pointers;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        return ;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        png_destroy_write_struct(&png_ptr, NULL);
        return ;
    }

     if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return ;
    }

    png_set_IHDR(png_ptr, info_ptr,
                 size, size, // width and height
                 16, // bit depth
                 PNG_COLOR_TYPE_GRAY, // color type
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    // Initialize rows of PNG.
    row_pointers = (png_bytepp)png_malloc(png_ptr,
        size*png_sizeof(png_bytep));

    for (int i=0; i<size; i++)
       row_pointers[i]=NULL;

    for (int i=0; i<size; i++)
       row_pointers[i]=png_malloc(png_ptr, size*2);

    //set row data
    for (y = 0; y < size; ++y) {
        png_bytep row = row_pointers[y];
        for (x = 0; x < size; ++x) {
                short color = x+y;
                *row++ = (png_byte)(color & 0xFF);
                *row++ = (png_byte)(color >> 8);
        }
    }

    // Actually write the image data.
    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    //png_write_image(png_ptr, row_pointers);

    // Cleanup.
    for (y = 0; y < size; y++) {
        png_free(png_ptr, row_pointers[y]);
    }
    png_free(png_ptr, row_pointers);
    png_destroy_write_struct(&png_ptr, &info_ptr);
	*/
}

void translate_file_to_pgm(char* full_path, int latitude, int longitude)
{
    const int image_width = 1024*30;
    const int image_height = 1024*15;
    const size_t num_pixels = static_cast<size_t>(image_width) * image_height;

    // Use vector instead of malloc
    std::vector<short> file_data(num_pixels);

    std::ifstream image_file(full_path, std::ios::binary);
    if (!image_file)
    {
        throw std::runtime_error("Could not open input file");
    }

    image_file.read(reinterpret_cast<char*>(file_data.data()), num_pixels * sizeof(short));
    image_file.close();

    printf("bytes read = %ld\n", num_pixels * sizeof(short));

    // Convert to unsigned
    for (size_t i = 0; i < num_pixels; i++)
    {
        file_data[i] *= -1;
    }

    std::string file_name;
    if (latitude >= 0)
    {
        file_name = std::string("./lola_to_pgm/ldem_1024_") + 
                   std::to_string(latitude).substr(0, 2) + "n_" +
                   std::to_string(latitude + 15).substr(0, 2) + "n_" +
                   std::to_string(longitude).substr(0, 3) + "_" +
                   std::to_string(longitude + 30).substr(0, 3) + ".pgm";
    }
    else
    {
        file_name = std::string("./lola_to_pgm/ldem_1024_") + 
                   std::to_string(-latitude).substr(0, 2) + "s_" +
                   std::to_string(-1 * (latitude + 15)).substr(0, 2) + "s_" +
                   std::to_string(longitude).substr(0, 3) + "_" +
                   std::to_string(longitude + 30).substr(0, 3) + ".pgm";
    }

    printf("%s\n", file_name.c_str());
    std::ofstream pgm_file(file_name, std::ios::binary);
    if (!pgm_file)
    {
        throw std::runtime_error("Could not open output file");
    }

    // Write PGM header
    pgm_file << "P5\n" << image_width << " " << image_height << "\n65535\n";
    
    // Write image data
    pgm_file.write(reinterpret_cast<const char*>(file_data.data()), 
                   num_pixels * sizeof(short));

    pgm_file.close();
    // Vector will automatically clean up
}

float radius(float height)
{
	float radius = (height*SCALING_FACTOR + OFFSET)*1e-6; // units are in megameters
	return radius;
}

double deg2rad(double degrees)
{
    return degrees*4.0*atan(1.0)/180.0;
}