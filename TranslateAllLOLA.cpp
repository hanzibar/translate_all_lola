// A simple program that translates LARGE LOLA IMG to 1024 by 1024 sections

#include "TranslateAllLOLA.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <png.h>
#include "libbmp.h"
// #include <math.h>
// #include <boost/asio.hpp>
// #include <boost/asio/serial_port.hpp>

#define LINES			15360
#define LINE_SAMPLES	30720	// LSB_INTEGER
#define SAMPLE_BITS		16
#define TILE_SIZE		1024

int main(int argc, char* argv[])
{
	char file_path[1024];
	char file_name[1024];

	printf("main\n");
	
	// sprintf(file_path, "/home/jdneff/Projects/godot_projects/download_db_files/db_files/");
 	sprintf(file_path, "/home/jdneff/.var/app/org.godotengine.Godot/data/godot/app_userdata/Moon-13/lola_data/");
	
	for (int latitude=-90; latitude<90; latitude+=15)
	{
		for (int longitude=0; longitude<360; longitude+=30)
		{
			char full_path[1024];

			if (latitude >= 0)
			{
				sprintf(file_name, "ldem_1024_%02dn_%02dn_%03d_%03d.img", latitude, latitude+15, longitude, longitude+30);
			}
			else
			{
				sprintf(file_name, "ldem_1024_%02ds_%02ds_%03d_%03d.img", -latitude, -1*(latitude+15), longitude, longitude+30);
			}
			
			sprintf(full_path, "%s%s", file_path, file_name);
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

void translate_file(char* full_path, int latitude, int longitude)
{
	FILE* image_file;
	image_file = fopen(full_path,"rb");

	int x_tiles = LINE_SAMPLES/TILE_SIZE;
	int y_tiles = LINES/TILE_SIZE;
	long samples = LINES*LINE_SAMPLES;

	long file_array_size = samples*sizeof(float);
	unsigned short* file_array = (unsigned short*)malloc(file_array_size);

	fread(file_array, sizeof(unsigned short), samples, image_file);
	fclose(image_file);

	for (int y_tile=0; y_tile<y_tiles; y_tile++)
	{
		for (int x_tile=0; x_tile<x_tiles; x_tile++)
		{
			char file_name[512];
			sprintf(file_name, "./one_deg_grids/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
			FILE* tile_file = fopen(file_name, "wb");

			if (tile_file == NULL) 
			{
				printf("File not opened\n");
				exit(0);
			}

			for (int line=0; line<TILE_SIZE; line++)
			{
				long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

				for (int sample=0; sample<TILE_SIZE; sample++)
				{
					long sample_index = line_index + sample;
					float height = (float)file_array[sample_index];
					printf("%f\n", height);
					fwrite(&height, sizeof(float), 1, tile_file);
				}
			}

			fclose(tile_file);

			printf("TILE DONE --------------------------------------------------\n");

			std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(1));

			tile_file = fopen(file_name, "rb");

			for (int y=0; y<TILE_SIZE; y++)
				for (int x=0; x<TILE_SIZE; x++)
				{
					float height;
					fread(&height, sizeof(float), 1, tile_file);
					printf("%f\n", height);
				}
			
			fclose(tile_file);
		}
	}

	free(file_array);
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
	char file_path[1024];
	char tile_file_name[1024];
	char stl_file_name[1024];
	char full_path[1024];
	float* file_array;
	
	FILE* stl_file;
	FILE* tile_file;
	stl_header header;
	vertex t_vertex;
	triangle t_triangle;

	sprintf(tile_file_name, "%03d_%03d_Tile.img", longitude, latitude);
	sprintf(stl_file_name, "%03d_%03d_tile.stl", longitude, latitude);
	sprintf(file_path, "/home/jdneff/Projects/godot_projects/translate_all_lola/one_deg_grids/");
	
	sprintf(full_path, "%s%s", file_path, tile_file_name);
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
	FILE* image_file;
	image_file = fopen(full_path,"rb");

	int x_tiles = LINE_SAMPLES/TILE_SIZE;
	int y_tiles = LINES/TILE_SIZE;
	long samples = LINES*LINE_SAMPLES;

	long file_array_size = samples*sizeof(float);
	unsigned short* file_array = (unsigned short*)malloc(file_array_size);

	fread(file_array, sizeof(unsigned short), samples, image_file);
	fclose(image_file);

	for (int y_tile=0; y_tile<y_tiles; y_tile++)
	{
		for (int x_tile=0; x_tile<x_tiles; x_tile++)
		{
			char file_name[512];
			sprintf(file_name, "./one_deg_grids_uchar/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
			FILE* tile_file = fopen(file_name, "wb");

			if (tile_file == NULL) 
			{
				printf("File not opened\n");
				exit(0);
			}

			for (int line=0; line<TILE_SIZE; line++)
			{
				long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

				for (int sample=0; sample<TILE_SIZE; sample++)
				{
					long sample_index = line_index + sample;
					unsigned short height = file_array[sample_index];
					// unsigned char uc_height = height >> 8;
					unsigned char uc_height = ((float)height/32767)*255;
					// printf("%u\t%u\n", height, (unsigned char)uc_height);
					fwrite(&uc_height, sizeof(unsigned char), 1, tile_file);
				}
			}

			fclose(tile_file);

			printf("x%i tile of %i, x%i tile of %i done\n", x_tile, x_tiles, y_tile, y_tiles);
		}
	}

	free(file_array);
}

void translate_file_to_ushort(char* full_path, int latitude, int longitude)
{
	FILE* image_file;
	image_file = fopen(full_path,"rb");

	int x_tiles = LINE_SAMPLES/TILE_SIZE;
	int y_tiles = LINES/TILE_SIZE;
	long samples = LINES*LINE_SAMPLES;

	long file_array_size = samples*sizeof(float);
	unsigned short* file_array = (unsigned short*)malloc(file_array_size);

	fread(file_array, sizeof(unsigned short), samples, image_file);
	fclose(image_file);

	for (int y_tile=0; y_tile<y_tiles; y_tile++)
	{
		for (int x_tile=0; x_tile<x_tiles; x_tile++)
		{
			char file_name[512];
			sprintf(file_name, "./one_deg_grids_ushort/%03d_%03d_Tile.img", longitude+x_tile, (latitude+90)+y_tile);
			FILE* tile_file = fopen(file_name, "wb");

			if (tile_file == NULL) 
			{
				printf("File not opened\n");
				exit(0);
			}

			for (int line=0; line<TILE_SIZE; line++)
			{
				long line_index = (y_tile*TILE_SIZE + line)*LINE_SAMPLES + x_tile*TILE_SIZE;

				for (int sample=0; sample<TILE_SIZE; sample++)
				{
					long sample_index = line_index + sample;
					unsigned short height = file_array[sample_index];
					fwrite(&height, sizeof(unsigned short), 1, tile_file);
				}
			}

			fclose(tile_file);

			printf("x%i tile of %i, x%i tile of %i done\n", x_tile, x_tiles, y_tile, y_tiles);
		}
	}

	free(file_array);
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
	FILE* image_file;
	int image_width = 1024*30;
	int image_height = 1024*15;
	size_t block_size = image_width*image_height*sizeof(short);
	printf("block_size = %ld\n", block_size);

	short* file_data = (short*)malloc(block_size);
	if (file_data == NULL) printf("file_data == NULL\n");

	image_file = fopen(full_path,"rb");
	if (image_file == NULL) printf("image_file == NULL\n");

	long num_read = fread(file_data, sizeof(short), image_width*image_height, image_file);
	fclose(image_file);
	printf("bytes read = %ld\n", num_read*sizeof(short));

	//convert to unsigned
	for (long i=0; i<image_width*image_height; i++)
	{
		file_data[i] *= -1;
	}

	char file_name[1024];

	if (latitude >= 0)
	{
		sprintf(file_name, "./lola_to_pgm/ldem_1024_%02dn_%02dn_%03d_%03d.pgm", latitude, latitude+15, longitude, longitude+30);
	}
	else
	{
		sprintf(file_name, "./lola_to_pgm/ldem_1024_%02ds_%02ds_%03d_%03d.pgm", -latitude, -1*(latitude+15), longitude, longitude+30);
	}
	
	printf("%s\n", file_name);
	FILE* bpm_file = fopen(file_name, "wb");

	fprintf(bpm_file, "P5\n");
	fprintf(bpm_file, "%i %i\n", image_width, image_height);
	fprintf(bpm_file, "%i\n", 65535);
	fwrite(file_data, sizeof(short), image_width*image_height, bpm_file);

	fclose(bpm_file);
	free(file_data);
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