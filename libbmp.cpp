#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <stdexcept>
#include "libbmp.h"

// BMP_HEADER
void bmp_header_init_df(bmp_header* header, const int width, const int height)
{
    header->bfSize = (sizeof(bmp_pixel) * width + BMP_GET_PADDING(width)) * abs(height);
    header->bfReserved = 0;
    header->bfOffBits = 54;
    header->biSize = 40;
    header->biWidth = width;
    header->biHeight = height;
    header->biPlanes = 1;
    header->biBitCount = 24;
    header->biCompression = 0;
    header->biSizeImage = 0;
    header->biXPelsPerMeter = 0;
    header->biYPelsPerMeter = 0;
    header->biClrUsed = 0;
    header->biClrImportant = 0;
}

bmp_error bmp_header_write(const bmp_header* header, FILE* img_file)
{
    if (header == nullptr)
    {
        return BMP_HEADER_NOT_INITIALIZED;
    }
    if (img_file == nullptr)
    {
        return BMP_FILE_NOT_OPENED;
    }

    const unsigned short magic = BMP_MAGIC;
    if (fwrite(&magic, sizeof(magic), 1, img_file) != 1)
    {
        return BMP_ERROR;
    }

    if (fwrite(header, sizeof(bmp_header), 1, img_file) != 1)
    {
        return BMP_ERROR;
    }
    return BMP_OK;
}

bmp_error bmp_header_read(bmp_header* header, FILE* img_file)
{
    if (img_file == nullptr)
    {
        return BMP_FILE_NOT_OPENED;
    }

    unsigned short magic;
    if (fread(&magic, sizeof(magic), 1, img_file) != 1 || magic != BMP_MAGIC)
    {
        return BMP_INVALID_FILE;
    }

    if (fread(header, sizeof(bmp_header), 1, img_file) != 1)
    {
        return BMP_ERROR;
    }

    return BMP_OK;
}

// BMP_PIXEL
void bmp_pixel_init(bmp_pixel* pxl, const unsigned char red,
                   const unsigned char green, const unsigned char blue)
{
    pxl->red = red;
    pxl->green = green;
    pxl->blue = blue;
}

// BMP_IMG
void bmp_img_alloc(bmp_img* img)
{
    const size_t h = abs(img->img_header.biHeight);
    const size_t w = img->img_header.biWidth;

    try
    {
        // Use vector of vectors for automatic memory management
        img->img_pixels_vec.resize(h);
        for (size_t y = 0; y < h; y++)
        {
            img->img_pixels_vec[y].resize(w);
        }
        
        // Update the raw pointer array to point to the vector data
        img->img_pixels = new bmp_pixel*[h];
        for (size_t y = 0; y < h; y++)
        {
            img->img_pixels[y] = img->img_pixels_vec[y].data();
        }
    }
    catch (const std::bad_alloc&)
    {
        throw std::runtime_error("Failed to allocate memory for BMP image");
    }
}

void bmp_img_init_df(bmp_img* img, const int width, const int height)
{
    bmp_header_init_df(&img->img_header, width, height);
    bmp_img_alloc(img);
}

void bmp_img_free(bmp_img* img)
{
    if (img->img_pixels != nullptr)
    {
        delete[] img->img_pixels;
        img->img_pixels = nullptr;
    }
    // Vectors will clean up automatically
    img->img_pixels_vec.clear();
}

bmp_error bmp_img_write(const bmp_img* img, const char* filename)
{
    if (img == nullptr)
    {
        return BMP_HEADER_NOT_INITIALIZED;
    }

    FILE* img_file = fopen(filename, "wb");
    if (img_file == nullptr)
    {
        return BMP_FILE_NOT_OPENED;
    }

    const bmp_error err = bmp_header_write(&img->img_header, img_file);
    if (err != BMP_OK)
    {
        fclose(img_file);
        return err;
    }

    const int h = abs(img->img_header.biHeight);
    const int offset = (img->img_header.biHeight > 0 ? 0 : h - 1);
    const size_t padding = BMP_GET_PADDING(img->img_header.biWidth);

    for (int y = 0; y < h; y++)
    {
        const int idx = (offset < 0 ? y : h - 1 - y);

        if (fwrite(img->img_pixels[idx], sizeof(bmp_pixel),
                  img->img_header.biWidth, img_file) != (size_t)img->img_header.biWidth)
        {
            fclose(img_file);
            return BMP_ERROR;
        }
        if (padding > 0)
        {
            const unsigned char padding_data[3] = {0, 0, 0};
            if (fwrite(padding_data, 1, padding, img_file) != padding)
            {
                fclose(img_file);
                return BMP_ERROR;
            }
        }
    }

    fclose(img_file);
    return BMP_OK;
}

bmp_error bmp_img_read(bmp_img* img, const char* filename)
{
    FILE* img_file = fopen(filename, "rb");
    if (img_file == nullptr)
    {
        return BMP_FILE_NOT_OPENED;
    }

    const bmp_error err = bmp_header_read(&img->img_header, img_file);
    if (err != BMP_OK)
    {
        fclose(img_file);
        return err;
    }

    bmp_img_alloc(img);

    const int h = abs(img->img_header.biHeight);
    const int offset = (img->img_header.biHeight > 0 ? 0 : h - 1);
    const size_t padding = BMP_GET_PADDING(img->img_header.biWidth);

    for (int y = 0; y < h; y++)
    {
        const int idx = (offset < 0 ? y : h - 1 - y);

        if (fread(img->img_pixels[idx], sizeof(bmp_pixel),
                 img->img_header.biWidth, img_file) != (size_t)img->img_header.biWidth)
        {
            fclose(img_file);
            return BMP_ERROR;
        }
        
        if (padding > 0)
        {
            unsigned char padding_data[3];
            if (fread(padding_data, 1, padding, img_file) != padding)
            {
                fclose(img_file);
                return BMP_ERROR;
            }
        }
    }

    fclose(img_file);
    return BMP_OK;
}
