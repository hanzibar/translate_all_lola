# Translate All LOLA

A C++ utility for processing Lunar Orbiter Laser Altimeter (LOLA) data into various file formats. This tool allows for the conversion of lunar topographic data based on specific latitude and longitude coordinates.

## Features

- Convert LOLA data to multiple output formats:
  - STL (3D model format)
  - BMP (bitmap image)
  - PGM (portable graymap)
  - Raw data formats (unsigned char, unsigned short)
- Geographic coordinate processing with latitude and longitude support
- Scaling and offset adjustments for proper data representation
- Handles vertex and triangle mesh generation for 3D models

## Building the Project

The project uses CMake as its build system. To build:

```bash
mkdir build
cd build
cmake ..
make
```

The minimum required CMake version is 2.6.

## Usage

The program provides several translation functions for different output formats:

- `translate_file`: Main translation function
- `translate_file_to_stl`: Converts to STL format
- `translate_file_to_bmp`: Converts to BMP format
- `translate_file_to_pgm`: Converts to PGM format
- `translate_file_to_uchar`: Converts to unsigned char format
- `translate_file_to_ushort`: Converts to unsigned short format

## Dependencies

- C++11 compatible compiler
- CMake (>= 2.6)
- Standard C++ libraries

## Installation

After building, the executable will be installed to the system's bin directory, and headers will be installed to the include directory using CMake's install targets.
