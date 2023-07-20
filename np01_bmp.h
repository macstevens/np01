/* np01_bmp.h 

Reference: 
    https://engineering.purdue.edu/ece264/17au/hw/HW15


Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#include <vector>
#include <stdint.h>

namespace np01 {

#pragma pack(push,2) 
struct np01_bmp_header{         // Total: 54 bytes
public:
    uint16_t  type;             // Magic identifier: 0x4d42
    uint32_t  size;             // File size in bytes
    uint16_t  reserved1;        // Not used
    uint16_t  reserved2;        // Not used
    uint32_t  offset;           // Offset to image data in bytes from beginning of file (54 bytes)
    uint32_t  dib_header_size;  // DIB Header size in bytes (40 bytes)
    int32_t   width_px;         // Width of the image
    int32_t   height_px;        // Height of image
    uint16_t  num_planes;       // Number of color planes
    uint16_t  bits_per_pixel;   // Bits per pixel
    uint32_t  compression;      // Compression type
    uint32_t  image_size_bytes; // Image size in bytes
    int32_t   x_resolution_ppm; // Pixels per meter
    int32_t   y_resolution_ppm; // Pixels per meter
    uint32_t  num_colors;       // Number of colors  
    uint32_t  important_colors; // Important colors 
};
#pragma pack(pop)

struct np01_bmp_file_init_params{
public:
    int32_t   width_px;         // Width of the image
    int32_t   height_px;        // Height of image
};

struct np01_bmp_color{
public:
    uint8_t m_red;
    uint8_t m_green;
    uint8_t m_blue;
    bool    m_blend;
public:
    np01_bmp_color(uint8_t r, uint8_t g, uint8_t b, bool blend=false):
        m_red(r),m_green(g),m_blue(b),m_blend(blend){}
};


class np01_bmp_file{
private:
    typedef std::vector<uint8_t> uint8_vec;
    typedef uint8_vec::iterator uint8_vec_itr;
    typedef uint8_vec::const_iterator uint8_vec_citr;
    typedef std::vector<uint8_vec> uint8_vec_vec;
    typedef uint8_vec_vec::iterator uint8_vec_vec_itr;
    typedef uint8_vec_vec::const_iterator uint8_vec_vec_citr;
private:
    np01_bmp_header m_header;
    uint8_vec_vec m_rows;
public:
    np01_bmp_file();
    np01_bmp_file(const np01_bmp_file_init_params& init_params);
   ~np01_bmp_file();
    void init(const np01_bmp_file_init_params& init_params);
    void draw_pixel(const int32_t& i,const int32_t& j, const np01_bmp_color& color);
    void draw_line(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1,const int32_t& j1, const np01_bmp_color& color);
    void draw_box(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1,const int32_t& j1, const np01_bmp_color& color);
    void draw_diamond(const int32_t& i_ctr,const int32_t& j_ctr, 
        const int32_t& w, const np01_bmp_color& color);
    void write_file(const char *file_name) const;
};


}
