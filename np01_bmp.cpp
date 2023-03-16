/* np01_bmp.cpp 

Copyright (c) 2021 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "np01_bmp.h"

namespace np01 {


np01_bmp_file::np01_bmp_file(){
}

np01_bmp_file::np01_bmp_file(const np01_bmp_file_init_params& init_params){
init(init_params);
}

np01_bmp_file::~np01_bmp_file(){
}

void np01_bmp_file::init(const np01_bmp_file_init_params& init_params){
/* initialize rows */
size_t row_idx;
size_t row_sz = static_cast<size_t>(init_params.width_px) * 3;
if( (row_sz % 4) != 0){
    row_sz += ( 4 - (row_sz % 4) );}
assert( (row_sz % 4) == 0);
m_rows.clear();
m_rows.resize(init_params.height_px);
for(row_idx=0;row_idx!=static_cast<size_t>(init_params.height_px);++row_idx){
    (m_rows.at(row_idx)).resize(row_sz, 0); }

/* initialize header */
memset(&m_header, 0, sizeof(m_header));
const uint32_t img_sz = (static_cast<uint32_t>(init_params.height_px) *
    static_cast<uint32_t>(row_sz));
m_header.type=0x4d42;               // Magic identifier: 0x4d42
m_header.size = 54 + img_sz;        // File size in bytes
m_header.reserved1=0;               // Not used
m_header.reserved2=0;               // Not used
m_header.offset=54;                 // Offset to image data in bytes from beginning of file (54 bytes)
m_header.dib_header_size=40;        // DIB Header size in bytes (40 bytes)
m_header.width_px=init_params.width_px;     // Width of the image
m_header.height_px=init_params.height_px;   // Height of image
m_header.num_planes=1;              // Number of color planes
m_header.bits_per_pixel=24;         // Bits per pixel
m_header.compression=0;             // Compression type
m_header.image_size_bytes=img_sz;   // Image size in bytes
m_header.x_resolution_ppm=0;        // Pixels per meter
m_header.y_resolution_ppm=0;        // Pixels per meter
m_header.num_colors=0;              // Number of colors  
m_header.important_colors=0;        // Important colors 
}

void np01_bmp_file::draw_pixel(const int32_t& i,const int32_t& j,
    const np01_bmp_color& color){
if((i>=0) && (i<m_header.width_px) && (j>=0) && (j<m_header.height_px)){
    uint8_vec *row = &(m_rows.at(j));
    uint8_vec_itr byte_itr = (row->begin()) + (3*i);
    if(color.m_blend){
        const uint16_t b0 = static_cast<uint16_t>(*byte_itr);
        const uint16_t g0 = static_cast<uint16_t>(*(byte_itr+1));
        const uint16_t r0 = static_cast<uint16_t>(*(byte_itr+2));
        const uint16_t b1 = b0 + static_cast<uint16_t>(color.m_blue);
        const uint16_t g1 = g0 + static_cast<uint16_t>(color.m_green);
        const uint16_t r1 = r0 + static_cast<uint16_t>(color.m_red);
        static const uint16_t max_byte_val = 255;
        const uint16_t b_extra = (b1 > max_byte_val) ? (b1 - max_byte_val) : 0;
        const uint16_t g_extra = (g1 > max_byte_val) ? (g1 - max_byte_val) : 0;
        const uint16_t r_extra = (r1 > max_byte_val) ? (r1 - max_byte_val) : 0;
        const uint16_t extra = (b_extra + g_extra + r_extra) / 16;
        const uint16_t b2 = b1 + extra;
        const uint16_t g2 = g1 + extra;
        const uint16_t r2 = r1 + extra;
        const uint16_t b = (b2 < max_byte_val) ? b2 : max_byte_val;
        const uint16_t g = (g2 < max_byte_val) ? g2 : max_byte_val;
        const uint16_t r = (r2 < max_byte_val) ? r2 : max_byte_val;
        *byte_itr = static_cast<uint8_t>(b);
        *(byte_itr+1) = static_cast<uint8_t>(g);
        *(byte_itr+2) = static_cast<uint8_t>(r);
        }
    else{
        *byte_itr = color.m_blue;
        *(byte_itr+1) = color.m_green;
        *(byte_itr+2) = color.m_red;
        }
    }
}

void np01_bmp_file::draw_line(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1, const int32_t& j1, const np01_bmp_color& color){
int32_t i, j;
const int32_t i_step = (i1 > i0) ?  1 : -1;
const int32_t j_step = (j1 > j0) ?  1 : -1;
if((i0 == i1) || (j0 == j1)){
    if(i0 != i1){
        for(i = i0; i != i1; i += i_step){
            draw_pixel(i, j0, color);
            }
        }
    else if(j0 != j1){
        for(j = j0; j != j1; j += j_step){
            draw_pixel(i0, j, color);
            }
        }
    draw_pixel(i1, j1, color);
    }
else {
    draw_pixel(i0, j0, color);
    draw_pixel(i1, j1, color);
    const int32_t di = i1 - i0;
    const int32_t dj = j1 - j0;
    if(labs(di) > labs(dj)){
        const double m = static_cast<double>(dj) / static_cast<double>(di);
        for(int32_t ii = i_step; ii != di; ii += i_step){
            double jj = static_cast<double>(ii) * m;
            j = j0 + lrint(jj);
            draw_pixel(i0+ii, j, color);
            }
        }
    else{
        const double m = static_cast<double>(di) / static_cast<double>(dj);
        for(int32_t jj = j_step; jj != dj; jj += j_step){
            double ii = static_cast<double>(jj) * m;
            i = i0 + lrint(ii);
            draw_pixel(i, j0+jj, color);
            }
        }
    }
}

void np01_bmp_file::draw_box(const int32_t& i0,const int32_t& j0, 
        const int32_t& i1, const int32_t& j1, const np01_bmp_color& color){
const int32_t i_min = (i0 < i1) ? i0 : i1;
const int32_t i_max = (i0 > i1) ? i0 : i1;
const int32_t j_min = (j0 < j1) ? j0 : j1;
const int32_t j_max = (j0 > j1) ? j0 : j1;
for(int32_t i = i_min; i <= i_max; ++i){
    for(int32_t j = j_min; j <= j_max; ++j){
        draw_pixel(i, j, color);
        }
    }
}

void np01_bmp_file::draw_diamond(const int32_t& i_ctr,const int32_t& j_ctr, 
        const int32_t& w, const np01_bmp_color& color){
assert(w>=0);
const int32_t hw = w/2;
for(int32_t i = -hw; i <= hw; ++i){
    const int32_t j_min = labs(i) - hw;
    const int32_t j_max = -j_min;
    for(int32_t j = j_min; j <= j_max; ++j){
        draw_pixel(i_ctr+i, j_ctr+j, color);
        }
    }
}


void np01_bmp_file::write_file(const char *file_name) const{
std::ofstream os(file_name, std::ios::out | std::ios::binary);

const uint8_t *hdr_byte_ptr = (const uint8_t *)&m_header;
const uint8_t * const hdr_byte_end_ptr = hdr_byte_ptr + 54;
for(; (hdr_byte_ptr != hdr_byte_end_ptr) && os.good(); ++hdr_byte_ptr){
    const uint8_t& byte = *hdr_byte_ptr;
    os << byte;
    }

uint8_vec_vec_citr row_itr = m_rows.begin();
for(; (row_itr != m_rows.end()) && os.good(); ++row_itr){
    const uint8_vec& row = *row_itr;
    uint8_vec_citr byte_itr = row.begin();
    for(; (byte_itr != row.end()) && os.good(); ++byte_itr){
        const uint8_t& byte = *byte_itr;
        os << byte;
        }
    }
os.close();
}


}
