/* np01.h Newport Algorithms version 1

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/
#ifndef NP01_H
#define NP01_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <string>
#include <set>
#include <utility>
#include <vector>

namespace np01 {

typedef uint8_t np01_uint8;
typedef uint16_t np01_uint16;
typedef uint32_t np01_uint32;
typedef uint64_t np01_uint64;
typedef double np01_float64;

class np01_xy_node;
class np01_xy_node_map;
class np01_xy_edge;
class np01_xy_edge_loc_grid_node;
class np01_xy_edge_map;
class np01_xy_pair_map;

typedef std::pair<np01_uint16,np01_uint16> np01_uint16_pair;

typedef std::vector<np01_uint32> np01_uint32_vec;
typedef np01_uint32_vec::const_iterator np01_uint32_vec_citr;
typedef np01_uint32_vec::iterator np01_uint32_vec_itr;

typedef std::pair<np01_uint32,np01_uint32> np01_uint32_pair;
typedef std::vector<np01_uint32_pair> np01_uint32_pair_vec;
typedef np01_uint32_pair_vec::const_iterator np01_uint32_pair_vec_citr;
typedef np01_uint32_pair_vec::iterator np01_uint32_pair_vec_itr;

typedef std::vector<np01_float64> np01_float64_vec;
typedef np01_float64_vec::const_iterator np01_float64_vec_citr;
typedef np01_float64_vec::iterator np01_float64_vec_itr;

typedef std::pair<np01_float64,np01_float64> np01_xy;
typedef std::vector<np01_xy> np01_xy_vec;
typedef np01_xy_vec::const_iterator np01_xy_vec_citr;
typedef np01_xy_vec::iterator np01_xy_vec_itr;


enum np01_xy_group_type{ NP01_XY_GROUP_TYPE_UNKNOWN, NP01_XY_GROUP_TYPE_A,
    NP01_XY_GROUP_TYPE_B };

/* (x,y) point 

   +-----------------+    +-----------------+    +-----------------+
   |np01_xy_edge_map |    |np01_xy_edge_map |    |np01_xy_edge_map |
   |m_xy_group_type=A|    |                 |    |m_xy_group_type=B|  
   +-----------------+    +-----------------+    +-----------------+
                 ^                      ^                      ^
                 |                      |                      |
   +-------------|---+                  |        +-------------|---+
   |np01_xy_node |   |    +-------------|---+    |np01_xy_node |   |
   |m_owner------+   |    |np01_xy_edge |   |    |m_owner------+   |
   |m_edge--------------->|m_owner------+   |<----m_edge           |
   |m_xy_group_type=A|<----m_node_a         |    |m_xy_group_type=B|
   |m_x, m_y         |    |m_node_b------------->|m_x, m_y         |   
   +-----------------+    +-----------------+    +-----------------+

*/
class np01_xy_node{
private:
    np01_uint32 m_idx;
    np01_xy_node_map *m_owner;
    np01_xy_edge *m_edge;
    np01_uint16 m_loc_grid_i, m_loc_grid_j;
    np01_xy_node *m_loc_grid_prev, *m_loc_grid_next;
    np01_float64 m_x, m_y;
    np01_xy_group_type m_xy_group_type;
public:
    np01_xy_node();
   ~np01_xy_node();

    void set_idx(const np01_uint32& i){ m_idx=i; }
    void set_owner(np01_xy_node_map *owner){ m_owner=owner; }
    void set_edge(np01_xy_edge *edge){ m_edge=edge; }
    void set_loc_grid_i(const np01_uint16& i){ m_loc_grid_i=i; }
    void set_loc_grid_j(const np01_uint16& j){ m_loc_grid_j=j; }
    void set_loc_grid_prev(np01_xy_node *p){ m_loc_grid_prev = p; }
    void set_loc_grid_next(np01_xy_node *n){ m_loc_grid_next = n; }
    void set_x(const np01_float64& x){ m_x=x; }
    void set_y(const np01_float64& y){ m_y=y; }
    void set_xy(const np01_float64& x, const np01_float64& y){ m_x=x; m_y=y; }
    void set_xy_group_type( const np01_xy_group_type& t ){ m_xy_group_type = t; }
    void disconnect_edge();

    np01_uint32 get_idx() const { return m_idx; }
    np01_xy_node_map *get_owner() const { return m_owner; }
    np01_xy_edge *get_edge() const { return m_edge; }
    np01_uint16 get_loc_grid_i() const { return m_loc_grid_i; }
    np01_uint16 get_loc_grid_j() const { return m_loc_grid_j; }
    np01_xy_node *get_loc_grid_prev() const { return m_loc_grid_prev; }
    np01_xy_node *get_loc_grid_next() const { return m_loc_grid_next; }
    const np01_float64& get_x() const {return m_x; }
    const np01_float64& get_y() const {return m_y; }
    np01_xy_group_type get_xy_group_type() const { return m_xy_group_type; }

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;

    static np01_uint32 get_invalid_idx(){
        return std::numeric_limits<np01_uint32>::max(); }
};


/* initialize from points */
struct np01_xy_loc_grid_dim_init_params{
    size_t point_count; /* number of points */
    np01_xy bb_min_xy, bb_max_xy;  /* bounding box of points */
    np01_float64 loc_grid_density; /* target density = points / grid square */
    size_t max_loc_grid_sq_count; /* target max number of grid squares */
};


/* locator grid dimensions 

           < - - - - - - - - - w - - - - - - - - - >
                      (number of squares)
           +-------+-------+-------+-------+-------+-   ^
           |       |       |       |       |       |
           |       |       |       |       |       |    |
           |       |       |       |       |       |
           +-------+-------+-------+-------+-------+-   |
           |       |       |       |       |       |     
           |       |       |       |       |       |    h (number of squares)
           |       |       |       |       |       |
       ^   +-------+-------+-------+-------+-------+-   |
       |   |       |       |       |       |       |
     sq_sz |       |       |       |       |       |    |
       |   |       |       |       |       |       |
       v   +-------+-------+-------+-------+-------+-   v
           (x_min, y_min)                  <-sq_sz->
*/
class np01_xy_loc_grid_dim {
private:
    np01_uint16 m_w, m_h; /* width, height */
    np01_float64 m_x_min, m_y_min; /* lower left corner */
    np01_float64 m_sq_size; /* square size */
public:
    np01_xy_loc_grid_dim(): m_w(0), m_h(0), m_x_min(0.0), m_y_min(0.0),
        m_sq_size(0.0){}
    np01_xy_loc_grid_dim(const np01_xy_loc_grid_dim& d): m_w(d.m_w), 
        m_h(d.m_h), m_x_min(d.m_x_min), m_y_min(d.m_y_min),
        m_sq_size(d.m_sq_size){}
   ~np01_xy_loc_grid_dim() {}

    np01_xy_loc_grid_dim& operator=(const np01_xy_loc_grid_dim& d)
        { m_w=d.m_w; m_h=d.m_h; m_x_min=d.m_x_min; m_y_min=d.m_y_min;
        m_sq_size=d.m_sq_size; return *this; }
    bool operator==(const np01_xy_loc_grid_dim& other) const {
        return ((m_w == other.m_w) && (m_h == other.m_h) &&
            (m_x_min == other.m_x_min) && (m_y_min == other.m_y_min) &&
            (m_sq_size == other.m_sq_size)); }
    bool operator!=(const np01_xy_loc_grid_dim& other) const {
        return ((m_w != other.m_w) && (m_h != other.m_h) &&
            (m_x_min != other.m_x_min) && (m_y_min != other.m_y_min) &&
            (m_sq_size != other.m_sq_size)); }

    void reset(){ m_w=0; m_h=0; m_x_min=0.0; m_y_min=0.0; m_sq_size=0.0; }
    void init(const np01_xy_loc_grid_dim_init_params *init_params);
    void set_w(const np01_uint16& w){ m_w = w; }
    void set_h(const np01_uint16& h){ m_h = h; }
    void set_x_min(const np01_float64& x_min){ m_x_min = x_min; }
    void set_y_min(const np01_float64& y_min){ m_y_min = y_min; }
    void set_sq_size(const np01_float64& sq_size){ m_sq_size = sq_size; }

    np01_uint16 get_w() const { return m_w; }
    np01_uint16 get_h() const { return m_h; }
    const np01_float64& get_x_min() const { return m_x_min; }
    const np01_float64& get_y_min() const { return m_y_min; }
    const np01_float64& get_sq_size() const { return m_sq_size; }

    np01_uint16 get_i(const np01_float64& x) const {
        np01_uint16 i = 0;
        if((m_sq_size>0.0) && (m_w > 0)){
            np01_float64 dx = x - m_x_min;  np01_float64 ii = dx / m_sq_size;
            i= (ii<1.0)? 0 : (ii>=(np01_float64)m_w)? m_w-1: (np01_uint16)ii; }
        return i; }
    np01_uint16 get_j(const np01_float64& y) const {
        np01_uint16 j = 0;
        if((m_sq_size>0.0) && (m_h > 0)){
            np01_float64 dy = y - m_y_min;  np01_float64 jj = dy / m_sq_size;
            j= (jj<1.0)? 0 : (jj>=(np01_float64)m_h)? m_h-1: (np01_uint16)jj; }
        return j; }
    void get_bb_indices( const np01_xy& xy_min, const np01_xy& xy_max,
        np01_uint16_pair *ij_min, np01_uint16_pair *ij_max) const;
    np01_float64 get_sq_ctr_x( const np01_uint16& i ) const{
        np01_float64 ctr_x;
        if(m_w>0){
            np01_float64 ii=static_cast<np01_float64>((i>=m_w)?(m_w-1):i);
            ctr_x = m_x_min + ( m_sq_size * ( ii + 0.5 ) ); }
        else{ ctr_x = m_x_min; }
        return ctr_x; }
    np01_float64 get_sq_ctr_y( const np01_uint16& j ) const{
        np01_float64 ctr_y;
        if(m_h>0){
            np01_float64 jj=static_cast<np01_float64>((j>=m_h)?(m_h-1):j);
            ctr_y = m_y_min + ( m_sq_size * ( jj + 0.5 ) ); }
        else{ ctr_y = m_y_min; }
        return ctr_y; }

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
};

/* container for set of (x,y) points */
class np01_xy_node_map {
public:
    typedef std::vector<np01_xy_node *> node_vec;
    typedef node_vec::iterator node_vec_itr;
    typedef node_vec::const_iterator node_vec_citr;
private:
    np01_xy_pair_map *m_owner;
    np01_xy_loc_grid_dim m_loc_grid_dim;
    node_vec m_loc_grid;
    node_vec m_node_vec;
    np01_xy_group_type m_xy_group_type;
public:
    np01_xy_node_map();
   ~np01_xy_node_map();

    void set_owner( np01_xy_pair_map *owner ){ m_owner = owner; }
    void init_loc_grid( np01_xy_loc_grid_dim& d );
    void reserve( const np01_uint32& sz ) { m_node_vec.reserve(sz); }
    np01_xy_node *create_node( const np01_xy& xy );
    void set_xy_group_type( const np01_xy_group_type& t ){ m_xy_group_type = t;}

    np01_xy_pair_map *get_owner() const { return m_owner; }
    const np01_xy_loc_grid_dim& get_loc_grid_dim()const{return m_loc_grid_dim;}
    void get_nodes_in_bb(const np01_xy& xy_min, const np01_xy& xy_max,
        node_vec *nodes) const;
    void get_nodes_in_circle(const np01_xy& xy_ctr, const np01_float64& r,
        node_vec *nodes) const;
    void get_all_nodes(node_vec *nodes) const;
    np01_xy_node *get_near_node(const np01_xy& xy, np01_float64 *d,
        node_vec *utility_vec) const;
    np01_xy_node *get_near_unconnected_node(const np01_xy& xy,
        const node_vec *usable_node_lookup_table, node_vec *utility_vec)const;
    np01_uint32 get_node_count() const {
        return static_cast<np01_uint32>(m_node_vec.size()); }
    np01_xy_node *get_node_by_idx( const np01_uint32& i ) const {
        return (( i < m_node_vec.size() ) ? m_node_vec[i] : NULL); }
    np01_xy_group_type get_xy_group_type() const { return m_xy_group_type; }

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
private:
    void insert_node_in_loc_grid(np01_xy_node *n);
    void clear_loc_grid();
    void free_nodes();
};


/* edge from (xa,ya) to (xb,yb) */
class np01_xy_edge{
private:
    np01_uint32 m_idx;
    np01_xy_edge_map *m_owner;
    np01_xy_node *m_node_a, *m_node_b;
    np01_xy_edge_loc_grid_node *m_loc_grid_head_node;
    np01_float64 m_ab_fwd_x, m_ab_fwd_y; /* unit vector A->B */
    np01_float64 m_fwd_dot_a, m_fwd_dot_b; /* unit vector dot (xa,ya) and (xb, yb) */
    np01_float64 m_fwd_cross_ab; /* unit vector cross ((xa+xb)/2,(ya+yb)/2) */

public:
    np01_xy_edge();
   ~np01_xy_edge();

    void set_idx( const np01_uint32& i) { m_idx = i; }
    void set_owner( np01_xy_edge_map *owner ) { m_owner = owner; }
    void set_node_a( np01_xy_node *a ) { m_node_a = a; }
    void set_node_b( np01_xy_node *b ) { m_node_b = b; }
    void set_loc_grid_head_node( np01_xy_edge_loc_grid_node *n ){
        m_loc_grid_head_node = n; }
    void set_ab_fwd_x( const np01_float64& x ){ m_ab_fwd_x = x; }
    void set_ab_fwd_y( const np01_float64& y ){ m_ab_fwd_y = y; }
    void set_fwd_dot_a( const np01_float64& da ){ m_fwd_dot_a = da; }
    void set_fwd_dot_b( const np01_float64& db ){ m_fwd_dot_b = db; }
    void set_fwd_cross_ab( const np01_float64& cab ){ m_fwd_cross_ab = cab; }
    void update_fwd_dot_cross_ab();

    np01_uint32 get_idx() const { return m_idx; }
    np01_xy_edge_map *get_owner() const { return m_owner; };
    np01_xy_node *get_node_a() const { return m_node_a; }
    np01_xy_node *get_node_b() const { return m_node_b; }
    np01_xy_edge_loc_grid_node *get_loc_grid_head_node() const {
        return m_loc_grid_head_node; }
    const np01_float64& get_ab_fwd_x() const { return m_ab_fwd_x; }
    const np01_float64& get_ab_fwd_y() const { return m_ab_fwd_y; }
    const np01_float64& get_fwd_dot_a() const { return m_fwd_dot_a; }
    const np01_float64& get_fwd_dot_b() const { return m_fwd_dot_b; }
    const np01_float64& get_fwd_cross_ab() const { return m_fwd_cross_ab; }
    bool is_in_bb(const np01_xy& xy_min, const np01_xy& xy_max) const;
    bool is_in_circle(const np01_xy& xy_ctr, const np01_float64& r) const;
    bool intersects(const np01_xy_edge* other) const;

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
};

/*
  +---------------------------------------------------------------------------------------+
  |np01_xy_edge                                                                           |
  |m_loc_grid_head_node+                                                                  |
  +--------------------|------------------------------------------------------------------+
                       |  ^                   ^                   ^                   ^
                       v  |                   |                   |                   |
                  +-------|-----+     +-------|-----+     +-------|-----+     +-------|-----+
                  |m_owner+     |     |m_owner+     |     |m_owner+     |     |m_owner+     |
                  |m_e_prev=NULL|<-----m_e_prev     |<-----m_e_prev     |<-----m_e_prev     |
                  |m_e_next---------->|m_e_next---------->|m_e_next---------->|m_e_next=NULL|
                  +-------------+     +-------------+     +-------------+     +-------------+




  +----------------------------------------------------------------------------------------+
  |np01_xy_edge_map                                                                        |
  |            +-----------------+-----------------+-----------------+-----------------+   |                                                   |
  |            |   i=0     j=0   |   i=0     j=1   |   i=1     j=0   |   i=1     j=1   |   |
  |            +-----------------+-----------------+-----------------+-----------------+   |                                                   |
  | m_loc_grid | 0    *          | 1    *          | 2    *          | 3    *          |   |
  |            +------|----------+------|----------+------|----------+------|----------+   |                                                                         |
  +-------------------|-----------------|-----------------|-----------------|--------------+
                      |                 |                 |                 |
                      v                 v                 v                 v
           +------------+    +------------+    +------------+    +------------+
           |i=0 j=0     |    |i=0 j=1     |    |i=1 j=0     |    |i=1 j=1     |  
           |m_owner=1   |    |m_owner=1   |    |m_owner=2   |    |m_owner=1   |  
           |m_prev=NULL |    |m_prev=NULL |    |m_prev=NULL |    |m_prev=NULL |
           |m_next-+    |    |m_next-+    |    |m_next-+    |    |m_next-+    |
           +-------|----+    +-------|----+    +-------|----+    +-------|----+
                   |  ^              |  ^              |  ^              |  ^
                   v  |              v  |              v  |              v  |
           +----------|-+    +----------|-+    +----------|-+    +----------|-+
           |i=0 j=0   | |    |i=0 j=1   | |    |i=1 j=0   | |    |i=1 j=1   | |  
           |m_owner=2 | |    |m_owner=2 | |    |m_owner=3 | |    |m_owner=3 | |  
           |m_prev----+ |    |m_prev----+ |    |m_prev----+ |    |m_prev----+ |
           |m_next=NULL |    |m_next-+    |    |m_next-+    |    |m_next-+    |
           +------------+    +-------|----+    +-------|----+    +-------|----+
                                     |  ^              |  ^              |  ^  
                                     v  |              v  |              v  |  
                             +----------|-+    +----------|-+    +----------|-+
                             |i=0 j=1   | |    |i=1 j=0   | |    |i=1 j=1   | |
                             |m_owner=3 | |    |m_owner=4 | |    |m_owner=5 | |
                             |m_prev----+ |    |m_prev----+ |    |m_prev----+ |
                             |m_next-+    |    |m_next=NULL |    |m_next=NULL |
                             +-------|----+    +------------+    +------------+
                                     |  ^     
                                     v  |     
                             +----------|-+   
                             |i=0 j=1   | |   
                             |m_owner=4 | |   
                             |m_prev----+ |   
                             |m_next=NULL |   
                             +------------+   
                 
*/
class np01_xy_edge_loc_grid_node{
private:
    np01_xy_edge *m_owner;
    np01_xy_edge_map *m_edge_map;
    np01_uint16 m_i, m_j;
    np01_xy_edge_loc_grid_node *m_prev, *m_next; /* same grid square */
    np01_xy_edge_loc_grid_node *m_e_prev, *m_e_next; /* same edge */
public:
    np01_xy_edge_loc_grid_node(): m_owner(NULL), m_edge_map(NULL), m_i(0), m_j(0),
        m_prev(NULL), m_next(NULL), m_e_prev(NULL), m_e_next(NULL) {}
   ~np01_xy_edge_loc_grid_node(){}

    void set_owner( np01_xy_edge *owner ){ m_owner = owner; }
    void set_edge_map( np01_xy_edge_map *m ){m_edge_map = m; }
    void set_i( const np01_uint16& i){ m_i = i; }
    void set_j( const np01_uint16& j){ m_j = j; }
    void set_prev( np01_xy_edge_loc_grid_node *p){ m_prev = p; }
    void set_next( np01_xy_edge_loc_grid_node *n){ m_next = n; }
    void set_e_prev( np01_xy_edge_loc_grid_node *p){ m_e_prev = p; }
    void set_e_next( np01_xy_edge_loc_grid_node *n){ m_e_next = n; }

    np01_xy_edge *get_owner() const { return m_owner; }
    np01_xy_edge_map *get_edge_map() const { return m_edge_map; }
    np01_uint16 get_i() const { return m_i; }
    np01_uint16 get_j() const { return m_j; }
    np01_xy_edge_loc_grid_node *get_prev() const { return m_prev; }
    np01_xy_edge_loc_grid_node *get_next() const { return m_next; }
    np01_xy_edge_loc_grid_node *get_e_prev() const { return m_e_prev; }
    np01_xy_edge_loc_grid_node *get_e_next() const { return m_e_next; }

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
};


class np01_xy_edge_map{
public:
    typedef std::vector<np01_xy_edge *> edge_vec;
    typedef edge_vec::iterator edge_vec_itr;
    typedef edge_vec::const_iterator edge_vec_citr;
    typedef std::vector<np01_xy_edge_loc_grid_node *> loc_grid_node_vec;
    typedef loc_grid_node_vec::iterator loc_grid_node_vec_itr;
    typedef loc_grid_node_vec::const_iterator loc_grid_node_vec_citr;
    typedef std::pair<np01_uint32, np01_xy_edge *> idx_edge_ptr_pair;
    typedef std::vector<idx_edge_ptr_pair> idx_edge_ptr_pair_vec;
    typedef idx_edge_ptr_pair_vec::const_iterator idx_edge_ptr_pair_vec_citr;
    typedef idx_edge_ptr_pair_vec::iterator idx_edge_ptr_pair_vec_itr;
private:
    np01_xy_pair_map *m_owner;
    edge_vec m_edge_vec;
    np01_xy_loc_grid_dim m_loc_grid_dim;
    loc_grid_node_vec m_loc_grid;
    np01_xy_edge_loc_grid_node *m_loc_grid_node_free_chain;
public:
    np01_xy_edge_map();
   ~np01_xy_edge_map();

    void set_owner( np01_xy_pair_map *owner ){m_owner = owner; }
    void reserve( const np01_uint32& sz ) { m_edge_vec.reserve(sz); }
    void init_loc_grid( np01_xy_loc_grid_dim& d );
    np01_xy_edge *create_edge( np01_xy_node *node_a, np01_xy_node *node_b );
    void update_edge_ab(np01_xy_edge *edge, np01_xy_node *node_a,
        np01_xy_node *node_b);
    void insert_edge_in_loc_grid(np01_xy_edge *edge);
    void remove_edge_from_loc_grid(np01_xy_edge *edge);

    np01_xy_pair_map *get_owner() const { return m_owner; }
    const np01_xy_loc_grid_dim& get_loc_grid_dim()const{return m_loc_grid_dim;}
    const np01_xy_edge_loc_grid_node *get_loc_grid_head_node(
        const np01_uint16 i, const np01_uint16 j){
        const size_t loc_grid_idx = static_cast<size_t>(j) +
           (static_cast<size_t>(i)*static_cast<size_t>(m_loc_grid_dim.get_h()));
        return (loc_grid_idx < m_loc_grid.size()) ?
            m_loc_grid[loc_grid_idx] : NULL; }
    void get_edges_near_edge(const np01_xy_edge *e, edge_vec *edges,
        idx_edge_ptr_pair_vec *utility_vec ) const;
    void get_edges_in_bb(const np01_xy& xy_min, const np01_xy& xy_max,
        edge_vec *edges, idx_edge_ptr_pair_vec *utility_vec ) const;
    void get_edges_in_circle(const np01_xy& xy_ctr, const np01_float64& r,
        edge_vec *edges, idx_edge_ptr_pair_vec *utility_vec) const;
    void get_intersecting_edges(const np01_xy_edge *e, edge_vec *edges,
        idx_edge_ptr_pair_vec *utility_vec ) const;
    np01_uint32 get_edge_count() const {
        return static_cast<np01_uint32>(m_edge_vec.size()); }
    np01_xy_edge *get_edge_by_idx( const np01_uint32& i ) const {
        return (( i < m_edge_vec.size() ) ? m_edge_vec[i] : NULL); }
    void get_edge_crossover_count( np01_uint32_vec *edge_crossover_count_vec,
        np01_uint32 *total_edge_crossover_count) const;

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    int verify_edges_fully_connected( char *err_msg,
        const size_t err_msg_capacity, size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    std::ostream& ostream_output_brief_table(std::ostream& os) const;
private:
    np01_xy_edge_loc_grid_node *alloc_loc_grid_node();
    void free_loc_grid_node(np01_xy_edge_loc_grid_node *n);
    void clear_loc_grid();
    void free_all_edges();
    void free_free_chain();
};


struct np01_xy_pair_map_init_ab_params{
    np01_xy_vec xy_vec_a;
    np01_xy_vec xy_vec_b;
    np01_float64 loc_grid_density; /* points (a+b) / grid square */
    size_t max_loc_grid_sq_count; /* target max number of grid squares */
};


class np01_xy_pair_map{
private:
    np01_xy_node_map *m_node_map_a;
    np01_xy_node_map *m_node_map_b;
    np01_xy_edge_map *m_edge_map;
public:
    np01_xy_pair_map();
   ~np01_xy_pair_map();
    void init_ab( const np01_xy_pair_map_init_ab_params *init_params);
    void init_edges_ctr_out( const np01_float64 dsq_noise_ampl );
    void init_edges_far_near();

    np01_xy_node_map *get_node_map_a() const { return m_node_map_a; }
    np01_xy_node_map *get_node_map_b() const { return m_node_map_b; }
    np01_xy_edge_map *get_edge_map() const { return m_edge_map; }

    int verify_data( char *err_msg, const size_t err_msg_capacity,
                     size_t *err_msg_pos ) const;
    std::ostream& ostream_output(std::ostream& os) const;
    void write_bmp_file(const char *file_name) const;
};

/* cost (criteria) for genetic algorithm sub-solution */
struct np01_xy_pair_sub_soln_cost{
public:
    np01_uint32 m_unpaired_count; /* unpaired A+B. zero if only A or only B */
    np01_float64 m_viol_edge_len_sq_sum; /* sum(len^2) for violating edges */
    np01_float64 m_edge_len_sum; /* sum of edge lengths */
public:
    np01_xy_pair_sub_soln_cost():m_unpaired_count(0),
        m_viol_edge_len_sq_sum(0.0), m_edge_len_sum(0.0){}
    np01_xy_pair_sub_soln_cost(const np01_xy_pair_sub_soln_cost& master):
        m_unpaired_count(master.m_unpaired_count), 
        m_viol_edge_len_sq_sum(master.m_viol_edge_len_sq_sum),
        m_edge_len_sum(master.m_edge_len_sum){}
   ~np01_xy_pair_sub_soln_cost(){}
    np01_xy_pair_sub_soln_cost& operator=(const np01_xy_pair_sub_soln_cost& 
        master){
        m_unpaired_count=master.m_unpaired_count;
        m_viol_edge_len_sq_sum=master.m_viol_edge_len_sq_sum;
        m_edge_len_sum=master.m_edge_len_sum;
        return *this; }
    bool operator==(const np01_xy_pair_sub_soln_cost& other) const{
        return ((m_unpaired_count==other.m_unpaired_count) &&
            (m_viol_edge_len_sq_sum==other.m_viol_edge_len_sq_sum) &&
            (m_edge_len_sum==other.m_edge_len_sum)); }
    bool operator!=(const np01_xy_pair_sub_soln_cost& other) const{
        return ((m_unpaired_count!=other.m_unpaired_count) ||
            (m_viol_edge_len_sq_sum!=other.m_viol_edge_len_sq_sum) ||
            (m_edge_len_sum!=other.m_edge_len_sum)); }
    bool operator<(const np01_xy_pair_sub_soln_cost& other) const{
        bool result = false;
        if(m_unpaired_count<other.m_unpaired_count){result=true;}
        else if(other.m_unpaired_count<m_unpaired_count){result=false;}
        else if(m_viol_edge_len_sq_sum<other.m_viol_edge_len_sq_sum){
            result=true;}
        else if(other.m_viol_edge_len_sq_sum<m_viol_edge_len_sq_sum){
            result=false;}
        else if(m_edge_len_sum<other.m_edge_len_sum){result=true;}
        else{result=false;}
        return result; }
    void reset(){m_unpaired_count=0; m_viol_edge_len_sq_sum=0.0;
        m_edge_len_sum=0.0;}
};

/* genetic algorithm sub-solution (or complete solution) */
class np01_xy_pair_sub_soln{
private:
    np01_uint32_pair_vec m_ab_idx_pair_vec;
    np01_float64_vec m_edge_len_vec;
    np01_xy_pair_sub_soln_cost m_cost;
    bool m_cost_valid;
    np01_xy_pair_sub_soln *m_free_chain_next;
public:
    np01_xy_pair_sub_soln(): m_ab_idx_pair_vec(),m_edge_len_vec(),m_cost(),
        m_cost_valid(false),m_free_chain_next(NULL){}
   ~np01_xy_pair_sub_soln(){}                  
    void reset(){m_ab_idx_pair_vec.clear(); m_edge_len_vec.clear();
        m_cost.reset(); m_cost_valid=false; m_free_chain_next=NULL; }
    void reserve(const size_t sz){m_ab_idx_pair_vec.reserve(sz);
        m_edge_len_vec.reserve(sz); }
    void add_ab_idx_pair(const np01_uint32_pair& ab_idx_pair){
       m_ab_idx_pair_vec.push_back(ab_idx_pair); }
    void set_ab_idx_pair(const size_t i, const np01_uint32_pair& ab_idx_pair){
         if(i < m_ab_idx_pair_vec.size()){
             m_ab_idx_pair_vec.at(i) = ab_idx_pair; }
         }
    void sort_unique_ab_idx_pairs();
    void set_len(const np01_uint32& i, const np01_float64& edge_len){
        if(m_edge_len_vec.size() >= i){m_edge_len_vec.resize(i+1);}
        m_edge_len_vec.at(i) = edge_len; }
    void set_cost(const np01_xy_pair_sub_soln_cost& cost){m_cost=cost;}
    void set_cost_valid(const bool& valid){m_cost_valid = valid;}
    void set_free_chain_next(np01_xy_pair_sub_soln *n){m_free_chain_next=n;}

    size_t get_ab_idx_pair_count() const{ return m_ab_idx_pair_vec.size(); }
    const np01_uint32_pair& get_ab_idx_pair(const size_t& i) const{
        return m_ab_idx_pair_vec.at(i); }
    size_t get_edge_len_count() const{ return m_edge_len_vec.size(); }
    const np01_float64& get_edge_len(const size_t& i) const{
        return m_edge_len_vec.at(i); }
    const np01_xy_pair_sub_soln_cost& get_cost() const { return m_cost; }
    bool is_cost_valid() const{return m_cost_valid;}
    np01_xy_pair_sub_soln *get_free_chain_next() const {
        return m_free_chain_next; }
    std::ostream& ostream_output(std::ostream& os) const;
};

/* genetic algorithm initialization parameters */
struct np01_xy_pair_genetic_wksp_init_params{
public:
    np01_xy_pair_map *m_xy_pair_map;
    np01_float64 m_target_max_edge_len;
    np01_uint32 m_max_iteration_count;
    np01_uint32 m_no_change_max_iteration_count;
    np01_uint32 m_inner_loop_max_iteration_count;
public:
    void set_default(np01_xy_pair_map *xy_pair_map);
};

struct np01_xy_pair_genetic_wksp_result{
public:
    np01_xy_pair_sub_soln m_solution;
    np01_uint32_vec m_edge_crossover_count_vec;
    np01_uint32 m_total_edge_crossover_count;
    np01_uint32 m_len_viol_count;
    np01_float64 m_max_result_length; /* actual maximum A->B length */
};

class np01_xy_pair_genetic_wksp{
private:
    typedef std::vector<np01_xy_pair_sub_soln *>soln_ptr_vec;
    typedef soln_ptr_vec::iterator soln_ptr_vec_itr;
    typedef soln_ptr_vec::const_iterator soln_ptr_vec_citr;
    typedef std::set<np01_uint32> np01_uint32_set;
    typedef np01_uint32_set::iterator np01_uint32_set_itr;
    typedef np01_uint32_set::const_iterator np01_uint32_set_citr;
private:
    np01_xy_pair_map *m_xy_pair_map; /* target data set */

    /* control parameters */
    np01_float64 m_target_max_edge_len;
    np01_uint32 m_max_iteration_count;
    np01_uint32 m_no_change_max_iteration_count;
    np01_uint32 m_inner_loop_max_iteration_count;
    
    np01_uint32_vec m_usage_count_a_vec; /* count: node A[i] was processed */
    np01_uint32_vec m_usage_count_b_vec; /* count: node B[i] was processed */
    np01_uint32 m_a_idx;
    np01_uint32 m_b_idx;
    np01_uint32 m_node_a_count;
    np01_uint32 m_node_b_count;
    np01_uint32 m_a_offset; /* randomized offset */
    np01_uint32 m_b_offset; /* randomized offset */
    np01_uint32 m_target_usage_count;  /* looking for node with */

    /*utility containers */
    np01_xy_node_map::node_vec m_circle_a_nodes;
    np01_xy_node_map::node_vec m_circle_b_nodes;
    np01_xy_edge_map::edge_vec m_near_edges;
    np01_xy_edge_map::idx_edge_ptr_pair_vec m_utility_vec;
    np01_uint32_set m_unused_idx_set_a;
    np01_uint32_set m_unused_idx_set_b;

    soln_ptr_vec m_soln_vec; /* population of solutions of subset of nodes */
    np01_xy_pair_sub_soln *m_soln_free_chain; /* soln objects ready for use */

    bool m_changed; /* data set was changed during middle loop */

    /* loop counts */
    np01_uint32 m_iteration_count;
    np01_uint32 m_no_change_iteration_count;
    np01_uint32 m_inner_loop_iteration_count;

    np01_uint32 m_rand_uint32; /* pseudo-random number*/

public:
    np01_xy_pair_genetic_wksp();
    np01_xy_pair_genetic_wksp(const np01_xy_pair_genetic_wksp_init_params&
        init_params);
   ~np01_xy_pair_genetic_wksp();
    void init(const np01_xy_pair_genetic_wksp_init_params& init_params);
    void run();
    void get_result(np01_xy_pair_genetic_wksp_result *result) const;
private:
    void advance_a_idx();
    void advance_b_idx();
    void seed_soln_vec();
    void mutate();
    void mutate2(np01_xy_pair_sub_soln *s);
    void mutate3(np01_xy_pair_sub_soln *s);
    void mutate4(np01_xy_pair_sub_soln *s);
    void gen_crossover();
    void cull();
    void save_best_soln();
    np01_xy_pair_sub_soln *alloc_sub_soln();
    void free_sub_soln(np01_xy_pair_sub_soln *s);

    void free_soln_vec();
    void free_free_chain();
    void compute_len_cost(np01_xy_pair_sub_soln *s) const;
    void advance_rand();
    int get_rand_int() const { return static_cast<int>(
        (static_cast<unsigned int>(m_rand_uint32 >> 16)) & RAND_MAX); }
    
#if 0

Genetic Min Distance Assignment Algorithm 
input:
  vector<xy> points_a  /* smaller set of points */
  vector<xy> points_b  /* larger set of points */
  double     max_allowed_distance
output
  vector<pair<int, int> > assign_idx_pairs
  double total_dist_passing  /* sum of distances between paired points, each distance <= max_allowed_distance */
  double total_dist_failing  /* sum of distances between paired points, each distance > max_allowed_distance */
  int unassigned_count

Cost Function
  input vector<pair<int, int> > assign_idx_pairs   /* full set or partial set */
output
  double total_dist_passing  /* sum of distances between paired points, each distance <= max_allowed_distance */
  double total_dist_failing  /* sum of distances between paired points, each distance > max_allowed_distance */
  int unassigned_count

Workspace
  vector<xy> points_a  /* smaller set of points */
  vector<xy> points_b  /* larger set of points */
  double     max_allowed_distance

  stopping criteria

  double locator_grid_square_size
  xy locator_grid_center
  vector< node *> locator_grid_a
  vector< node *> locator_grid_b
  vector<int> processed_count a

  /* result = best complete solution */
  vector<pair<int, int> > assign_idx_pairs
  double total_dist_passing  /* sum of distances between paired points, each distance <= max_allowed_distance */
  double total_dist_failing  /* sum of distances between paired points, each distance > max_allowed_distance */
  int unassigned_count

  linked list (free chain) of unused candidate solutions
 
  linked list candidate solutions in current genetic algorithm

  xy, radius of current genetic algorithm 


Procedure
  Initialize locator grid
  Initial Assignment
    for each point a, assign to nearest unassignd point b
  While !done
    find least processed point a
    select random radius
    get all pairs<a,b> such that one or both points are within circle
    run genetic algorithm on set
      create a population of candidate solutions from the set
        copy of initial set
        random mutations (swaps) of initial set 
      score each solution
      for n iterations
         remove solutions
           keep best solution
           select other solution randomly, preferring to remove worse solutions
         create new solutions by mutation  
      merge best solution into the complete solution
        if no change
      update processed count for each point a
    if iteration count >= max iterations,
      done
    if iteration count >= no_change threshold A and no-change iteration count B
       done
#endif
};


class np01_main{
public:
    static int run(int argc, char *argv[]);
private:
    static const char *m_prog_name;
    static const char *m_version_str;
    int m_argc;
    const char **m_argv;
    bool m_test_option;
    int m_test_number;
    bool m_pair_points_option;
    bool m_xy_a_file_option;
    std::string m_xy_a_filename;
    bool m_xy_b_file_option;
    std::string m_xy_b_filename;
    bool m_out_ab_file_option;
    std::string m_out_ab_filename;
    bool m_max_len_option;
    double m_max_len;
    bool m_iterations_option;
    int m_iterations;
    bool m_test_iterations_option;
    int m_test_iterations;
    bool m_test_rand_seed_option;
    int m_test_rand_seed;
public:
    np01_main(int argc, char *argv[]);
   ~np01_main();
    int execute();
private:
    void parse_cmd_line();
    void execute_pp();
    void execute_test_0();
    void execute_test_1();
    void execute_test_2();
    void execute_test_3();
    void execute_test_4();
    void execute_test_5();
    void test_init_xy_vec_method_0(const np01_uint32& approx_point_count,
        const np01_float64& typical_point_spacing, 
        const np01_uint32& rand_seed, np01_xy_vec *xy_vec);
    void test_init_xy_vec_method_1(const np01_uint32& approx_point_count,
        const np01_float64& typical_point_spacing, 
        const np01_uint32& rand_seed, np01_xy_vec *xy_vec);
    void test_init_xy_vec_method_2(const np01_uint32& approx_point_count,
        const np01_float64& typical_point_spacing, 
        const np01_uint32& rand_seed, np01_xy_vec *xy_vec);


};


void np01_snprintf( char *buf, const size_t buf_capacity, size_t *buf_pos,
               const char *fmt, ... );

} /* namespace np01 */

#endif /* NP01_H */
