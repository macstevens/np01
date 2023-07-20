/* np01.cpp Newport Algorithms version 1

Reference: 
    https://www.ifte.de/mitarbeiter/lienig/date2008pp837-842.pdf
      Pin Assignment Algorithm minimizing total length and crossovers.

Reference: 
    https://ieeexplore.ieee.org/document/5351234
    https://www.researchgate.net/publication/224089030_Pin_assignment_for_wire_length_minimization_after_floorplanning_phase
    https://www.semanticscholar.org/paper/Pin-assignment-for-wire-length-minimization-after-He-Dong/cbd8b3de85ee95d91a8bd77e6924ae2da7fe032d
    Pin assignment for wire length minimization after floorplanning phase
    DOI:10.1109/ASICON.2009.5351234
    Authors: Xu He, Sheqin Dong
    Minimization of Total Wire Length, simulated annealing

Reference:
    http://www.gstitt.ece.ufl.edu/courses/fall19/eel4720_5721/reading/Routing.pdf
    "minimize the total length of interconnect required."

Reference:
    https://en.wikipedia.org/wiki/Grid_(spatial_index)

Reference:
    https://gistbok.ucgis.org/bok-topics/spatial-indexing
    "decompose 2D plane into a list of cells ... Fixed grid index is an n×n 
    array of equal-size cells. Each one is associated with a list of spatial
    objects which intersect or overlap with the cell."

Copyright (c) 2023 Mac Stevens <stevensm@earthlink.net> <www.macstevens.net>

Permission to use, copy, modify, and distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

Reference: https://opensource.org/licenses/ISC
*/

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <limits>
#include <set>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "np01.h"
#include "np01_bmp.h"


#define CF01_SUPPORT ( 1 )
#if defined( CF01_SUPPORT )
#include "cf01.h"
    #define AA_INCR_CALL_DEPTH()        CF01_AA_INCR_CALL_DEPTH()
    #define AA_DECR_CALL_DEPTH()        CF01_AA_DECR_CALL_DEPTH()
    #define AUTO_ASSERT( _condition )   CF01_AUTO_ASSERT(_condition)
    #define AA_ALWAYS_ASSERT( _condition ) \
                  CF01_AA_XDBG_ASSERT((_condition), CF01_AA_DEBUG_LEVEL_0)
    #define AA_SHOULD_RUN_XDBG(_dbg_lvl) CF01_AA_SHOULD_RUN_XDBG(_dbg_lvl)
    #define AA_XDBG_ASSERT( _condition, _dbg_lvl ) \
                            CF01_AA_XDBG_ASSERT( (_condition), _dbg_lvl )
    #define AA_ERR_BUF()                CF01_AA_ERR_BUF()
    #define AA_ERR_BUF_CAPACITY()       CF01_AA_ERR_BUF_CAPACITY()
    #define AA_ERR_BUF_POS_PTR()        CF01_AA_ERR_BUF_POS_PTR()
#else
    #define AA_INCR_CALL_DEPTH()        
    #define AA_DECR_CALL_DEPTH()        
    #define AUTO_ASSERT( _condition )   assert(_condition)
    #define AA_ALWAYS_ASSERT( _condition )   assert(_condition)
    #define AA_SHOULD_RUN_XDBG(_dbg_lvl)  (0)
    #define AA_XDBG_ASSERT( _condition, _dbg_lvl ) 
    #define AA_ERR_BUF()                (NULL)
    #define AA_ERR_BUF_CAPACITY()       (0)
    #define AA_ERR_BUF_POS_PTR()        (NULL)
#endif


namespace np01 {

#define NP01_DEFAULT_MAX_DOT_CROSS_ERR (1e-8)
#define NP01_DEFAULT_MAX_FWD_ERR (1e-8)
#define NP01_GEN_MAX_POPULATION_SIZE (12)


np01_xy_node::np01_xy_node(): m_idx(0), m_owner(NULL), m_edge(NULL),
    m_loc_grid_i(0), m_loc_grid_j(0), m_loc_grid_prev(NULL),
    m_loc_grid_next(NULL), m_x(0.0), m_y(0.0),
    m_xy_group_type(NP01_XY_GROUP_TYPE_UNKNOWN)
{}

np01_xy_node::~np01_xy_node(){
AA_INCR_CALL_DEPTH();
if(NULL != m_edge){
    disconnect_edge();}
if(NULL != m_owner){
    AA_ALWAYS_ASSERT(false);} /* update owner locator grid, etc. */
AA_DECR_CALL_DEPTH();
}

int np01_xy_node::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
const np01_xy_edge_map *edge_owner=(NULL==m_edge) ? NULL : m_edge->get_owner();
const np01_xy_pair_map *edge_owner_owner = (NULL==edge_owner) ?
    NULL : edge_owner->get_owner();
const np01_xy_pair_map *owner_owner = (NULL==m_owner) ?
    NULL : m_owner->get_owner();
if( NULL != m_owner ) {
    const np01_xy_loc_grid_dim& loc_grid_dim = m_owner->get_loc_grid_dim();
    const np01_float64 loc_grid_i_check = loc_grid_dim.get_i( m_x );
    const np01_float64 loc_grid_j_check = loc_grid_dim.get_j( m_y );

    if( m_idx >= m_owner->get_node_count() ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_idx(%i) >= m_owner->get_node_count(%i)\n",
            m_idx, m_owner->get_node_count()); }
    else if( this != m_owner->get_node_by_idx(m_idx) ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != m_owner->get_node_by_idx(%i)=%x\n",
            this, m_idx, m_owner->get_node_by_idx(m_idx)); }
    if(m_owner->get_xy_group_type() != m_xy_group_type) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_owner->get_xy_group_type()=%i != m_xy_group_type=%i\n",
            m_owner->get_xy_group_type(), m_xy_group_type); }
    if((loc_grid_dim.get_w() == 0) || (loc_grid_dim.get_h() == 0) ||
        (loc_grid_dim.get_sq_size() <= 0.0)){
        /*m_loc_grid_i should be 0*/
        /*m_loc_grid_j should be 0*/
        /*loc_grid_prev should be NULL*/
        /*loc_grid_next should be NULL*/
        }
    else{
        if(m_loc_grid_i >= loc_grid_dim.get_w()) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "m_loc_grid_i=%i >= loc_grid_dim.get_w()=%i\n",
                m_loc_grid_i, loc_grid_dim.get_w()); }
        if(m_loc_grid_i != loc_grid_i_check) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "m_loc_grid_i=%i>=loc_grid_i_check=%i:x=%f:x_min=%f:sqsz=%f\n",
                m_loc_grid_i, loc_grid_i_check, m_x, loc_grid_dim.get_x_min(),
                loc_grid_dim.get_sq_size()); }
        if(m_loc_grid_j >= loc_grid_dim.get_h()) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "m_loc_grid_j=%i >= loc_grid_dim.get_h()=%i\n",
                m_loc_grid_j, loc_grid_dim.get_h()); }
        if(m_loc_grid_j != loc_grid_j_check) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "m_loc_grid_j=%i>=loc_grid_j_check=%i:y=%f:y_min=%f:sqsz=%f\n",
                m_loc_grid_j, loc_grid_j_check, m_y, loc_grid_dim.get_y_min(),
                loc_grid_dim.get_sq_size()); }
        }

    if((owner_owner != edge_owner_owner) && (NULL != m_edge)) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "((this=%x)->owner=%x)->owner=%x != "
            "((m_edge=%x)->owner=%x)->owner=%x\n",  this, m_owner, owner_owner,
            m_edge, edge_owner, edge_owner_owner); }
    }

if(NULL != m_edge) {
    if(NP01_XY_GROUP_TYPE_A == m_xy_group_type ) {
        if(this != m_edge->get_node_a()) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "this=%x != (edge=%x)->node_a=%x\n",
                this, m_edge, m_edge->get_node_a());
            }
        }
    else if(NP01_XY_GROUP_TYPE_B == m_xy_group_type ) {
        if(this != m_edge->get_node_b()) {
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "this=%x != (edge=%x)->node_b=%x\n",
                this, m_edge, m_edge->get_node_b());
            }
        }
    else{
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x  edge=%x  xy_group_type=%i\n",
            this, m_edge, m_xy_group_type);
        }
    }

if(NULL != m_loc_grid_prev)
    {
    if(this == m_loc_grid_prev) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this == m_loc_grid_prev=%x\n", this); }
    if(this != m_loc_grid_prev->m_loc_grid_next) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != m_loc_grid_prev(%x)->m_loc_grid_next=%x\n", 
            this, m_loc_grid_prev, m_loc_grid_prev->m_loc_grid_next); }
    if(m_idx == m_loc_grid_prev->m_idx) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_idx == m_loc_grid_prev->m_idx=%i\n", m_idx); }
    if(m_owner != m_loc_grid_prev->m_owner) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_owner=%x != m_loc_grid_prev->m_owner=%x\n", 
            m_owner, m_loc_grid_prev->m_owner); }
    if((NULL!=m_edge)&&(NULL!=m_loc_grid_prev->m_edge)&&
       (m_edge==m_loc_grid_prev->m_edge)){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_edge=%x==m_loc_grid_prev(%x)->m_edge\n", 
            m_edge, m_loc_grid_prev->m_edge); }
    if(m_loc_grid_i != m_loc_grid_prev->m_loc_grid_i) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_loc_grid_i=%i != m_loc_grid_prev->m_loc_grid_i=%i\n", 
            m_loc_grid_i, m_loc_grid_prev->m_loc_grid_i); }
    if(m_loc_grid_j != m_loc_grid_prev->m_loc_grid_j) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_loc_grid_j=%i != m_loc_grid_prev->m_loc_grid_j=%i\n", 
            m_loc_grid_j, m_loc_grid_prev->m_loc_grid_j); }
    }

if(NULL != m_loc_grid_next)
    {
    if(this == m_loc_grid_next) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this == m_loc_grid_next=%x\n", this); }
    if(this != m_loc_grid_next->m_loc_grid_prev) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != m_loc_grid_next(%x)->m_loc_grid_prev=%x\n", 
            this, m_loc_grid_next, m_loc_grid_next->m_loc_grid_prev); }
    if(m_idx == m_loc_grid_next->m_idx) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_idx == m_loc_grid_next->m_idx=%i\n", m_idx); }
    if(m_owner != m_loc_grid_next->m_owner) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_owner=%x != m_loc_grid_next->m_owner=%x\n", 
            m_owner, m_loc_grid_next->m_owner); }
    if((NULL!=m_edge)&&(NULL!=m_loc_grid_next->m_edge)&&
       (m_edge==m_loc_grid_next->m_edge)){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_edge=%x==m_loc_grid_next(%x)->m_edge\n", 
            m_edge, m_loc_grid_next->m_edge); }
    if(m_loc_grid_i != m_loc_grid_next->m_loc_grid_i) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_loc_grid_i=%i != m_loc_grid_next->m_loc_grid_i=%i\n", 
            m_loc_grid_i, m_loc_grid_next->m_loc_grid_i); }
    if(m_loc_grid_j != m_loc_grid_next->m_loc_grid_j) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_loc_grid_j=%i != m_loc_grid_next->m_loc_grid_j=%i\n", 
            m_loc_grid_j, m_loc_grid_next->m_loc_grid_j); }
    }

return err_cnt;
}

std::ostream& np01_xy_node::ostream_output(std::ostream& os) const{
os << "<xy_node=" << std::hex << this << ">\n";
os << "<idx>" << std::dec << m_idx << "</idx>\n";
os << "<owner>" << std::hex << m_owner << "</owner>\n";
os << "<edge>" << m_edge << std::dec << "</edge>\n";
os << "<i>" << m_loc_grid_i << "</i><j>" << m_loc_grid_j << "</j>\n";
os << "<loc_grid_prev>" << std::hex << m_loc_grid_prev << "</loc_grid_prev>\n";
os << "<loc_grid_next>" << m_loc_grid_next << std::dec << "</loc_grid_next>\n";
os << "<x>" << m_x << "</x><y>" << m_x << "</y>\n";
os << "<xy_group_type>" << m_xy_group_type << "</xy_group_type>\n";
os << "</xy_node>\n";
return os;
}

void np01_xy_node::disconnect_edge(){
AA_INCR_CALL_DEPTH();
if( NULL != m_edge ){
    switch(m_xy_group_type){
        case NP01_XY_GROUP_TYPE_A:
            AUTO_ASSERT(this == m_edge->get_node_a());
            m_edge->set_node_a(NULL);
            break;
        case NP01_XY_GROUP_TYPE_B:
            AUTO_ASSERT(this == m_edge->get_node_b());
            m_edge->set_node_b(NULL);
            break;
        case NP01_XY_GROUP_TYPE_UNKNOWN: /* fall through */
        default:
            AA_ALWAYS_ASSERT(false);
            if(this == m_edge->get_node_a()){
                m_edge->set_node_a(NULL);}
            if(this == m_edge->get_node_b()){
                m_edge->set_node_b(NULL);}
            break;
        }
    set_edge(NULL);
    }
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
}


void np01_xy_loc_grid_dim::init(
    const np01_xy_loc_grid_dim_init_params *init_params){
AA_INCR_CALL_DEPTH();
if(NULL == init_params){
    AA_ALWAYS_ASSERT(false);
    reset();
    }
else if((init_params->point_count == 0) ||
    (init_params->loc_grid_density <= 0.0) ||
    (init_params->max_loc_grid_sq_count == 0) ||
    ((init_params->bb_min_xy).first > (init_params->bb_max_xy).first) ||
    ((init_params->bb_min_xy).second > (init_params->bb_max_xy).second)){
    AA_ALWAYS_ASSERT(init_params->point_count > 0);
    AA_ALWAYS_ASSERT(init_params->loc_grid_density > 0.0);
    AA_ALWAYS_ASSERT(init_params->max_loc_grid_sq_count > 0);
    AA_ALWAYS_ASSERT((init_params->bb_min_xy).first <= (init_params->bb_max_xy).first);
    AA_ALWAYS_ASSERT((init_params->bb_min_xy).second <= (init_params->bb_max_xy).second);
    reset();
    }
else{
    /* bounding box */
    const np01_float64& x_min = (init_params->bb_min_xy).first;
    const np01_float64& x_max = (init_params->bb_max_xy).first;
    const np01_float64& y_min = (init_params->bb_min_xy).second;
    const np01_float64& y_max = (init_params->bb_max_xy).second;
    const np01_float64 bb_dx = x_max-x_min;
    const np01_float64 bb_dy = y_max-y_min;
    const np01_float64 bb_area = bb_dx*bb_dy;
    AA_ALWAYS_ASSERT(bb_dx >= 0.0);
    AA_ALWAYS_ASSERT(bb_dy >= 0.0);
    AA_ALWAYS_ASSERT(bb_area >= 0.0);

    /* estimate number of grid squares */
    np01_float64 grid_sq_count_float =
       static_cast<np01_float64>(init_params->point_count)/
       (init_params->loc_grid_density);
    const np01_float64 max_grid_sq_count_float =
        static_cast<np01_float64>(init_params->max_loc_grid_sq_count);
    if(grid_sq_count_float > max_grid_sq_count_float){
        grid_sq_count_float = max_grid_sq_count_float; }
    if( grid_sq_count_float < 1.0 ){
        grid_sq_count_float = 1.0; }

    /* estimate grid square size */
    np01_float64 sq_sz_approx = 1.0;
    if( bb_area > 0.0 ){
        const np01_float64 sq_area_approx = bb_area/grid_sq_count_float;
        AA_ALWAYS_ASSERT(sq_area_approx > 0.0);
        sq_sz_approx = sqrt(sq_area_approx);
        }
    else if( bb_dx > 0.0 ){
        sq_sz_approx = bb_dx / grid_sq_count_float;
        }
    else if( bb_dy > 0.0 ){
        sq_sz_approx = bb_dy / grid_sq_count_float;
        }
    else{
        sq_sz_approx = 1.0;
        }

    /* round up to an integer/(2^k) fraction */
    int n, f33;
    double f, ff;
    f = frexp( sq_sz_approx, &n );
    AA_ALWAYS_ASSERT(f >= 0.5);
    f33 = 1 + static_cast<int>(32.0 * f);
    ff = static_cast<double>(f33)/32.0;
    AA_ALWAYS_ASSERT(ff <= 1.0);
    AA_ALWAYS_ASSERT(ff >= 0.5);
    set_sq_size(ldexp(ff, n));
    AA_ALWAYS_ASSERT(get_sq_size() >= sq_sz_approx);
    AA_ALWAYS_ASSERT(get_sq_size() > 0.0);

    /* loc grid width, height*/
    static const np01_float64 max_wh_float = static_cast<np01_float64>(
        std::numeric_limits<np01_uint16>::max() - 1);
    np01_float64 w_float = 1.0 + (bb_dx / get_sq_size());
    np01_float64 h_float = 1.0 + (bb_dy / get_sq_size());
    if( w_float > max_wh_float ){ w_float = max_wh_float; }
    if( h_float > max_wh_float ){ h_float = max_wh_float; }
    set_w( static_cast<np01_uint16>(w_float));
    set_h( static_cast<np01_uint16>(h_float));

    /* loc grid lower left corner*/
    const np01_float64 bb_ctr_x = (x_min+x_max)/2.0;
    const np01_float64 bb_ctr_y = (y_min+y_max)/2.0;
    const np01_float64 gw = static_cast<np01_float64>(get_w())*get_sq_size();
    const np01_float64 gh = static_cast<np01_float64>(get_h())*get_sq_size();
    AUTO_ASSERT(gw >= bb_dx);
    AUTO_ASSERT(gh >= bb_dy);
    set_x_min(bb_ctr_x-(gw/2.0));
    set_y_min(bb_ctr_y-(gh/2.0));
    AUTO_ASSERT(get_x_min() <= x_min);
    AUTO_ASSERT(get_y_min() <= y_min);
    AUTO_ASSERT((get_x_min()+gw) >= x_max);
    AUTO_ASSERT((get_y_min()+gh) >= y_max);
    }
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
}

void np01_xy_loc_grid_dim::get_bb_indices( const np01_xy& xy_min, 
    const np01_xy& xy_max, np01_uint16_pair *ij_min,
    np01_uint16_pair *ij_max) const{
if( NULL != ij_min ) {
    ij_min->first = get_i(xy_min.first);
    ij_min->second = get_j(xy_min.second);}
if( NULL != ij_max ) {
    ij_max->first = get_i(xy_max.first);
    ij_max->second = get_j(xy_max.second);}
}

int np01_xy_loc_grid_dim::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
if(((m_w==0) && (m_h!=0)) || ((m_w!=0) && (m_h==0))){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "w=%i != h=%i\n", m_w, m_h); }
if((m_sq_size < 0.0) || ((m_w!=0) && (m_h!=0) && (m_sq_size <= 0.0))){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "w=%i  h=%i  sq_sz=%f\n", m_w, m_h, m_sq_size); }
return err_cnt;
}

std::ostream& np01_xy_loc_grid_dim::ostream_output(std::ostream& os) const{
os << "<xy_loc_grid_dim>\n";
os << "<w>" << m_w << "</w><h>" << m_h << "</h>\n";
os << "<x_min>" << m_x_min << "</x_min><y_min>" << m_y_min << "</y_min>\n";
os << "<sq_size>" << m_sq_size << "</sq_size>\n";
os << "</xy_loc_grid_dim>\n";
return os;
}

np01_xy_node_map::np01_xy_node_map(): m_owner(NULL), m_loc_grid_dim(),
    m_loc_grid(), m_node_vec(), m_xy_group_type(NP01_XY_GROUP_TYPE_UNKNOWN){}

np01_xy_node_map::~np01_xy_node_map(){
AA_INCR_CALL_DEPTH();
clear_loc_grid();
free_nodes();
AA_DECR_CALL_DEPTH();
}


void np01_xy_node_map::init_loc_grid( np01_xy_loc_grid_dim& d ){
AA_INCR_CALL_DEPTH();
size_t loc_grid_sz;
if(!m_loc_grid.empty()){
    clear_loc_grid();}
m_loc_grid_dim = d;
loc_grid_sz = static_cast<size_t>(d.get_w()) * static_cast<size_t>(d.get_h());
m_loc_grid.resize(loc_grid_sz, NULL);
if(!m_node_vec.empty()){
    node_vec_citr n_itr=m_node_vec.begin();
    for(;n_itr!= m_node_vec.end(); ++n_itr){
        np01_xy_node *n = *n_itr;
        AA_ALWAYS_ASSERT(NULL != n);
        insert_node_in_loc_grid(n);
        }
    }
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_DECR_CALL_DEPTH();
}

np01_xy_node *np01_xy_node_map::create_node( const np01_xy& xy ){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
np01_xy_node *n = new np01_xy_node();
n->set_idx(static_cast<np01_uint32>(m_node_vec.size()));
m_node_vec.push_back(n);
n->set_owner(this);
AA_ALWAYS_ASSERT(n->get_edge() == NULL);
n->set_xy(xy.first, xy.second);
n->set_xy_group_type(m_xy_group_type);
insert_node_in_loc_grid(n);
AA_DECR_CALL_DEPTH();
AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
return n;
}

void np01_xy_node_map::get_nodes_in_bb(const np01_xy& xy_min,
    const np01_xy& xy_max, node_vec *nodes) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AUTO_ASSERT(xy_min.first <= xy_max.first);
AUTO_ASSERT(xy_min.second <= xy_max.second);
AA_ALWAYS_ASSERT(NULL != nodes);
size_t lg_idx;
np01_uint16 i, j;
np01_uint16_pair ij_min, ij_max;
m_loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
for(i = ij_min.first; i <= ij_max.first; ++i){
    for(j = ij_min.second; j <= ij_max.second; ++j){
        lg_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
            static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid.size());
        if(lg_idx < m_loc_grid.size()){
            np01_xy_node *node = m_loc_grid[lg_idx];
            while(NULL != node){
                const np01_float64& x = node->get_x();
                const np01_float64& y = node->get_y();
                if((xy_min.first <= x) && (x <= xy_max.first) &&
                    (xy_min.second <= y) && (y <= xy_max.second) ){
                    nodes->push_back(node);
                    }
                node = node->get_loc_grid_next();
                }
            }
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_node_map::get_nodes_in_circle(const np01_xy& xy_ctr, const np01_float64& r,
        node_vec *nodes) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(r >= 0.0);
AA_ALWAYS_ASSERT(NULL != nodes);
const np01_xy xy_min(xy_ctr.first - r, xy_ctr.second - r);
const np01_xy xy_max(xy_ctr.first + r, xy_ctr.second + r);
const np01_float64 r_sq = r * r;
size_t lg_idx;
np01_uint16 i, j;
np01_uint16_pair ij_min, ij_max;
m_loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
for(i = ij_min.first; i <= ij_max.first; ++i){
    for(j = ij_min.second; j <= ij_max.second; ++j){
        lg_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
            static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid.size());
        if(lg_idx < m_loc_grid.size()){
            np01_xy_node *node = m_loc_grid[lg_idx];
            while(NULL != node){
                const np01_float64& x = node->get_x();
                const np01_float64& y = node->get_y();
                const np01_float64 dx = x - xy_ctr.first;
                const np01_float64 dy = y - xy_ctr.second;
                const np01_float64 d_sq = (dx * dx) + (dy * dy);
                if( d_sq <= r_sq ){
                    nodes->push_back(node); }
                node = node->get_loc_grid_next();
                }
            }
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_node_map::get_all_nodes(node_vec *nodes) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != nodes);
nodes->reserve(nodes->size() + m_node_vec.size());
nodes->insert(nodes->end(), m_node_vec.begin(), m_node_vec.end());
AA_DECR_CALL_DEPTH();
}

/* xy: search center,  d: output parameter, distance from xy to near node
utility_vec: used during search, then left empty. */
np01_xy_node *np01_xy_node_map::get_near_node(const np01_xy& xy, np01_float64 *d,
    node_vec *utility_vec) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);

np01_xy_node *near_node = NULL;
np01_float64 near_dsq = 0.0;
if(!m_node_vec.empty()){
    /* utility vector */
    node_vec *search_vec = (NULL==utility_vec) ?  new node_vec() : utility_vec;
    search_vec->clear();

    /* search radius */
    const np01_float64 w = static_cast<np01_float64>(m_loc_grid_dim.get_w());
    const np01_float64 h = static_cast<np01_float64>(m_loc_grid_dim.get_h());
    const np01_float64& sq_sz = m_loc_grid_dim.get_sq_size();
    const np01_float64 min_search_r = sq_sz / 2.0;
    const np01_float64 max_search_r = sq_sz * (w + h) / 2.0;
    np01_float64 search_r = min_search_r;

    /* search */
    enum search_mode_type {srch_md_circle, srch_md_all, srch_md_done};
    search_mode_type search_mode = ((min_search_r>0.0) && (max_search_r>0.0)) ?
        srch_md_circle : srch_md_all;
    while(srch_md_done != search_mode){
        AA_ALWAYS_ASSERT(NULL != search_vec);
        AUTO_ASSERT(search_vec->empty());
        AUTO_ASSERT(NULL == near_node);
        /* search in circle or get all nodes */
        switch(search_mode){
            case srch_md_circle:
                get_nodes_in_circle(xy,search_r,search_vec);
                break;
            case srch_md_all:
                get_all_nodes(search_vec);
                break;
            case srch_md_done:
            default:
                AA_ALWAYS_ASSERT(false);
                break;
            }

        /* find closest node */
        near_dsq = std::numeric_limits<np01_float64>::max();
        node_vec_citr node_itr = search_vec->begin();
        for(; node_itr != search_vec->end(); ++node_itr){
            np01_xy_node *node = *node_itr;
            AA_ALWAYS_ASSERT(NULL != node);
            const np01_float64 dx = node->get_x() - xy.first;
            const np01_float64 dy = node->get_y() - xy.second;
            const np01_float64 dsq = (dx*dx) + (dy*dy);
            if((NULL == near_node) || (dsq < near_dsq) ||
              ((dsq==near_dsq)&&(node->get_idx()<near_node->get_idx()))){
                near_node = node;
                near_dsq = dsq;
                }
            }
        search_vec->clear();

        if(NULL == near_node){
            switch(search_mode){
                case srch_md_circle:
                    search_r *= 2.0; /* double search radius */
                    if(search_r > max_search_r){
                        search_mode = srch_md_all; }
                    break;
                case srch_md_all:
                    search_mode=srch_md_done;
                    break;
                case srch_md_done:
                default:
                    AA_ALWAYS_ASSERT(false);
                    search_mode=srch_md_done;
                    break;
                }
            }
        else{
            search_mode=srch_md_done; }
        }

    /* free utility vector */
    if(NULL == utility_vec){
        delete search_vec; }
    }

AA_ALWAYS_ASSERT(near_dsq>=0.0);
if(NULL != d){ *d=sqrt(near_dsq); }

AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
return near_node;
}


/* find nearest node that is not connected to an edge.
xy: search center,  usable_node_lookup_table: must be sorted.  if non-NULL,
near node must also be found in this table.
utility_vec: used during search, then left empty. */
np01_xy_node *np01_xy_node_map::get_near_unconnected_node(const np01_xy& xy,
        const node_vec *usable_node_lookup_table, node_vec *utility_vec)const
{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);

np01_xy_node *near_node = NULL;
np01_float64 near_dsq = 0.0;
if(!m_node_vec.empty()){
    /* utility vector */
    node_vec *search_vec = (NULL==utility_vec) ?  new node_vec() : utility_vec;
    search_vec->clear();

    /* search radius */
    const np01_float64 w = static_cast<np01_float64>(m_loc_grid_dim.get_w());
    const np01_float64 h = static_cast<np01_float64>(m_loc_grid_dim.get_h());
    const np01_float64& sq_sz = m_loc_grid_dim.get_sq_size();
    const np01_float64 min_search_r = sq_sz / 2.0;
    const np01_float64 max_search_r = sq_sz * (w + h) / 2.0;
    np01_float64 search_r = min_search_r;

    /* search */
    enum search_mode_type {srch_md_circle, srch_md_all, srch_md_done};
    search_mode_type search_mode = ((min_search_r>0.0) && (max_search_r>0.0)) ?
        srch_md_circle : srch_md_all;
    while(srch_md_done != search_mode){
        AA_ALWAYS_ASSERT(NULL != search_vec);
        AUTO_ASSERT(search_vec->empty());
        AUTO_ASSERT(NULL == near_node);
        /* search in circle or get all nodes */
        switch(search_mode){
            case srch_md_circle:
                get_nodes_in_circle(xy,search_r,search_vec);
                break;
            case srch_md_all:
                get_all_nodes(search_vec);
                break;
            case srch_md_done:
            default:
                AA_ALWAYS_ASSERT(false);
                break;
            }

        /* find closest node that is not connected to an edge and is 
        also found in the usable lookup table */
        near_dsq = std::numeric_limits<np01_float64>::max();
        node_vec_citr node_itr = search_vec->begin();
        for(; node_itr != search_vec->end(); ++node_itr){
            np01_xy_node *node = *node_itr;
            AA_ALWAYS_ASSERT(NULL != node);
            const np01_float64 dx = node->get_x() - xy.first;
            const np01_float64 dy = node->get_y() - xy.second;
            const np01_float64 dsq = (dx*dx) + (dy*dy);
            if((NULL == near_node) || (dsq < near_dsq) ||
              ((dsq==near_dsq)&&(node->get_idx()<near_node->get_idx()))){
                if(node->get_edge() == NULL){
                    if(NULL == usable_node_lookup_table){
                        near_node = node;
                        near_dsq = dsq;
                        }
                    else{
                        const node_vec_citr usable_node_itr = std::lower_bound(
                            usable_node_lookup_table->begin(),
                            usable_node_lookup_table->end(), node);
                        if((usable_node_itr != usable_node_lookup_table->end())
                            && (*usable_node_itr == node)){
                            near_node = node;
                            near_dsq = dsq;
                            }
                        }
                    }
                }
            }
        search_vec->clear();

        if(NULL == near_node){
            switch(search_mode){
                case srch_md_circle:
                    search_r *= 2.0; /* double search radius */
                    if(search_r > max_search_r){
                        search_mode = srch_md_all; }
                    break;
                case srch_md_all:
                    search_mode=srch_md_done;
                    break;
                case srch_md_done:
                default:
                    AA_ALWAYS_ASSERT(false);
                    search_mode=srch_md_done;
                    break;
                }
            }
        else{
            search_mode=srch_md_done; }
        }

    /* free utility vector */
    if(NULL == utility_vec){
        delete search_vec; }
    }


AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
return near_node;
}

int np01_xy_node_map::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
const np01_xy_node *node, *next_node, *search_node;
const np01_uint16 w = m_loc_grid_dim.get_w();
const np01_uint16 h = m_loc_grid_dim.get_h();
const size_t loc_grid_size_check = 
    static_cast<size_t>(w) * static_cast<size_t>(h);
np01_uint16 i,j;
size_t idx, loc_grid_idx, grid_sq_node_count, loc_grid_total_node_count, found_count;

if((NULL != m_owner) && (this!=m_owner->get_node_map_a()) &&
    (this!=m_owner->get_node_map_b())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "this=%x  owner=%x  owner->node_map_a=%x  owner->node_map_a=%x\n", 
        this, m_owner, m_owner->get_node_map_a(), m_owner->get_node_map_b());
    }

err_cnt += m_loc_grid_dim.verify_data(err_msg, err_msg_capacity, err_msg_pos);

if (m_loc_grid.size() != loc_grid_size_check){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "m_loc_grid.size()=%i != (w=%i)*(h=%i)=%i\n", 
        m_loc_grid.size(), w, h, loc_grid_size_check);
    }
if( m_loc_grid.size() >= loc_grid_size_check){
    loc_grid_total_node_count=0;
    for (i=0; i<w; ++i) {
        for (j=0; j<h; ++j) {
            loc_grid_idx=(static_cast<size_t>(i)*static_cast<size_t>(h))+
                static_cast<size_t>(j);
            node=m_loc_grid.at(loc_grid_idx);
            grid_sq_node_count=0;
            while((grid_sq_node_count<=m_node_vec.size()) && (NULL != node)){
                ++grid_sq_node_count;
                ++loc_grid_total_node_count;
                if(this != node->get_owner()) {
                    ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "this=%x != (node=%x)->owner=%x\n", 
                        this, node, node->get_owner());
                    }
                if(node->get_idx() >= m_node_vec.size()) {
                    ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "(node=%x)->idx=%i >= m_node_vec.size=%i\n", 
                        node, node->get_idx(), m_node_vec.size());
                    }
                else if(m_node_vec.at(node->get_idx()) != node) {
                     ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "m_node_vec.at(node->idx=%i)=%x != node=%x\n", 
                        m_node_vec.at(node->get_idx()), node->get_idx(), node);
                    }
                if(node->get_loc_grid_i() != i) {
                    ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "(node=%x)->loc_grid_i=%i >= i=%i\n", 
                        node, node->get_loc_grid_i(), i);
                    }
                if(node->get_loc_grid_j() != j) {
                    ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "(node=%x)->loc_grid_j=%i >= j=%i\n", 
                        node, node->get_loc_grid_j(), j);
                    }

                /* advance */
                next_node = node->get_loc_grid_next();
                if(NULL != next_node){
                    if( node != next_node->get_loc_grid_prev() ) {
                        ++err_cnt;
                        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                            "node=%x != (next_node=%x)->loc_grid_prev=%x\n", 
                            node, next_node, next_node->get_loc_grid_prev() );
                        }
                    if(grid_sq_node_count >= m_node_vec.size()){
                        ++err_cnt;
                        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                            "grid_sq_node_count=%i >= m_node_vec.size=%i\n", 
                            grid_sq_node_count, m_node_vec.size() );
                        }
                    }
                node = next_node;
                }
            }
        }
    if(m_node_vec.size() != loc_grid_total_node_count){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "node count=%i != nodes in loc grid count=%i\n", 
            m_node_vec.size(), loc_grid_total_node_count );
        }
    }

for(idx=0; idx < m_node_vec.size(); ++idx){
    node = m_node_vec.at(idx);
    if( NULL == node ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            " m_node_vec.at(%i) == NULL\n", idx );
        }
    else{
        err_cnt += node->verify_data(err_msg, err_msg_capacity, err_msg_pos);
        if(this != node->get_owner()){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "idx=%i  this=%x != (node=%x)->owner=%x\n", 
                idx, this, node, node->get_owner() );
            }
        if(node->get_idx() != idx){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "(node=%x)->get_idx()=%i != idx=%i\n", 
                node, node->get_idx(), idx);
            }
        if((node->get_loc_grid_i() >= w) && 
           ((w>0) || node->get_loc_grid_i()> 0) ){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "(node=%x)->loc_grid_i=%i >= loc grid w=%i\n", 
                node, node->get_loc_grid_i(), w);
            }
        if((node->get_loc_grid_j() >= h) && 
           ((h>0) || node->get_loc_grid_j()> 0) ){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "(node=%x)->loc_grid_j=%i >= loc grid h=%i\n", 
                node, node->get_loc_grid_j(), h);
            }

        /* search for same node in locator grid at 
           get_loc_grid_i(), get_loc_grid_j()*/
        if( !m_loc_grid.empty()){
            loc_grid_idx=(static_cast<size_t>(node->get_loc_grid_i())*
              static_cast<size_t>(h))+ static_cast<size_t>(node->get_loc_grid_j());
            if((loc_grid_idx > 0) && ( loc_grid_idx >= m_loc_grid.size())) {
                ++err_cnt;
                np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                    "i=%i  j=%i  loc_grid_idx=%i >= loc_grid.size=%i\n", 
                    node->get_loc_grid_i(), node->get_loc_grid_j(),
                    loc_grid_idx, m_loc_grid.size());
                }
            else {
                search_node = m_loc_grid.at(loc_grid_idx);
                grid_sq_node_count=0;
                found_count=0;
                while((NULL != search_node) &&
                    (grid_sq_node_count < m_node_vec.size()) ){
                    if(search_node == node) {
                        ++grid_sq_node_count;
                        ++found_count; }
                    search_node = search_node->get_loc_grid_next();
                    }
                if(1 != found_count){
                    ++err_cnt;
                    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                        "i=%i  j=%i  loc_grid_idx=%i (node=%x) found count=%i\n", 
                        node->get_loc_grid_i(), node->get_loc_grid_j(),
                        loc_grid_idx, node, found_count);
                    }
                }
            }

        if(node->get_xy_group_type() != m_xy_group_type){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "(node=%x)->xy_group_type=%i != xy_group_type=%i\n", 
                node, node->get_xy_group_type(), m_xy_group_type);
            }
        }

    }
return err_cnt;
}

std::ostream& np01_xy_node_map::ostream_output(std::ostream& os) const{
os << "<xy_node_map>\n";
os << "<owner>" << std::hex << m_owner << std::dec << "</owner>\n";
os << "<xy_group_type>" << m_xy_group_type << "</xy_group_type>\n";
m_loc_grid_dim.ostream_output(os);
os << "<loc_grid><count>" << m_loc_grid.size() << "</count>\n";
np01_uint16 i = 0, j = 0;
for(node_vec_citr lg_node_itr=m_loc_grid.begin();
    lg_node_itr!=m_loc_grid.end(); ++lg_node_itr){
    const np01_xy_node *lg_node = *lg_node_itr;
    os << "</i=" << i << "></j=" << j << "><xy_node=" 
        << std::hex << lg_node << std::dec << ">";
    if( NULL != lg_node ){
        os << "<idx=" << (lg_node->get_idx())
            << "></i=" << (lg_node->get_loc_grid_i()) << "></j="
            << (lg_node->get_loc_grid_j()) << ">";
        }
    os << "</xy_node>\n";
    ++j;
    if (j >= m_loc_grid_dim.get_h()){ j = 0; ++i; }
    }
os << "</loc_grid>\n";
os << "<node_vec>\n";
for(node_vec_citr node_itr=m_node_vec.begin(); node_itr!=m_node_vec.end();
    ++node_itr){
    const np01_xy_node *node = *node_itr;
    if(NULL == node){
        os << "</xy_node=NULL>\n";
        }
    else{
        node->ostream_output(os);
        }
    }
os << "</node_vec>\n";
os << "<xy_node_map>\n";

return os;
}

void np01_xy_node_map::insert_node_in_loc_grid(np01_xy_node *n){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
np01_uint16 i,j;
size_t loc_grid_idx;
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(this == n->get_owner());
AUTO_ASSERT(NULL == n->get_loc_grid_prev());
AUTO_ASSERT(NULL == n->get_loc_grid_next());
i = m_loc_grid_dim.get_i(n->get_x());
j = m_loc_grid_dim.get_j(n->get_y());
loc_grid_idx = static_cast<size_t>(j) + 
    (static_cast<size_t>(i) * static_cast<size_t>(m_loc_grid_dim.get_h()));
AA_ALWAYS_ASSERT(loc_grid_idx < m_loc_grid.size());
np01_xy_node *nn = m_loc_grid[loc_grid_idx];
m_loc_grid[loc_grid_idx] = n;
n->set_loc_grid_next(nn);
n->set_loc_grid_i(i);
n->set_loc_grid_j(j);
if(NULL != nn){
    AUTO_ASSERT(NULL == nn->get_loc_grid_prev());
    AUTO_ASSERT(nn->get_loc_grid_i() == i);
    AUTO_ASSERT(nn->get_loc_grid_j() == j);
    nn->set_loc_grid_prev(n);
    }
AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_node_map::clear_loc_grid(){
AA_INCR_CALL_DEPTH();
node_vec_citr lg_itr;
np01_xy_node *node;
np01_xy_node *next_lg_node;
for(lg_itr=m_loc_grid.begin(); lg_itr!=m_loc_grid.end(); ++lg_itr){
    node=*lg_itr;
    while(NULL != node){
        next_lg_node=node->get_loc_grid_next();
        node->set_loc_grid_i(0);
        node->set_loc_grid_j(0);
        node->set_loc_grid_prev(NULL);
        node->set_loc_grid_next(NULL);
        node=next_lg_node;
        }
    }
m_loc_grid.clear();
m_loc_grid_dim.reset();
AA_DECR_CALL_DEPTH();
}

void np01_xy_node_map::free_nodes(){
AA_INCR_CALL_DEPTH();
node_vec_citr node_itr;
np01_xy_node *node;
for(node_itr=m_node_vec.begin(); node_itr!=m_node_vec.end(); ++node_itr){
    node=*node_itr;
    AA_ALWAYS_ASSERT(NULL != node);
    if(NULL != node->get_edge()){
        node->disconnect_edge(); }
    AUTO_ASSERT(NULL == node->get_loc_grid_prev());
    AUTO_ASSERT(NULL == node->get_loc_grid_next());
    node->set_owner(NULL);
    delete( node );
    }
m_node_vec.clear();
AA_DECR_CALL_DEPTH();
}


np01_xy_edge::np01_xy_edge(): m_idx(0), m_owner(NULL), m_node_a(NULL),
    m_node_b(NULL), m_loc_grid_head_node(NULL), m_ab_fwd_x(1.0),
    m_ab_fwd_y(0.0), m_fwd_dot_a(0.0), m_fwd_dot_b(0.0),  m_fwd_cross_ab(0.0)
{}

np01_xy_edge::~np01_xy_edge(){
/* if node_a not NULL, disconnect.
   if node_b not NULL, disconnect.
   if m_loc_grid_head_node not NULL, remove from locator grid. */
}

void np01_xy_edge::update_fwd_dot_cross_ab(){
AA_INCR_CALL_DEPTH();
if((NULL == m_node_a) || (NULL == m_node_b)){
    m_ab_fwd_x = 1.0;
    m_ab_fwd_y = 0.0;
    m_fwd_dot_a = (NULL != m_node_a) ? m_node_a->get_x() :
        (NULL != m_node_b) ? m_node_b->get_x() : 0.0;
    m_fwd_dot_b = m_fwd_dot_a;
    m_fwd_cross_ab = (NULL != m_node_a) ? m_node_a->get_y() :
        (NULL != m_node_b) ? m_node_b->get_y() : 0.0;
    }
else{
    np01_float64 dx = m_node_b->get_x() - m_node_a->get_x();
    np01_float64 dy = m_node_b->get_y() - m_node_a->get_y();
    if((fabs(dx) < 1e-144) && (fabs(dy) < 1e-144)){
        dx = (dx > 0.0) ? 1.0 : -1.0;
        dy = (dy > 0.0) ? 1.0 : -1.0;}
    const np01_float64 d = sqrt((dx*dx) + (dy*dy));
    const np01_float64 x_ab = (m_node_a->get_x() + m_node_b->get_x()) / 2.0;
    const np01_float64 y_ab = (m_node_a->get_y() + m_node_b->get_y()) / 2.0;
    m_ab_fwd_x = dx / d;
    m_ab_fwd_y = dy / d;
    m_fwd_dot_a = (m_ab_fwd_x*(m_node_a->get_x())) + 
                  (m_ab_fwd_y*(m_node_a->get_y()));
    m_fwd_dot_b = (m_ab_fwd_x*(m_node_b->get_x())) + 
                  (m_ab_fwd_y*(m_node_b->get_y()));
    m_fwd_cross_ab = (m_ab_fwd_x*y_ab) - (m_ab_fwd_y*x_ab);
    }
AA_DECR_CALL_DEPTH();
}

bool np01_xy_edge::is_in_bb(const np01_xy& xy_min, const np01_xy& xy_max)const{
AA_INCR_CALL_DEPTH();
bool result = false;

if((NULL != m_node_a) && (NULL != m_node_b)){
    const np01_float64& x_min = xy_min.first;
    const np01_float64& y_min = xy_min.second;
    const np01_float64& x_max = xy_max.first;
    const np01_float64& y_max = xy_max.second;
    AA_ALWAYS_ASSERT(x_min <= x_max);
    AA_ALWAYS_ASSERT(y_min <= y_max);
    const np01_float64& xa = m_node_a->get_x();
    const np01_float64& ya = m_node_a->get_y();
    const np01_float64& xb = m_node_b->get_x();
    const np01_float64& yb = m_node_b->get_y();
    const np01_float64 x_ctr = (x_min + x_max)/2.0;
    const np01_float64 y_ctr = (y_min + y_max)/2.0;
    const np01_float64 fwd_dot_ctr = (m_ab_fwd_x * x_ctr) + (m_ab_fwd_y * y_ctr);
    if(fwd_dot_ctr < m_fwd_dot_a){
        /* Node A is the closest point to bb center */
        if((x_min <= xa) && (xa <= x_max) && (y_min <= ya) && (ya <= y_max)){
            result = true;
            }
        }
    else if(fwd_dot_ctr > m_fwd_dot_b){
        /* Node B is the closest point to bb center */
        if((x_min <= xb) && (xb <= x_max) && (y_min <= yb) && (yb <= y_max)){
            result = true;
            }
        }
    else{
        const np01_float64 fwd_d = fwd_dot_ctr - m_fwd_dot_a;
        const np01_float64 x_near = xa + (fwd_d * m_ab_fwd_x);
        const np01_float64 y_near = ya + (fwd_d * m_ab_fwd_y);
        if((x_min <= x_near) && (x_near <= x_max) &&
            (y_min <= y_near) && (y_near <= y_max)){
            result = true;
            }
        else{
            /* nearest bb corner */
            const np01_float64 xc = (x_near < x_ctr) ? x_min : x_max;
            const np01_float64 yc = (y_near < y_ctr) ? y_min : y_max;
            /* near point to nearest corner */
            np01_float64 x_near_c, y_near_c;

            const np01_float64 fwd_dot_c = (m_ab_fwd_x*xc) + (m_ab_fwd_y*yc);
            if(fwd_dot_c < m_fwd_dot_a){
                /* Node A is the closest point to corner */
                x_near_c = xa;
                y_near_c = ya;
                }
            else if(fwd_dot_c > m_fwd_dot_b){
                /* Node B is the closest point to corner */
                x_near_c = xb;
                y_near_c = yb;
                }
            else{
                const np01_float64 fwd_d_c = fwd_dot_c - m_fwd_dot_a;
                x_near_c = xa + (fwd_d_c * m_ab_fwd_x);
                y_near_c = ya + (fwd_d_c * m_ab_fwd_y);
                }
            if((x_min <= x_near_c) && (x_near_c <= x_max) &&
                (y_min <= y_near_c) && (y_near_c <= y_max)){
                result = true;
                }
            }
        }
    }
AA_DECR_CALL_DEPTH();
return false;
}

bool np01_xy_edge::is_in_circle(const np01_xy& xy_ctr,
    const np01_float64& r) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(r >= 0.0);
bool result = false;
if((NULL != m_node_a) && (NULL != m_node_b)){
    const np01_float64& x_ctr = xy_ctr.first;
    const np01_float64& y_ctr = xy_ctr.second;
    const np01_float64 fwd_dot_ctr = (m_ab_fwd_x * x_ctr) + (m_ab_fwd_y * y_ctr);
    if(fwd_dot_ctr<m_fwd_dot_a){
        /* Node A is the closest point to circle center */
        const np01_float64 dx = (m_node_a->get_x()) - x_ctr;
        const np01_float64 dy = (m_node_a->get_y()) - y_ctr;
        const np01_float64 dsq = (dx*dx) + (dy*dy);
        const np01_float64 rsq = r*r;
        if(dsq <= rsq){
            result = true; }
        }
    else if(fwd_dot_ctr>m_fwd_dot_b){
        /* Node B is the closest point to circle center */
        const np01_float64 dx = (m_node_b->get_x()) - x_ctr;
        const np01_float64 dy = (m_node_b->get_y()) - y_ctr;
        const np01_float64 dsq = (dx*dx) + (dy*dy);
        const np01_float64 rsq = r*r;
        if(dsq <= rsq){
            result = true; }
        }
    else{
        /* circle center is closest to a point on edge between A and B. */
        const np01_float64 fwd_cross_ctr = (m_ab_fwd_x*y_ctr) - (m_ab_fwd_y*x_ctr);
        const np01_float64 d = fabs(m_fwd_cross_ab - fwd_cross_ctr);
        if(d <= r){
            result = true;}
        }
    }

AA_DECR_CALL_DEPTH();
return result;
}

/* check if quadrilateral (A0,A1,B0,B1) is convex 

         A1......B0
        .  *   /   .
        .    /*     .
       .   /     *   .
       . /          * .
      A0...............B1

*/
bool np01_xy_edge::intersects(const np01_xy_edge* other) const{
bool result = false;
if((NULL != m_node_a) && (NULL != m_node_b) && (NULL != other) &&
   (NULL != other->get_node_a()) && (NULL != other->get_node_b())){
   const np01_float64& xa0 = m_node_a->get_x();
   const np01_float64& ya0 = m_node_a->get_y();
   const np01_float64& xb0 = m_node_b->get_x();
   const np01_float64& yb0 = m_node_b->get_y();
   const np01_float64& xa1 = other->get_node_a()->get_x();
   const np01_float64& ya1 = other->get_node_a()->get_y();
   const np01_float64& xb1 = other->get_node_b()->get_x();
   const np01_float64& yb1 = other->get_node_b()->get_y();
   const np01_float64 dxa0a1 = xa1 - xa0;
   const np01_float64 dya0a1 = ya1 - ya0;
   const np01_float64 dxa1b0 = xb0 - xa1;
   const np01_float64 dya1b0 = yb0 - ya1;
   const np01_float64 dxb0b1 = xb1 - xb0;
   const np01_float64 dyb0b1 = yb1 - yb0;
   const np01_float64 dxb1a0 = xa0 - xb1;
   const np01_float64 dyb1a0 = ya0 - yb1;
   const np01_float64 cross_a0a1b0 = (dxa0a1*dya1b0)-(dya0a1*dxa1b0);
   const np01_float64 cross_a1b0b1 = (dxa1b0*dyb0b1)-(dya1b0*dxb0b1);
   const np01_float64 cross_b0b1a0 = (dxb0b1*dyb1a0)-(dyb0b1*dxb1a0);
   const np01_float64 cross_b1a0a1 = (dxb1a0*dya0a1)-(dyb1a0*dxa0a1);
   result = ( ( ( cross_a0a1b0 <= 0.0) &&
                ( cross_a1b0b1 <= 0.0) &&
                ( cross_b0b1a0 <= 0.0) &&
                ( cross_b1a0a1 <= 0.0) ) ||
              ( ( cross_a0a1b0 >= 0.0) &&
                ( cross_a1b0b1 >= 0.0) &&
                ( cross_b0b1a0 >= 0.0) &&
                ( cross_b1a0a1 >= 0.0) ) );
   }

return result;
}

int np01_xy_edge::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
size_t max_loc_grid_node_e_count=0;
size_t loc_grid_node_e_count=0;
const np01_xy_edge_loc_grid_node *loc_grid_node;
const np01_xy_edge_loc_grid_node *loc_grid_node_e_next;
size_t max_loc_grid_node_count=0;
size_t loc_grid_node_count=0;
const np01_xy_edge_loc_grid_node *loc_grid_head_node;
const np01_xy_edge_loc_grid_node *loc_grid_node_prev;
const np01_xy_edge_loc_grid_node *loc_grid_head_node_check;
const np01_xy_loc_grid_dim loc_grid_dim = 
    (NULL == m_owner) ? np01_xy_loc_grid_dim() : m_owner->get_loc_grid_dim();
const np01_xy_node_map *node_a_owner=(NULL==m_node_a) ? 
    NULL : m_node_a->get_owner();
const np01_xy_pair_map *node_a_owner_owner = (NULL==node_a_owner) ?
    NULL : node_a_owner->get_owner();
const np01_xy_node_map *node_b_owner=(NULL==m_node_b) ? 
    NULL : m_node_b->get_owner();
const np01_xy_pair_map *node_b_owner_owner = (NULL==node_b_owner) ?
    NULL : node_b_owner->get_owner();
const np01_xy_pair_map *owner_owner = (NULL==m_owner) ?
    NULL : m_owner->get_owner();
const np01_float64 min_max_dot_cross_err = NP01_DEFAULT_MAX_DOT_CROSS_ERR *
    ((NULL == m_owner) ? 1.0 : (m_owner->get_loc_grid_dim()).get_sq_size());
const np01_float64 max_dot_cross_err = 
    (min_max_dot_cross_err > NP01_DEFAULT_MAX_DOT_CROSS_ERR) ?
    min_max_dot_cross_err : NP01_DEFAULT_MAX_DOT_CROSS_ERR;
if( NULL != m_owner ) {
    max_loc_grid_node_e_count = static_cast<size_t>(loc_grid_dim.get_w()) 
        * static_cast<size_t>(loc_grid_dim.get_h());
    max_loc_grid_node_count = m_owner->get_edge_count();
    if( m_idx >= m_owner->get_edge_count() ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_idx(%i) >= m_owner->get_edge_count(%i)\n",
            m_idx, m_owner->get_edge_count()); }
    else if( this != m_owner->get_edge_by_idx(m_idx) ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != m_owner->get_edge_by_idx(%i)=%x\n",
            this, m_idx, m_owner->get_edge_by_idx(m_idx)); }

    if((NULL != m_node_a) &&( owner_owner != node_a_owner_owner)) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "((this=%x)->owner=%x)->owner=%x != "
            "((m_node_a=%x)->owner=%x)->owner=%x\n",this, m_owner, owner_owner,
            m_node_a, node_a_owner, node_a_owner_owner); }
    if((NULL != m_node_b) &&( owner_owner != node_b_owner_owner)) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "((this=%x)->owner=%x)->owner=%x != "
            "((m_node_b=%x)->owner=%x)->owner=%x\n",this, m_owner, owner_owner,
            m_node_b, node_a_owner, node_b_owner_owner); }
    }

if( (NULL != m_node_a) && (NULL != m_node_b)){
    const np01_float64 dx_ab = m_node_b->get_x() - m_node_a->get_x();
    const np01_float64 dy_ab = m_node_b->get_y() - m_node_a->get_y();
    const np01_float64 dsq_ab = (dx_ab*dx_ab) + (dy_ab*dy_ab);
    const np01_float64 d_ab = sqrt(dsq_ab);
    if(dsq_ab > 1e-50){
        const np01_float64 ab_fwd_x = dx_ab / d_ab;
        const np01_float64 ab_fwd_y = dy_ab / d_ab;
        if((fabs(ab_fwd_x - m_ab_fwd_x) > NP01_DEFAULT_MAX_FWD_ERR ) ||
           (fabs(ab_fwd_y - m_ab_fwd_y) > NP01_DEFAULT_MAX_FWD_ERR ) ){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i (m_fwd=%f,%f) != (fwd=%f,%f)"
                " | A=(%f,%f) B=(%f,%f)\n",
                this, m_idx, m_ab_fwd_x, m_ab_fwd_y, ab_fwd_x, ab_fwd_y,
                m_node_a->get_x(), m_node_a->get_y(),
                m_node_b->get_x(), m_node_b->get_y());
            }
        }

    const np01_float64 fwd_dot_a = (m_ab_fwd_x * (m_node_a->get_x())) +
        (m_ab_fwd_y * (m_node_a->get_y()));
    const np01_float64 fwd_cross_a = (m_ab_fwd_x * (m_node_a->get_y())) -
        (m_ab_fwd_y * (m_node_a->get_x()));
    if(fabs(fwd_dot_a - m_fwd_dot_a) > max_dot_cross_err ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge=%x:i%i fwd_dot_a=%f != (fwd=%f,%f) dot (A=%f,%f)=%f"
            " | (B=%f,%f)\n",
            this, m_idx, m_fwd_dot_a, m_ab_fwd_x, m_ab_fwd_y,
            m_node_a->get_x(), m_node_a->get_y(), fwd_dot_a,
            m_node_b->get_x(), m_node_b->get_y());
        }
    if(fabs(fwd_cross_a - m_fwd_cross_ab) > max_dot_cross_err ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge=%x:i%i fwd_cross_a=%f != (fwd=%f,%f) cross (A=%f,%f)=%f\n",
            this, m_idx, m_fwd_cross_ab, m_ab_fwd_x, m_ab_fwd_y,
            m_node_a->get_x(), m_node_a->get_y(), fwd_cross_a);
        }
    if(m_node_a->get_xy_group_type() != NP01_XY_GROUP_TYPE_A){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "(m_node_a=%x)->get_xy_group_type()=%i != A=%i\n",
            this, m_node_a, m_node_a->get_xy_group_type(),
            NP01_XY_GROUP_TYPE_A);
        }
    if(this != m_node_a->get_edge()) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != (node_a=%x)->edge=%x\n",
            this, m_node_a, m_node_a->get_edge());
        }


    const np01_float64 fwd_dot_b = (m_ab_fwd_x * (m_node_b->get_x())) +
        (m_ab_fwd_y * (m_node_b->get_y()));
    const np01_float64 fwd_cross_b = (m_ab_fwd_x * (m_node_b->get_y())) -
        (m_ab_fwd_y * (m_node_b->get_x()));
    if(fabs(fwd_dot_b - m_fwd_dot_b) > max_dot_cross_err ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge=%x:i%i fwd_dot_b=%f != (fwd=%f,%f) dot (B=%f,%f)=%f\n",
            this, m_idx, m_fwd_dot_b, m_ab_fwd_x, m_ab_fwd_y,
            m_node_b->get_x(), m_node_b->get_y(), fwd_dot_b);
        }
    if(fabs(fwd_cross_b - m_fwd_cross_ab) > max_dot_cross_err ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge=%x:i%i fwd_cross_b=%f != (fwd=%f,%f) cross (B=%f,%f)=%f\n",
            this, m_idx, m_fwd_cross_ab, m_ab_fwd_x, m_ab_fwd_y,
            m_node_b->get_x(), m_node_b->get_y(), fwd_cross_b);
        }
    if(m_node_b->get_xy_group_type() != NP01_XY_GROUP_TYPE_B){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "(m_node_b=%x)->get_xy_group_type()=%i != B=%i\n",
            this, m_node_b, m_node_b->get_xy_group_type(),
            NP01_XY_GROUP_TYPE_B);
        }
    if(this != m_node_b->get_edge()) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != (node_b=%x)->edge=%x\n",
            this, m_node_b, m_node_b->get_edge());
        }

    const np01_float64 ab_dx = m_node_b->get_x() - m_node_a->get_x();
    const np01_float64 ab_dy = m_node_b->get_y() - m_node_a->get_y();
    const np01_float64 ab_fwd_len_sq = 
        (m_ab_fwd_x*m_ab_fwd_x) + (m_ab_fwd_y*m_ab_fwd_y);
    const np01_float64 np01_ab_fwd_cross_ab =
        (m_ab_fwd_x * ab_dy)-(m_ab_fwd_y * ab_dx);
    if( fabs(ab_fwd_len_sq - 1.0) > 1e-8){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x:i%i  ab_fwd(%f,%f)len_sq=%f != 1.0\n",
            this, m_idx, m_ab_fwd_x, m_ab_fwd_y, ab_fwd_len_sq);
        }
    if( fabs(np01_ab_fwd_cross_ab) > 1e-8){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x:i%i  ab_fwd(%f,%f) cross A->B(%f,%f)=%f != 0.0\n",
            this, m_idx, m_ab_fwd_x, m_ab_fwd_y, ab_dx, ab_dy,
            np01_ab_fwd_cross_ab);
        }
    if(m_fwd_dot_a > m_fwd_dot_b){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x:i%i  fwd_dot_a=%f > fwd_dot_b=%f\n",
            this, m_idx, m_fwd_dot_a, m_fwd_dot_b);
        }
    }

if((NULL != m_loc_grid_head_node) && 
    (NULL != m_loc_grid_head_node->get_e_prev())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "this=%x:i%i  (m_loc_grid_head_node=%x)->get_e_prev=%x != NULL\n",
        this, m_idx, m_loc_grid_head_node, m_loc_grid_head_node->get_e_prev());
    }
loc_grid_node_e_count=0;
loc_grid_node = m_loc_grid_head_node;
while((NULL != loc_grid_node) &&
    (loc_grid_node_e_count < max_loc_grid_node_e_count)){
    ++loc_grid_node_e_count;
    err_cnt += loc_grid_node->verify_data(err_msg, err_msg_capacity,
        err_msg_pos);
    if(this != loc_grid_node->get_owner()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x != (loc_grid_node=%x)->owner=%x\n",
            this, loc_grid_node, loc_grid_node->get_owner());
        }
    if(NULL == loc_grid_node->get_edge_map()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x   (loc_grid_node=%x)->edge_map NULL\n",
            this, loc_grid_node);
        }
    else if(loc_grid_node->get_edge_map() != m_owner) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "this=%x   (loc_grid_node=%x)->edge_map=%x != owner=%x\n",
            this, loc_grid_node, loc_grid_node->get_edge_map(), m_owner);
        }
    else{
        /* max distance from edge to square center: sqrt(0.5) + tolerance */
        const np01_float64 max_sq_ctr_d = 0.8 * loc_grid_dim.get_sq_size(); 
        /* square center (x,y) */
        const np01_float64 sq_ctr_x =
            loc_grid_dim.get_sq_ctr_x(loc_grid_node->get_i());
        const np01_float64 sq_ctr_y =
            loc_grid_dim.get_sq_ctr_y(loc_grid_node->get_j());
        /* fwd_ab dot sq_ctr  */
        const np01_float64 fwd_dot_sq_ctr = (m_ab_fwd_x * sq_ctr_x) +
            (m_ab_fwd_y * sq_ctr_y);
        /* distance, measured parallel to A->B from square center to edge */
        const np01_float64 d_dot_sq_ctr =
          (m_fwd_dot_a>fwd_dot_sq_ctr)?(m_fwd_dot_a-fwd_dot_sq_ctr):/*before*/
          (fwd_dot_sq_ctr>m_fwd_dot_b)?(fwd_dot_sq_ctr-m_fwd_dot_b):/*after*/
          0.0; /* sq center lies within edge */
        /* fwd_ab cross sq_ctr  */
        const np01_float64 fwd_cross_sq_ctr = (m_ab_fwd_x * sq_ctr_y) -
            (m_ab_fwd_y * sq_ctr_x);
        /* distance, measured perpendicular to  A->B from sq center to edge */
        const np01_float64 d_cross_sq_ctr = fwd_cross_sq_ctr - m_fwd_cross_ab;

        if(d_dot_sq_ctr > max_sq_ctr_d){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i  aligned dist to loc grid sq ctr"
                "(x=%f,y=%f)  d_dot_sq_ctr=%f > max_sq_ctr_d=%f\n",
                this, m_idx, sq_ctr_x, sq_ctr_y, d_dot_sq_ctr, max_sq_ctr_d);
             }
        if(fabs(d_cross_sq_ctr) > max_sq_ctr_d){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i  perp dist to loc grid sq ctr"
                "(x=%f,y=%f)  fabs(d_cross_sq_ctr=%f) > max_sq_ctr_d=%f\n",
                this, m_idx, sq_ctr_x, sq_ctr_y, d_cross_sq_ctr, max_sq_ctr_d);
             }
        if( loc_grid_node->get_i() > loc_grid_dim.get_w()){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i  (loc_grid_node=%x)->i=%i > w=%i\n",
                this, m_idx, loc_grid_node, loc_grid_node->get_i(),
                loc_grid_dim.get_w());
             }
        if( loc_grid_node->get_j() > loc_grid_dim.get_h()){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i  (loc_grid_node=%x)->j=%i > h=%i\n",
                this, m_idx, loc_grid_node, loc_grid_node->get_j(),
                loc_grid_dim.get_h());
             }
 
        /* get head node from owner, search backward to find head node*/
        loc_grid_head_node = (NULL == m_owner) ? NULL :
            m_owner->get_loc_grid_head_node(loc_grid_node->get_i(), 
                loc_grid_node->get_j());
        loc_grid_node_count = 0;
        loc_grid_head_node_check = loc_grid_node;
        loc_grid_node_prev = loc_grid_node->get_prev();
        while( (NULL != loc_grid_node_prev) && 
            (loc_grid_node_count <= max_loc_grid_node_count)){
            ++loc_grid_node_count;
            loc_grid_head_node_check = loc_grid_node_prev;
            loc_grid_node_prev = loc_grid_node_prev->get_prev();
            }
        if(loc_grid_head_node_check->get_prev() != NULL){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i i=%i j=%i (loc_grid_head_node_check=%x)"
                "->get_prev()=%x != NULL\n",
                this, m_idx, loc_grid_node->get_i(), loc_grid_node->get_j(),
                loc_grid_head_node_check,loc_grid_head_node_check->get_prev());
            }
        if(loc_grid_head_node != loc_grid_head_node_check){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge=%x:i%i i=%i j=%i (loc_grid_head_node=%x)"
                "!= (loc_grid_head_node_check=%x)\n",
                this, m_idx, loc_grid_node->get_i(), loc_grid_node->get_j(),
                loc_grid_head_node, loc_grid_head_node_check);
            }
        }

    /* advance */
    loc_grid_node_e_next = loc_grid_node->get_e_next();
    if((NULL != loc_grid_node_e_next) &&
        (loc_grid_node_e_next->get_e_prev() != loc_grid_node)){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge=%x:i%i  ((loc_grid_node=%x)->e_next=%x)->e_prev=%x"
            " mismatch\n", this, m_idx, loc_grid_node, loc_grid_node_e_next,
            loc_grid_node_e_next->get_e_prev());
        }
    loc_grid_node = loc_grid_node_e_next;
    if( (NULL != loc_grid_node) && 
        (loc_grid_node_e_count >= max_loc_grid_node_e_count)){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
          "edge=%x:i%i loc_grid_node_e_count=%i>=max_loc_grid_node_e_count=%i\n",
          this, m_idx, loc_grid_node_e_count, max_loc_grid_node_e_count);
        }
    }

return err_cnt;
}

std::ostream& np01_xy_edge::ostream_output(std::ostream& os) const{
os << "<xy_edge=" << std::hex << this << std::dec << ">\n";
os << "<idx>" << m_idx << "</idx>\n";
os << "<owner>" << std::hex << m_owner << "</owner>\n";
os << "<node_a>" << m_node_a << "</node_a>\n";
os << "<node_b>" << m_node_b << "</node_b>\n";
os << "<loc_grid_head_node>"<<m_loc_grid_head_node<<"</loc_grid_head_node>\n";
if(NULL != m_loc_grid_head_node){
    const np01_xy_edge_loc_grid_node *lg_node = 
        m_loc_grid_head_node->get_e_next();
    size_t lg_count = 1;
    size_t lg_count_max = (NULL == m_owner) ? 32 : 
        ((m_owner->get_loc_grid_dim().get_w()) * 
        (m_owner->get_loc_grid_dim().get_h()));
    while((NULL != lg_node) && (lg_count < lg_count_max)){
        os << "  <lg_node>"<<lg_node<<"</lg_node>\n";
        lg_node = lg_node->get_e_next();
        ++lg_count;
        }    
    }
os << std::dec << "<ab_fwd_x>" << m_ab_fwd_x << "</ab_fwd_x>"
               << "<ab_fwd_y>" << m_ab_fwd_y << "</ab_fwd_y>\n";
os << "<fwd_dot_a>" << m_fwd_dot_a << "</fwd_dot_a>"
   << "<fwd_dot_b>" << m_fwd_dot_b << "</fwd_dot_b>\n";
os << "<fwd_cross_ab>" << m_fwd_cross_ab << "</fwd_cross_ab>\n";
os << "</xy_edge>\n";
return os;
}
 
int np01_xy_edge_loc_grid_node::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
size_t loc_grid_node_e_count;
size_t max_loc_grid_node_e_count=0;
size_t loc_grid_node_count;
size_t max_loc_grid_node_count=0;
const np01_xy_edge_loc_grid_node *loc_grid_node;
size_t found_count;
np01_xy_loc_grid_dim loc_grid_dim;

if( NULL != m_edge_map ) {
    loc_grid_dim = m_edge_map->get_loc_grid_dim();
    max_loc_grid_node_e_count = static_cast<size_t>(loc_grid_dim.get_w()) 
        * static_cast<size_t>(loc_grid_dim.get_h());
    max_loc_grid_node_count = m_edge_map->get_edge_count();
    }

if(NULL != m_owner){
    if(m_owner->get_owner() != m_edge_map){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node:%x (owner=%x)->owner=%x != edge_map=%x\n",
            this, m_owner, m_owner->get_owner(), m_edge_map);
        }

    /* m_owner->loc_grid_head->e_next->e_next ... this*/
    loc_grid_node = m_owner->get_loc_grid_head_node();
    loc_grid_node_e_count = 0;
    found_count = 0;
    while((NULL != loc_grid_node) &&
        (loc_grid_node_e_count <= max_loc_grid_node_e_count)){
        if(this == loc_grid_node){
            ++found_count; }
        ++loc_grid_node_e_count;
        loc_grid_node = loc_grid_node->get_e_next();
        }
    if(1 != found_count){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node:%x m_owner->loc_grid_head->"
            "e_next found_count=%i\n", this, found_count);
        }
    }

if(NULL != m_edge_map && NULL != m_owner){
    if(m_i >= loc_grid_dim.get_w()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node:%x  i=%i >= w=%i\n", this, m_i,
            loc_grid_dim.get_w());
        }
    if(m_j >= loc_grid_dim.get_h()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node:%x  j=%i >= h=%i\n", this, m_j,
            loc_grid_dim.get_h());
        }

    /* m_edge_map->get_loc_grid_head(m_i,m_j)->next->next-> ... this */
    loc_grid_node = m_edge_map->get_loc_grid_head_node(m_i, m_j);
    loc_grid_node_count = 0;
    found_count = 0;
    while((NULL != loc_grid_node) &&
        (loc_grid_node_count <= max_loc_grid_node_count)){
        if(this == loc_grid_node){
            ++found_count; }
        ++loc_grid_node_count;
        loc_grid_node = loc_grid_node->get_next();
        }
    if(1 != found_count){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node:%x edge_map->loc_grid_head(i=%i,j=%i)->next"
            "->next->... found_count=%i\n", this, m_i, m_j, found_count);
        }
    }

if(NULL != m_edge_map && NULL == m_owner){
    /* should be found on m_edge_map->m_loc_grid_node_free_chain */
    }
    
if(NULL != m_prev){
    if(this != m_prev->m_next){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node: this=%x != (prev=%x)->next=%x\n",
            this, m_prev, m_prev->m_next);
        }
    if((m_prev->m_owner == m_owner) && (NULL != m_owner)){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (prev=%x)->owner == owner=%x\n",
            this, m_prev, m_owner);
        }
    if(m_prev->m_edge_map != m_edge_map ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (prev=%x)->edge_map=%x != edge_map=%x\n",
            this, m_prev, m_prev->m_edge_map, m_edge_map);
        }
    }
             
if(NULL != m_next){
    if(this != m_next->m_prev){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node: this=%x != (next=%x)->prev=%x\n",
            this, m_next, m_next->m_prev);
        }
    if((m_next->m_owner == m_owner) && (NULL != m_owner) ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (next=%x)->owner == owner=%x\n",
            this, m_next, m_owner);
        }
    if(m_next->m_edge_map != m_edge_map ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (next=%x)->edge_map=%x != edge_map=%x\n",
            this, m_next, m_next->m_edge_map, m_edge_map);
        }
    }
    
if(NULL != m_e_prev){
    if(this != m_e_prev->m_e_next){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node: this=%x != (e_prev=%x)->e_next=%x\n",
            this, m_e_prev, m_e_prev->m_e_next);
        }
    if(m_e_prev->m_owner != m_owner ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (e_prev=%x)->owner=%x != owner=%x\n",
            this, m_e_prev, m_e_prev->m_owner, m_owner);
        }
    if(m_e_prev->m_edge_map != m_edge_map ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (e_prev=%x)->edge_map=%x != edge_map=%x\n",
            this, m_e_prev, m_e_prev->m_edge_map, m_edge_map);
        }
    }

if(NULL != m_e_next){
    if(this != m_e_next->m_e_prev){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node: this=%x != (e_next=%x)->e_prev=%x\n",
            this, m_e_next, m_e_next->m_e_prev);
        }
    if(m_e_next->m_owner != m_owner ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (e_next=%x)->owner=%x != owner=%x\n",
            this, m_e_next, m_e_next->m_owner, m_owner);
        }
    if(m_e_next->m_edge_map != m_edge_map ){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_loc_grid_node=%x: (e_next=%x)->edge_map=%x != edge_map=%x\n",
            this, m_e_next, m_e_next->m_edge_map, m_edge_map);
        }
    }

return err_cnt;
}

std::ostream& np01_xy_edge_loc_grid_node::ostream_output(std::ostream& os) const{
os << "<xy_edge_loc_grid_node=" << std::hex << this << ">\n";
os << "<owner>" << m_owner << "</owner>\n";
os << "<edge_map>" << m_edge_map << "</edge_map>\n";
os << std::dec << "<i>" << m_i << "</i><j>" << m_j << "</j>\n";
os << std::hex << "<prev>" << m_prev << "</prev>"
               << "<next>" << m_next << "</next>\n";
os << "<e_prev>" << m_e_prev << "</e_prev>"
   << "<e_next>" << m_e_next << "</e_next>\n";
os << std::dec << "</xy_edge_loc_grid_node>\n";
return os;
}

np01_xy_edge_map::np01_xy_edge_map(): m_owner(NULL), m_edge_vec(),
    m_loc_grid_dim(), m_loc_grid(), m_loc_grid_node_free_chain(NULL){}

np01_xy_edge_map::~np01_xy_edge_map(){
AA_INCR_CALL_DEPTH();
clear_loc_grid();
free_all_edges();
free_free_chain();
AA_DECR_CALL_DEPTH();
}

void np01_xy_edge_map::init_loc_grid( np01_xy_loc_grid_dim& d ){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
size_t loc_grid_sz;
if(!m_loc_grid.empty()){
    clear_loc_grid();}
AUTO_ASSERT(m_loc_grid.empty());
m_loc_grid_dim = d;
loc_grid_sz = static_cast<size_t>(d.get_w()) * static_cast<size_t>(d.get_h());
m_loc_grid.resize(loc_grid_sz, NULL);
if(!m_edge_vec.empty()){
    edge_vec_citr edge_itr=m_edge_vec.begin();
    for(;edge_itr!= m_edge_vec.end(); ++edge_itr){
        np01_xy_edge *edge = *edge_itr;
        AA_ALWAYS_ASSERT(NULL != edge);
        insert_edge_in_loc_grid(edge);
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

np01_xy_edge *np01_xy_edge_map::create_edge( np01_xy_node *node_a,
  np01_xy_node *node_b ) {
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(NULL != node_a);
AA_ALWAYS_ASSERT(NULL != node_b);
AUTO_ASSERT(node_a != node_b);
AUTO_ASSERT(node_a->get_xy_group_type() == NP01_XY_GROUP_TYPE_A);
AUTO_ASSERT(node_b->get_xy_group_type() == NP01_XY_GROUP_TYPE_B);
AUTO_ASSERT(node_a->get_edge() == NULL);
AUTO_ASSERT(node_b->get_edge() == NULL);
np01_xy_edge *edge = new np01_xy_edge();
edge->set_idx(static_cast<np01_uint32>(m_edge_vec.size()));
m_edge_vec.push_back(edge);
edge->set_owner(this);
edge->set_node_a(node_a);
node_a->set_edge(edge);
edge->set_node_b(node_b);
node_b->set_edge(edge);
edge->update_fwd_dot_cross_ab();
AA_ALWAYS_ASSERT(NULL == edge->get_loc_grid_head_node());
insert_edge_in_loc_grid(edge);
AUTO_ASSERT(NULL != edge->get_loc_grid_head_node());
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
return edge;
}

/* reassign new nodes A and B to edge */
void np01_xy_edge_map::update_edge_ab(np01_xy_edge *edge, np01_xy_node *node_a,
    np01_xy_node *node_b){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(NULL != edge);
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT(this == edge->get_owner());
AA_ALWAYS_ASSERT(NULL != node_a);
AA_ALWAYS_ASSERT(NULL != node_b);
remove_edge_from_loc_grid(edge);
np01_xy_node *na = edge->get_node_a();
if(na != node_a){
    if(NULL != na){
        na->disconnect_edge();}
    edge->set_node_a(node_a);
    node_a->set_edge(edge);
    }
np01_xy_node *nb = edge->get_node_b();
if(nb != node_b){
    if(NULL != nb){
        nb->disconnect_edge();}
    edge->set_node_b(node_b);
    node_b->set_edge(edge);
    }
edge->update_fwd_dot_cross_ab();
insert_edge_in_loc_grid(edge);
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}


/* prerequisites: 
  nodes A&B linked to edge.
  edge fwd,dot,cross AB updated
  m_loc_grid initialized


   +-------+-------+-------+-------+-------+-------+-------+
   |       |       |       |       |       | B     |       |
   |       |       |       |       |       | *     |       |
   |       |       |       |       |      *|       |       |
   +-------+-------+-------+-------+---*---+-------+-------+
   |       |       |       |       |*      |       |       |
   |       |       |       |     * |       |       |       |
   |       |       |       |  *    |       |       |       |
   +-------+-------+-------*-------+-------+-------+-------+
   |       |       |    *  |       |       |       |       |
   |       |       | *     |   +   |       |       |       |
   |       |      *|       |       |       |       |       |
   +-------+---*---+-------+-------+-------+-------+-------+
   |    A  |*      |       |       |       |       |       |
   |    *  |       |       |       |       |       |       |
   |       |       |       |       |       |       |       |
   +-------+-------+-------+-------+-------+-------+-------+
   |       |       |       |       |       |       |       |
   |       |       |       |       |       |       |       |
   |       |       |       |       |       |       |       |
   +-------+-------+-------+-------+-------+-------+-------+
   distance from edge AB to square center, when edge touches
   square boundary = sq_size * max( fwd cross (-0.5,-0.5),
                                    fwd cross (-0.5,+0.5),
                                    fwd cross (+0.5,-0.5),
                                    fwd cross (+0.5,+0.5))

*/
void np01_xy_edge_map::insert_edge_in_loc_grid(np01_xy_edge *edge){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != edge);
AA_ALWAYS_ASSERT(this == edge->get_owner());
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT(!m_loc_grid.empty());
AUTO_ASSERT(m_loc_grid_dim.get_w() > 0);
AUTO_ASSERT(m_loc_grid_dim.get_h() > 0);
AUTO_ASSERT(m_loc_grid_dim.get_sq_size() > 0.0);
if(NULL!=edge->get_loc_grid_head_node()){
    remove_edge_from_loc_grid(edge);}
AUTO_ASSERT(NULL==edge->get_loc_grid_head_node());
const np01_float64 fwd_x = edge->get_ab_fwd_x();
const np01_float64 fwd_y = edge->get_ab_fwd_y();
const np01_float64 fwd_xy_abs_sum = fabs(fwd_x + fwd_y);
const np01_float64 fwd_xy_abs_diff = fabs(fwd_x - fwd_y);
const np01_float64 fwd_xy_abs_max_sum_diff =
    (fwd_xy_abs_sum > fwd_xy_abs_diff) ? fwd_xy_abs_sum : fwd_xy_abs_diff;
const np01_float64 d_threshold = 0.5625 * (m_loc_grid_dim.get_sq_size())
    * fwd_xy_abs_max_sum_diff;
const np01_float64 cross_min_threshold= edge->get_fwd_cross_ab() - d_threshold;
const np01_float64 cross_max_threshold= edge->get_fwd_cross_ab() + d_threshold;
const np01_xy_node * const node_a = edge->get_node_a();
const np01_xy_node * const node_b = edge->get_node_b();
AA_ALWAYS_ASSERT(NULL != node_a);
AA_ALWAYS_ASSERT(NULL != node_b);
AUTO_ASSERT(node_a != node_b);
AUTO_ASSERT(node_a->get_xy_group_type() == NP01_XY_GROUP_TYPE_A);
AUTO_ASSERT(node_b->get_xy_group_type() == NP01_XY_GROUP_TYPE_B);
AUTO_ASSERT(node_a->get_edge() == edge);
AUTO_ASSERT(node_b->get_edge() == edge);
const np01_float64& xa = node_a->get_x();
const np01_float64& ya = node_a->get_y();
const np01_float64& xb = node_b->get_x();
const np01_float64& yb = node_b->get_y();
const np01_float64 x_min = (xa < xb) ? xa : xb;
const np01_float64 y_min = (ya < yb) ? ya : yb;
const np01_float64 x_max = (xa > xb) ? xa : xb;
const np01_float64 y_max = (ya > yb) ? ya : yb;
const np01_float64 grid_sq_sz_tol = m_loc_grid_dim.get_sq_size() * 0.0625;
const np01_uint16 i_min = m_loc_grid_dim.get_i(x_min - grid_sq_sz_tol);
const np01_uint16 j_min = m_loc_grid_dim.get_j(y_min - grid_sq_sz_tol);
const np01_uint16 i_max = m_loc_grid_dim.get_i(x_max + grid_sq_sz_tol);
const np01_uint16 j_max = m_loc_grid_dim.get_j(y_max + grid_sq_sz_tol);
np01_uint16 i,j;
size_t loc_grid_idx;
np01_float64 sq_ctr_x, sq_ctr_y, fwd_cross_sq_ctr;
for(i = i_min; i <= i_max; ++i){
    sq_ctr_x = m_loc_grid_dim.get_sq_ctr_x(i);
    for(j = j_min; j <= j_max; ++j){
        bool use_grid_sq = false;
        if((i_min==i_max) || (j_min==j_max)){
            use_grid_sq = true; }
        else{
            sq_ctr_y = m_loc_grid_dim.get_sq_ctr_y(j);
            fwd_cross_sq_ctr = (fwd_x * sq_ctr_y) - (fwd_y * sq_ctr_x);
            if( (cross_min_threshold <= fwd_cross_sq_ctr) &&
                (fwd_cross_sq_ctr <= cross_max_threshold) ){
                use_grid_sq = true;
                }
            }
        if( use_grid_sq ){
            loc_grid_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
                static_cast<size_t>(m_loc_grid_dim.get_h()));
            np01_xy_edge_loc_grid_node *lg_next = m_loc_grid[loc_grid_idx];
            np01_xy_edge_loc_grid_node *lg_node = alloc_loc_grid_node();
            lg_node->set_owner(edge);
            AUTO_ASSERT(lg_node->get_edge_map() == this);
            lg_node->set_i(i);
            lg_node->set_j(j);
            AUTO_ASSERT(lg_node->get_prev() == NULL);
            lg_node->set_next(lg_next);
            if(NULL != lg_next){
                lg_next->set_prev(lg_node); }
            m_loc_grid[loc_grid_idx] = lg_node;
            AUTO_ASSERT(lg_node->get_e_prev() == NULL);
            np01_xy_edge_loc_grid_node *lg_e_next =
                edge->get_loc_grid_head_node();
            if(NULL != lg_e_next){
                lg_e_next->set_e_prev(lg_node); }
            lg_node->set_e_next(lg_e_next);
            edge->set_loc_grid_head_node(lg_node);
            }
        }
    }
AUTO_ASSERT(NULL != edge->get_loc_grid_head_node());
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_edge_map::remove_edge_from_loc_grid(np01_xy_edge *edge){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(NULL != edge);
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AUTO_ASSERT(this == edge->get_owner());
np01_xy_edge_loc_grid_node *lg_node = edge->get_loc_grid_head_node();
while(NULL != lg_node){
    np01_xy_edge_loc_grid_node *lg_prev = lg_node->get_prev();
    np01_xy_edge_loc_grid_node *lg_next = lg_node->get_next();
    np01_xy_edge_loc_grid_node *lg_e_next = lg_node->get_e_next();
    if(NULL == lg_prev){
        const size_t lg_idx = static_cast<size_t>(lg_node->get_j()) +
            (static_cast<size_t>(lg_node->get_i())* 
             static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid.size());
        if(lg_idx < m_loc_grid.size()){
            AUTO_ASSERT(m_loc_grid.at(lg_idx) == lg_node);
            m_loc_grid[lg_idx] = lg_next;
            }
        }
    else{
        lg_prev->set_next(lg_next);
        }
    if(NULL != lg_next){
        lg_next->set_prev(lg_prev);
        }
    free_loc_grid_node(lg_node);
    lg_node=lg_e_next;
    }
edge->set_loc_grid_head_node(NULL);
AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_edge_map::get_edges_near_edge(const np01_xy_edge *e,
    edge_vec *edges, idx_edge_ptr_pair_vec *utility_vec ) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
if((NULL != e) && (NULL != edges)){
    AUTO_ASSERT(this == e->get_owner());
    const np01_xy_edge_loc_grid_node *lg_node = e->get_loc_grid_head_node();
    if(NULL == lg_node){
        /* get all edges except e */
        edge_vec_citr edge_itr = m_edge_vec.begin();
        for(; edge_itr != m_edge_vec.end(); ++edge_itr){
            np01_xy_edge *edge = *edge_itr;
            AUTO_ASSERT(NULL != edge);
            if(edge != e){
                edges->push_back(edge); }
            }
        }
    else{
        /* utility vector */
        idx_edge_ptr_pair_vec *search_vec = (NULL == utility_vec) ? 
            new idx_edge_ptr_pair_vec() : utility_vec;
        size_t search_vec_original_size = search_vec->size();

        /* search all grid squares occupied by e */
        while(NULL != lg_node){
            AUTO_ASSERT(0 == lg_node->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));

            /* search up in same grid square */
            const np01_xy_edge_loc_grid_node *lg_sq_node = lg_node->get_prev();
            while(NULL != lg_sq_node){
                np01_xy_edge *edge = lg_sq_node->get_owner();
                AUTO_ASSERT(edge != e);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                AA_XDBG_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()),
                    CF01_AA_DEBUG_LEVEL_2);
                lg_sq_node = lg_sq_node->get_prev();
                }

            /* search down in same grid square */
            lg_sq_node = lg_node->get_next();
            while(NULL != lg_sq_node){
                np01_xy_edge *edge = lg_sq_node->get_owner();
                AUTO_ASSERT(edge != e);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                AA_XDBG_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()),
                    CF01_AA_DEBUG_LEVEL_2);
                lg_sq_node = lg_sq_node->get_next();
                }

            lg_node = lg_node->get_e_next();
            }
    
        /* sort search vec, remove duplicates */
        idx_edge_ptr_pair_vec_itr sort_itr =
            search_vec->begin() + search_vec_original_size;
        std::sort(sort_itr, search_vec->end());
        idx_edge_ptr_pair_vec_itr unique_end_itr =
            std::unique(sort_itr, search_vec->end());
    
        /* copy nearby edges to output vector */
        idx_edge_ptr_pair_vec_itr unique_itr = sort_itr;
        for(; unique_itr != unique_end_itr; ++unique_itr){
            np01_xy_edge *edge = unique_itr->second;
            edges->push_back(edge);
            }
    
        /* restore utility vector */
        if(NULL == utility_vec){
            delete search_vec; }
        else{
            utility_vec->resize(search_vec_original_size); }
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

/* get all edges within bounding box */
void np01_xy_edge_map::get_edges_in_bb(const np01_xy& xy_min,
    const np01_xy& xy_max, edge_vec *edges,
    idx_edge_ptr_pair_vec *utility_vec) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(xy_min.first <= xy_max.first);
AA_ALWAYS_ASSERT(xy_min.second <= xy_max.second);
AA_ALWAYS_ASSERT(NULL != edges);

/* utility vector */
idx_edge_ptr_pair_vec *search_vec = (NULL == utility_vec) ? 
    new idx_edge_ptr_pair_vec() : utility_vec;
size_t search_vec_original_size = search_vec->size();

/* search locator grid, fill search_vec */
size_t lg_idx;
np01_uint16 i, j;
np01_uint16_pair ij_min, ij_max;
m_loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);
for(i = ij_min.first; i <= ij_max.first; ++i){
    for(j = ij_min.second; j <= ij_max.second; ++j){
        lg_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
            static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid.size());
        if(lg_idx < m_loc_grid.size()){
            np01_xy_edge_loc_grid_node *lg_node = m_loc_grid[lg_idx];
            while(NULL != lg_node){
                np01_xy_edge *edge = lg_node->get_owner();
                AA_ALWAYS_ASSERT(NULL != edge);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                lg_node = lg_node->get_next();
                }
            }
        }
    }

/* sort search vec, remove duplicates */
idx_edge_ptr_pair_vec_itr sort_itr =
    search_vec->begin() + search_vec_original_size;
std::sort(sort_itr, search_vec->end());
idx_edge_ptr_pair_vec_itr unique_end_itr =
    std::unique(sort_itr, search_vec->end());

/* copy nearby edges to output vector */
idx_edge_ptr_pair_vec_itr unique_itr = sort_itr;
for(; unique_itr != unique_end_itr; ++unique_itr){
    np01_xy_edge *edge = unique_itr->second;
    if(edge->is_in_bb(xy_min, xy_max)){
        edges->push_back(edge);
        }
    }

/* restore utility vector */
if(NULL == utility_vec){
    delete search_vec; }
else{
    utility_vec->resize(search_vec_original_size); }

AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

/* get all edges within circle 
edges output container
utility_vec, if non-NULL, a helper vector that will be added to,
  then left the same size, essentially unchanged.  Use should 
  speed up algorithm by avoiding excessive memory allocation
  and deallocation.

*/
void np01_xy_edge_map::get_edges_in_circle(const np01_xy& xy_ctr,
    const np01_float64& r, edge_vec *edges,
    idx_edge_ptr_pair_vec *utility_vec) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(r >= 0.0);
AA_ALWAYS_ASSERT(NULL != edges);
const np01_xy xy_min(xy_ctr.first - r, xy_ctr.second - r);
const np01_xy xy_max(xy_ctr.first + r, xy_ctr.second + r);
size_t lg_idx;
np01_uint16 i, j;
np01_uint16_pair ij_min, ij_max;
m_loc_grid_dim.get_bb_indices(xy_min, xy_max, &ij_min, &ij_max);

/* utility vector */
idx_edge_ptr_pair_vec *search_vec = (NULL == utility_vec) ? 
    new idx_edge_ptr_pair_vec() : utility_vec;
size_t search_vec_original_size = search_vec->size();

/* search locator grid, fill search_vec */
for(i = ij_min.first; i <= ij_max.first; ++i){
    for(j = ij_min.second; j <= ij_max.second; ++j){
        lg_idx = static_cast<size_t>(j) + (static_cast<size_t>(i) *
            static_cast<size_t>(m_loc_grid_dim.get_h()));
        AA_ALWAYS_ASSERT(lg_idx < m_loc_grid.size());
        if(lg_idx < m_loc_grid.size()){
            np01_xy_edge_loc_grid_node *lg_node = m_loc_grid[lg_idx];
            while(NULL != lg_node){
                np01_xy_edge *edge = lg_node->get_owner();
                AA_ALWAYS_ASSERT(NULL != edge);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                lg_node = lg_node->get_next();
                }
            }
        }
    }

/* sort search vec, remove duplicates */
idx_edge_ptr_pair_vec_itr sort_itr =
    search_vec->begin() + search_vec_original_size;
std::sort(sort_itr, search_vec->end());
idx_edge_ptr_pair_vec_itr unique_end_itr =
    std::unique(sort_itr, search_vec->end());

/* copy nearby edges to output vector */
idx_edge_ptr_pair_vec_itr unique_itr = sort_itr;
for(; unique_itr != unique_end_itr; ++unique_itr){
    np01_xy_edge *edge = unique_itr->second;
    if(edge->is_in_circle(xy_ctr, r)){
        edges->push_back(edge);
        }
    }

/* restore utility vector */
if(NULL == utility_vec){
    delete search_vec; }
else{
    utility_vec->resize(search_vec_original_size); }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

/* get all edges that intersect edge e */
void np01_xy_edge_map::get_intersecting_edges(const np01_xy_edge *e,
    edge_vec *edges, idx_edge_ptr_pair_vec *utility_vec ) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
if((NULL != e) && (NULL != edges)){
    AUTO_ASSERT(this == e->get_owner());
    const np01_xy_edge_loc_grid_node *lg_node = e->get_loc_grid_head_node();
    if(NULL == lg_node){
        /* get all edges except e */
        edge_vec_citr edge_itr = m_edge_vec.begin();
        for(; edge_itr != m_edge_vec.end(); ++edge_itr){
            np01_xy_edge *edge = *edge_itr;
            AUTO_ASSERT(NULL != edge);
            if((edge != e) && (edge->intersects(e))){
                edges->push_back(edge); }
            }
        }
    else{
        /* utility vector */
        idx_edge_ptr_pair_vec *search_vec = (NULL == utility_vec) ? 
            new idx_edge_ptr_pair_vec() : utility_vec;
        size_t search_vec_original_size = search_vec->size();

        /* search all grid squares occupied by e */
        while(NULL != lg_node){
            AUTO_ASSERT(0 == lg_node->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));

            /* search up in same grid square */
            const np01_xy_edge_loc_grid_node *lg_sq_node = lg_node->get_prev();
            while(NULL != lg_sq_node){
                np01_xy_edge *edge = lg_sq_node->get_owner();
                AUTO_ASSERT(edge != e);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                AUTO_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
                lg_sq_node = lg_sq_node->get_prev();
                }

            /* search down in same grid square */
            lg_sq_node = lg_node->get_next();
            while(NULL != lg_sq_node){
                np01_xy_edge *edge = lg_sq_node->get_owner();
                AUTO_ASSERT(edge != e);
                search_vec->push_back(idx_edge_ptr_pair(edge->get_idx(),edge));
                AUTO_ASSERT(0 == lg_sq_node->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
                lg_sq_node = lg_sq_node->get_next();
                }

            lg_node = lg_node->get_e_next();
            }
    
        /* sort search vec, remove duplicates */
        idx_edge_ptr_pair_vec_itr sort_itr =
            search_vec->begin() + search_vec_original_size;
        std::sort(sort_itr, search_vec->end());
        idx_edge_ptr_pair_vec_itr unique_end_itr =
            std::unique(sort_itr, search_vec->end());
    
        /* copy nearby edges to output vector */
        idx_edge_ptr_pair_vec_itr unique_itr = sort_itr;
        for(; unique_itr != unique_end_itr; ++unique_itr){
            np01_xy_edge *edge = unique_itr->second;
            AUTO_ASSERT(edge != e);
            if(edge->intersects(e)){
                edges->push_back(edge);
                }
            }
    
        /* restore utility vector */
        if(NULL == utility_vec){
            delete search_vec; }
        else{
            utility_vec->resize(search_vec_original_size); }
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}




void np01_xy_edge_map::get_edge_crossover_count(
    np01_uint32_vec *edge_crossover_count_vec,
    np01_uint32 *total_edge_crossover_count) const{
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
if(NULL != edge_crossover_count_vec){
    edge_crossover_count_vec->clear();
    edge_crossover_count_vec->resize(m_edge_vec.size(), 0);
    }
np01_uint32 total_x_count = 0;
np01_uint32 edge_idx_i, edge_idx_j;
edge_vec local_edges;
idx_edge_ptr_pair_vec utility_vec;
edge_idx_i=0;
for(edge_vec_citr edge_itr_i = m_edge_vec.begin(); edge_itr_i != m_edge_vec.end();
    ++edge_itr_i, ++edge_idx_i){
    const np01_xy_edge *edge_i = *edge_itr_i;
    AA_ALWAYS_ASSERT(NULL != edge_i);
    AA_ALWAYS_ASSERT(edge_i->get_idx() == edge_idx_i);
    get_edges_near_edge(edge_i, &local_edges, &utility_vec);
    for(edge_vec_citr edge_itr_j = local_edges.begin();
        edge_itr_j != local_edges.end(); ++edge_itr_j ){
        const np01_xy_edge *edge_j = *edge_itr_j;
        AA_ALWAYS_ASSERT(NULL != edge_j);
        AUTO_ASSERT(edge_i != edge_j)
        edge_idx_j = edge_j->get_idx();
        AA_ALWAYS_ASSERT(edge_idx_j < m_edge_vec.size());
        AUTO_ASSERT(edge_idx_i != edge_idx_j);
        if((edge_idx_i < edge_idx_j) && (edge_i->intersects(edge_j))){           
            if(NULL != edge_crossover_count_vec){
                ++(edge_crossover_count_vec->at(edge_idx_i));
                ++(edge_crossover_count_vec->at(edge_idx_j));
                }
            ++total_x_count;
            }
        }
    local_edges.clear();
    AUTO_ASSERT(utility_vec.empty());
    }
if(NULL != total_edge_crossover_count){
    *total_edge_crossover_count = total_x_count; }
AA_DECR_CALL_DEPTH();
}

int np01_xy_edge_map::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt=0;
size_t edge_idx, count, loc_grid_node_vec_sz_check;
np01_uint16 w, h, i, j;
size_t gn_vec_idx;
np01_xy_edge *edge;
np01_xy_edge_loc_grid_node *loc_grid_node;
loc_grid_node_vec total_gn_vec_e;
loc_grid_node_vec total_gn_vec_g;
const np01_xy_edge_loc_grid_node *free_gn;

if((NULL != m_owner) && (this != m_owner->get_edge_map())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "this=%x != (owner=%x)->edge_map=%x\n",
        this, m_owner, m_owner->get_edge_map() ); }

for(edge_idx = 0; edge_idx < m_edge_vec.size(); ++edge_idx){
    edge = m_edge_vec.at(edge_idx);
    if(NULL==edge){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_map=%x  edge_idx=%i  NULL edge\n", this, edge_idx ); }
    else{
        err_cnt += edge->verify_data(err_msg, err_msg_capacity, err_msg_pos);
        if(edge->get_idx() != edge_idx){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "(edge=%x)->idx=%i != edge_idx=%i\n",
                edge, edge->get_idx(), edge_idx ); }
        if(edge->get_owner() != this){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge_map (edge=%x:i%i)->owner=%x != this=%x\n",
                edge, edge->get_idx(), edge->get_owner(), this ); }

        /* add every loc grid node to total_gn_vec_e */
        loc_grid_node = edge->get_loc_grid_head_node();
        count = 0;
        while((NULL != loc_grid_node) && (count < m_loc_grid.size())){
            ++count;
            total_gn_vec_e.push_back(loc_grid_node);
            loc_grid_node = loc_grid_node->get_e_next();
            }
        }
    }

w=m_loc_grid_dim.get_w();
h=m_loc_grid_dim.get_h();
err_cnt=m_loc_grid_dim.verify_data(err_msg,err_msg_capacity,err_msg_pos);
loc_grid_node_vec_sz_check = static_cast<size_t>(w) * static_cast<size_t>(h);

if(m_loc_grid.size() != loc_grid_node_vec_sz_check){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "edge_map=%x m_loc_grid.size=%i != check_size=%i\n",
        this ); }

for(gn_vec_idx = 0; gn_vec_idx<m_loc_grid.size(); ++gn_vec_idx){
    i = (h > 0) ? static_cast<np01_uint16>(gn_vec_idx / h) : 0;
    if(i > w){i = (w>0) ? (w-1) : 0; }
    j = (h > 0) ? static_cast<np01_uint16>(gn_vec_idx % h) : 0;

    /* add every loc grid node to total_gn_vec_g */
    loc_grid_node = m_loc_grid.at(gn_vec_idx);
    count = 0;
    while((NULL != loc_grid_node) && (count < m_edge_vec.size())){
        ++count;
        err_cnt += loc_grid_node->verify_data(err_msg, err_msg_capacity,
            err_msg_pos);
        if(this != loc_grid_node->get_edge_map()){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge_map=%x  this != (loc_grid_node=%x:i=%i:j=%i)"
                "->edge_map=%x\n",
                this, loc_grid_node, i, j, loc_grid_node->get_edge_map() ); }
        if(loc_grid_node->get_i() != i){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge_map=%x (loc_grid_node=%x)->i=%i != i=-%i\n",
                this, loc_grid_node, loc_grid_node->get_i(), i ); }
        if(loc_grid_node->get_j() != j){
            ++err_cnt;
            np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
                "edge_map=%x (loc_grid_node=%x)->j=%i != j=-%i\n",
                this, loc_grid_node, loc_grid_node->get_j(), j ); }
        total_gn_vec_g.push_back(loc_grid_node);
        loc_grid_node = loc_grid_node->get_next();
        }
    }

/* compare total_gn_vec_e, total_gn_vec_g */
std::sort(total_gn_vec_e.begin(), total_gn_vec_e.end());
std::sort(total_gn_vec_g.begin(), total_gn_vec_g.end());
if(total_gn_vec_e != total_gn_vec_g){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "edge_map=%x total_gn_vec_e(sz=%i) != total_gn_vec_g(sz=%i)\n",
        this, total_gn_vec_e.size(), total_gn_vec_g.size() ); }
if( std::unique(total_gn_vec_e.begin(), total_gn_vec_e.end()) != 
    total_gn_vec_e.end() ){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "edge_map=%x total_gn_vec_e duplicates found\n",
        this ); }
if( std::unique(total_gn_vec_g.begin(), total_gn_vec_g.end()) != 
    total_gn_vec_g.end() ){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "edge_map=%x total_gn_vec_g duplicates found\n",
        this ); }

count = 0;
free_gn = m_loc_grid_node_free_chain;
while((count < 0xFFFF) && (NULL != free_gn)){
    err_cnt += free_gn->verify_data(err_msg, err_msg_capacity,
        err_msg_pos);
    if(NULL != free_gn->get_owner()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_map=%x  NULL != (free_gn=%x)->owner=%x\n",
            this, free_gn, free_gn->get_owner() ); }
    if(this != free_gn->get_edge_map()){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_map=%x  NULL != (free_gn=%x)->edge_map=%x\n",
            this, free_gn, free_gn->get_edge_map() ); }
    ++count;
    free_gn = free_gn->get_next();
    }

return err_cnt;
}

int np01_xy_edge_map::verify_edges_fully_connected( char *err_msg,
    const size_t err_msg_capacity, size_t *err_msg_pos ) const{
int err_cnt=0;
size_t edge_idx;
for(edge_idx = 0; edge_idx < m_edge_vec.size(); ++edge_idx){
    const np01_xy_edge *edge = m_edge_vec.at(edge_idx);
    if(NULL == edge){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_map=%x  edge_idx=%i  edge=NULL\n", this, edge_idx );
        }
    else if((NULL == edge->get_node_a()) || (NULL == edge->get_node_b())){
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "edge_map=%x  edge_idx=%i  edge=%x  node_a=%x  node_b=%x\n",
            this, edge_idx, edge, edge->get_node_a(), edge->get_node_b() );
        }
    }

return err_cnt;
}


std::ostream& np01_xy_edge_map::ostream_output(std::ostream& os) const{
os << "<xy_edge_map>\n";
os << "<owner>" << std::hex << m_owner << std::dec << "</owner>\n";
m_loc_grid_dim.ostream_output(os);
os << "<loc_grid_node_free_chain>" << std::hex << m_loc_grid_node_free_chain
   << std::dec << "</loc_grid_node_free_chain>\n";
os << "<edge_vec><count>" << m_edge_vec.size() << "</count>\n";
for(edge_vec_citr edge_itr = m_edge_vec.begin(); edge_itr != m_edge_vec.end();
    ++edge_itr){
    const np01_xy_edge *edge = *edge_itr;
    if(NULL == edge){
        os << "</edge=NULL>\n";
        }
    else{
        edge->ostream_output(os);
        }
    }
os << "</edge_vec>\n";
os << "<loc_grid><count>" << m_loc_grid.size() << "</count>\n";
np01_uint16 i = 0, j = 0;
for(loc_grid_node_vec_citr lg_node_itr=m_loc_grid.begin();
    lg_node_itr!=m_loc_grid.end(); ++lg_node_itr){
    const np01_xy_edge_loc_grid_node *lg_node = *lg_node_itr;
    os << "</i=" << i << "></j=" << j << ">";
    if(NULL == lg_node){
        os << "</lg_node=NULL>\n";
        }
    else{
        const size_t lg_count_max = m_edge_vec.size();
        size_t lg_count = 0;
        while((lg_count < lg_count_max) && (NULL != lg_node)){
            ++lg_count;
            lg_node->ostream_output(os);
            lg_node = lg_node->get_next();
            }
        }
    ++j;
    if (j >= m_loc_grid_dim.get_h()){ j = 0; ++i; }
    }
os << "</loc_grid>\n";
os << "</edge_vec>\n";
os << "</xy_edge_map>\n";
return os;
}

std::ostream& np01_xy_edge_map::ostream_output_brief_table(std::ostream& os) const{
os << "<xy_edge_map=" << std::hex << this << std::dec << ">\n";
size_t edge_idx = 0;
for(edge_vec_citr edge_itr = m_edge_vec.begin(); edge_itr != m_edge_vec.end();
    ++edge_itr, ++edge_idx){
    os << "[" << edge_idx << "]";
    const np01_xy_edge *edge = *edge_itr;
    if(NULL == edge){
        os << "NULL\n";
        }
    else{
        const np01_xy_node *node_a = edge->get_node_a();
        const np01_xy_node *node_b = edge->get_node_b();
        const np01_uint32 node_a_idx = (NULL == node_a) ?
            np01_xy_node::get_invalid_idx() : node_a->get_idx();
        const np01_uint32 node_b_idx = (NULL == node_b) ?
            np01_xy_node::get_invalid_idx() : node_b->get_idx();
        os << " " << node_a_idx << "," << node_b_idx << "\n";
        }
    }
os << "</xy_edge_map>\n";
return os;
}


/* take node from free chain or allocate new object */
np01_xy_edge_loc_grid_node *np01_xy_edge_map::alloc_loc_grid_node(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
np01_xy_edge_loc_grid_node *n;
if(NULL == m_loc_grid_node_free_chain){
    n = new np01_xy_edge_loc_grid_node();
    n->set_edge_map(this);
    AUTO_ASSERT(NULL == n->get_owner());
    AUTO_ASSERT(0 == n->get_i());
    AUTO_ASSERT(0 == n->get_j());
    AUTO_ASSERT(NULL == n->get_prev());
    AUTO_ASSERT(NULL == n->get_next());
    AUTO_ASSERT(NULL == n->get_e_prev());
    AUTO_ASSERT(NULL == n->get_e_next());
    }
else{
    n = m_loc_grid_node_free_chain;
    m_loc_grid_node_free_chain = m_loc_grid_node_free_chain->get_next();
    if(NULL != m_loc_grid_node_free_chain){
        m_loc_grid_node_free_chain->set_prev(NULL); }
    n->set_next(NULL);
    AUTO_ASSERT(NULL == n->get_owner());
    AUTO_ASSERT(this == n->get_edge_map());
    AUTO_ASSERT(0 == n->get_i());
    AUTO_ASSERT(0 == n->get_j());
    AUTO_ASSERT(NULL == n->get_prev());
    AUTO_ASSERT(NULL == n->get_e_prev());
    AUTO_ASSERT(NULL == n->get_e_next());
    }
AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
return n;
}

/* add to free chain */
void np01_xy_edge_map::free_loc_grid_node(np01_xy_edge_loc_grid_node *n){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT((NULL == m_loc_grid_node_free_chain) ||
   (m_loc_grid_node_free_chain->get_owner() == NULL) );
AA_ALWAYS_ASSERT(NULL != n);
AUTO_ASSERT(this == n->get_edge_map());
n->set_owner(NULL);
n->set_i(0);
n->set_j(0);
n->set_prev(NULL);
n->set_next(m_loc_grid_node_free_chain);
if(NULL != m_loc_grid_node_free_chain){
    m_loc_grid_node_free_chain->set_prev(n); }
n->set_e_prev(NULL);
n->set_e_next(NULL);
m_loc_grid_node_free_chain = n;
AUTO_ASSERT(0 == n->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}


void np01_xy_edge_map::clear_loc_grid(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
edge_vec_citr edge_itr = m_edge_vec.begin();
for(; edge_itr != m_edge_vec.end(); ++edge_itr){
    np01_xy_edge *edge = *edge_itr;
    remove_edge_from_loc_grid(edge);
    }

if( AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1))
    {
    loc_grid_node_vec_citr lg_itr = m_loc_grid.begin();
    for(; lg_itr != m_loc_grid.end(); ++lg_itr){
        np01_xy_edge_loc_grid_node *lg_node = *lg_itr;
        AA_XDBG_ASSERT(NULL == lg_node,CF01_AA_DEBUG_LEVEL_1 );
        }
    }

m_loc_grid_dim.reset();
m_loc_grid.clear();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_edge_map::free_all_edges(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
edge_vec_citr edge_itr;
np01_xy_edge *edge;
np01_xy_node *na;
np01_xy_node *nb;
for(edge_itr = m_edge_vec.begin(); edge_itr != m_edge_vec.end(); ++edge_itr){
    edge = *edge_itr;
    AA_ALWAYS_ASSERT(NULL != edge);
    remove_edge_from_loc_grid(edge);
    na = edge->get_node_a();
    if(NULL != na) {
        na->disconnect_edge(); }
    nb = edge->get_node_b();
    if(NULL != nb) {
        nb->disconnect_edge(); }
    delete edge;
    }
m_edge_vec.clear();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

void np01_xy_edge_map::free_free_chain(){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
while(NULL != m_loc_grid_node_free_chain){
    np01_xy_edge_loc_grid_node *n = m_loc_grid_node_free_chain->get_next();
    AUTO_ASSERT(NULL == n || n->get_prev()==m_loc_grid_node_free_chain);
    delete m_loc_grid_node_free_chain;
    m_loc_grid_node_free_chain = n;
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}


np01_xy_pair_map::np01_xy_pair_map(): m_node_map_a(new np01_xy_node_map()),
  m_node_map_b(new np01_xy_node_map()), m_edge_map(new np01_xy_edge_map()){
m_node_map_a->set_owner(this);
m_node_map_a->set_xy_group_type(NP01_XY_GROUP_TYPE_A);
m_node_map_b->set_owner(this);
m_node_map_b->set_xy_group_type(NP01_XY_GROUP_TYPE_B);
m_edge_map->set_owner(this);
}

np01_xy_pair_map::~np01_xy_pair_map(){
AA_INCR_CALL_DEPTH();
/* disconnect and delete all edges */
/* delete all nodes */
delete m_edge_map;
delete m_node_map_b;
delete m_node_map_a;
AA_DECR_CALL_DEPTH();
}

/* initialize node maps A & B, including locator grir.  Initialize edge
locator grid */
void np01_xy_pair_map::init_ab(
    const np01_xy_pair_map_init_ab_params *init_params){
AA_INCR_CALL_DEPTH();
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(NULL != m_node_map_a);
AA_ALWAYS_ASSERT(NULL != m_node_map_b);
AA_ALWAYS_ASSERT(NULL != m_edge_map);
AA_ALWAYS_ASSERT(NULL != init_params);
AUTO_ASSERT(init_params->loc_grid_density > 0.0);
AUTO_ASSERT(init_params->max_loc_grid_sq_count > 0);
AUTO_ASSERT(m_node_map_a->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_A == m_node_map_a->get_xy_group_type());
AUTO_ASSERT(m_node_map_b->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_B == m_node_map_b->get_xy_group_type());
AUTO_ASSERT(m_edge_map->get_owner() == this);

const size_t point_count_a = (init_params->xy_vec_a).size();
const size_t point_count_b = (init_params->xy_vec_b).size();
const size_t point_count_ab = point_count_a + point_count_b;
if((point_count_ab > 0) && (init_params->loc_grid_density > 0.0)){
    /* bounding box of all points */
    np01_float64 x,y;
    np01_float64 x_min=std::numeric_limits<np01_float64>::max(); 
    np01_float64 y_min=std::numeric_limits<np01_float64>::max(); 
    np01_float64 x_max=-std::numeric_limits<np01_float64>::max(); 
    np01_float64 y_max=-std::numeric_limits<np01_float64>::max(); 
    np01_xy_vec_citr xy_itr;
    np01_xy_vec_citr xy_end_itr = (init_params->xy_vec_a).end();
    for(xy_itr = (init_params->xy_vec_a).begin(); xy_itr != xy_end_itr; ++xy_itr){
        x=xy_itr->first;
        y=xy_itr->second;
        if(x < x_min){x_min=x;}
        else if(x > x_max){x_max=x;}
        if(y < y_min){y_min=y;}
        else if(y > y_max){y_max=y;}
        }
    xy_end_itr = (init_params->xy_vec_b).end();
    for(xy_itr = (init_params->xy_vec_b).begin(); xy_itr != xy_end_itr; ++xy_itr){
        x=xy_itr->first;
        y=xy_itr->second;
        if(x < x_min){x_min=x;}
        else if(x > x_max){x_max=x;}
        if(y < y_min){y_min=y;}
        else if(y > y_max){y_max=y;}
        }

    /* initialize locator grids */
    np01_xy_loc_grid_dim_init_params loc_grid_dim_init_params = {};
    loc_grid_dim_init_params.point_count = point_count_ab;            
    loc_grid_dim_init_params.bb_min_xy.first = x_min; 
    loc_grid_dim_init_params.bb_min_xy.second = y_min; 
    loc_grid_dim_init_params.bb_max_xy.first = x_max; 
    loc_grid_dim_init_params.bb_max_xy.second = y_max; 
    loc_grid_dim_init_params.loc_grid_density = init_params->loc_grid_density;
    loc_grid_dim_init_params.max_loc_grid_sq_count =
        init_params->max_loc_grid_sq_count;  
    np01_xy_loc_grid_dim loc_grid_dim;
    loc_grid_dim.init(&loc_grid_dim_init_params);
    m_node_map_a->init_loc_grid(loc_grid_dim);
    m_node_map_b->init_loc_grid(loc_grid_dim);
    m_edge_map->init_loc_grid(loc_grid_dim);

    /* reserve size for nodes, edges */
    m_node_map_a->reserve(static_cast<np01_uint32>(point_count_a));
    m_node_map_b->reserve(static_cast<np01_uint32>(point_count_b));
    m_edge_map->reserve(static_cast<np01_uint32>(
        (point_count_a < point_count_b) ? point_count_a : point_count_b));

    /* create nodes */
    xy_end_itr = (init_params->xy_vec_a).end();
    for(xy_itr = (init_params->xy_vec_a).begin(); xy_itr != xy_end_itr; ++xy_itr){
        const np01_xy& xy_a = *xy_itr;
        m_node_map_a->create_node(xy_a);
        }
    xy_end_itr = (init_params->xy_vec_b).end();
    for(xy_itr = (init_params->xy_vec_b).begin(); xy_itr != xy_end_itr; ++xy_itr){
        const np01_xy& xy_b = *xy_itr;
        m_node_map_b->create_node(xy_b);
        }
    }
AA_XDBG_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_DECR_CALL_DEPTH();
}

/* for each point in smaller set, find nearset point in larger set. create edges. */
void np01_xy_pair_map::init_edges_ctr_out(const np01_float64 dsq_noise_ampl = 0.0){
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_ALWAYS_ASSERT(NULL != m_node_map_a);
AA_ALWAYS_ASSERT(NULL != m_node_map_b);
AA_ALWAYS_ASSERT(NULL != m_edge_map);
AUTO_ASSERT(m_node_map_a->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_A == m_node_map_a->get_xy_group_type());
AUTO_ASSERT(m_node_map_b->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_B == m_node_map_b->get_xy_group_type());
AUTO_ASSERT(m_edge_map->get_owner() == this);
AUTO_ASSERT(m_edge_map->get_edge_count() == 0);
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));

np01_float64 dx,dy,dsq;
np01_uint32 rand_uint32 = 4294967291ul; /* http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php */
int rand_int = 0;

/* distinguish small and large point sets */
np01_xy_node_map *node_map_sm, *node_map_lg;
np01_uint32 small_node_count;
np01_uint32 i;
np01_xy_node *node_sm=NULL;
if( m_node_map_a->get_node_count() <= m_node_map_b->get_node_count() ){
    node_map_sm = m_node_map_a;
    node_map_lg = m_node_map_b;
    small_node_count = m_node_map_a->get_node_count();
    }
else{
    node_map_sm = m_node_map_b;
    node_map_lg = m_node_map_a;
    small_node_count = m_node_map_b->get_node_count();
    }
AUTO_ASSERT(node_map_sm->get_node_count() <= node_map_lg->get_node_count());

/* reserve memory */
m_edge_map->reserve(small_node_count);

/* find centroid of small point set*/
np01_xy xy_sum_sm(0.0, 0.0);
for(i = 0; i < small_node_count; ++i){
    node_sm = node_map_sm->get_node_by_idx(i);
    xy_sum_sm.first += node_sm->get_x();
    xy_sum_sm.second += node_sm->get_y();
    }
np01_xy xy_ctr_sm(0.0, 0.0);
if(small_node_count > 0){
    xy_ctr_sm.first = xy_sum_sm.first/
        static_cast<np01_float64>(small_node_count);
    xy_ctr_sm.second = xy_sum_sm.second/
        static_cast<np01_float64>(small_node_count);
    }

/* sort points from centroid outward */
typedef std::pair<np01_float64, np01_uint32> flt_idx_pair;
typedef std::vector<flt_idx_pair> flt_idx_pair_vec; 
typedef flt_idx_pair_vec::const_iterator flt_idx_pair_vec_citr;
typedef flt_idx_pair_vec::iterator flt_idx_pair_vec_itr;

flt_idx_pair_vec dsq_idx_vec_sm;
dsq_idx_vec_sm.reserve(small_node_count);
for(i = 0; i < small_node_count; ++i){
    node_sm = node_map_sm->get_node_by_idx(i);
    dx = node_sm->get_x() - xy_ctr_sm.first;
    dy = node_sm->get_y() - xy_ctr_sm.second;
    dsq = (dx*dx) + (dy*dy);
    if(dsq_noise_ampl > 0.0){
        /* advance random number, pg 359, The Standard C Library (c)1991,
        P. J. Plauger, ISBN 0-13-131509-9 */
        rand_uint32 = (rand_uint32 * 1103515245) + 12345;
        rand_int = static_cast<int>(
           (static_cast<unsigned int>(rand_uint32 >> 16)) & RAND_MAX);
        dsq += (dsq_noise_ampl * static_cast<np01_float64>(rand_int))/
            static_cast<np01_float64>(RAND_MAX);
        }
    dsq_idx_vec_sm.push_back(flt_idx_pair(dsq,i));
    }
std::sort(dsq_idx_vec_sm.begin(), dsq_idx_vec_sm.end());

/* for each point in small set, from centroid outward, find nearest unused
point in large set, join with edge. */
np01_xy_node_map::node_vec nodes_lg;
const np01_xy_loc_grid_dim& lg_dim = node_map_lg->get_loc_grid_dim();
AUTO_ASSERT(lg_dim.get_w() >= 1);
AUTO_ASSERT(lg_dim.get_h() >= 1);
const np01_float64& sq_sz = lg_dim.get_sq_size();
const np01_float64 min_search_r = sq_sz / 2.0;
const np01_float64 max_search_r = sq_sz * (lg_dim.get_w() + lg_dim.get_h()) / 2.0;
flt_idx_pair_vec_citr dsq_idx_itr = dsq_idx_vec_sm.begin();
for(; dsq_idx_itr != dsq_idx_vec_sm.end(); ++dsq_idx_itr){
    node_sm = node_map_sm->get_node_by_idx(dsq_idx_itr->second);
    AA_ALWAYS_ASSERT(NULL != node_sm);
    np01_xy_node *near_node_lg = NULL;
    const np01_xy xy_node_sm(node_sm->get_x(), node_sm->get_y());
    np01_float64 search_r = min_search_r;
    enum search_mode_type {srch_md_circle, srch_md_all, srch_md_done};
    search_mode_type search_mode = srch_md_circle;
    while(srch_md_done != search_mode){
        AUTO_ASSERT(nodes_lg.empty());
        AUTO_ASSERT(NULL == near_node_lg);
        /* search in circle in large set, or get all the large set  */
        switch(search_mode){
            case srch_md_circle:
                node_map_lg->get_nodes_in_circle(xy_node_sm,search_r,
                    &nodes_lg);
                break;
            case srch_md_all:
                node_map_lg->get_all_nodes(&nodes_lg);
                break;
            case srch_md_done:
            default:
                AA_ALWAYS_ASSERT(false);
                break;
            }

        /* find closest node that is not already connected to an edge */
        np01_float64 near_dsq = std::numeric_limits<np01_float64>::max();
        np01_xy_node_map::node_vec_citr node_lg_itr = nodes_lg.begin();
        for(; node_lg_itr != nodes_lg.end(); ++node_lg_itr){
            np01_xy_node *node_lg = *node_lg_itr;
            AA_ALWAYS_ASSERT(NULL != node_lg);
            if(node_lg->get_edge() == NULL){
                dx = node_lg->get_x() - node_sm->get_x();
                dy = node_lg->get_y() - node_sm->get_y();
                dsq = (dx*dx) + (dy*dy);
                if((NULL == near_node_lg) || (dsq < near_dsq) ||
                  ((dsq==near_dsq)&&(node_lg->get_idx()<near_node_lg->get_idx()))){
                    near_node_lg = node_lg;
                    near_dsq = dsq;
                    }
                }
            }
        nodes_lg.clear();

        if(NULL == near_node_lg){
            switch(search_mode){
                case srch_md_circle:
                    search_r *= 2.0; /* double search radius */
                    if(search_r > max_search_r){
                        search_mode = srch_md_all; }
                    break;
                case srch_md_all:
                    search_mode=srch_md_done;
                    break;
                case srch_md_done:
                default:
                    AA_ALWAYS_ASSERT(false);
                    search_mode=srch_md_done;
                    break;
                }
            }
        else{
            search_mode=srch_md_done; }

        /* join nodes with edge */
        if(NULL != near_node_lg){
            np01_xy_node *node_a;
            np01_xy_node *node_b;
            if(node_sm->get_xy_group_type() == NP01_XY_GROUP_TYPE_A){
                node_a = node_sm;
                node_b = near_node_lg;
                }
            else{
                node_a = near_node_lg;
                node_b = node_sm;
                }
            np01_xy_edge *edge = m_edge_map->create_edge(node_a, node_b);
            AA_ALWAYS_ASSERT(NULL != edge);
            AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),
                AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
            }
        }
    }

AUTO_ASSERT(m_edge_map->get_edge_count() == small_node_count);
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()));
AA_DECR_CALL_DEPTH();
}

/* create edges. */
void np01_xy_pair_map::init_edges_far_near(){
AA_INCR_CALL_DEPTH();
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));
AA_ALWAYS_ASSERT(NULL != m_node_map_a);
AA_ALWAYS_ASSERT(NULL != m_node_map_b);
AA_ALWAYS_ASSERT(NULL != m_edge_map);
AUTO_ASSERT(m_node_map_a->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_A == m_node_map_a->get_xy_group_type());
AUTO_ASSERT(m_node_map_b->get_owner() == this);
AUTO_ASSERT(NP01_XY_GROUP_TYPE_B == m_node_map_b->get_xy_group_type());
AUTO_ASSERT(m_edge_map->get_owner() == this);
AUTO_ASSERT(m_edge_map->get_edge_count() == 0);
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR() ));

np01_float64 x,y,dx,dy,dsq,d;

/* distinguish small and large point sets */
np01_xy_node_map *node_map_sm, *node_map_lg;
np01_uint32 small_node_count;
np01_uint32 large_node_count;
np01_uint32 i,j;
np01_xy_node *node_sm, *node_lg;
np01_xy_node_map::node_vec utility_vec;
if( m_node_map_a->get_node_count() <= m_node_map_b->get_node_count() ){
    node_map_sm = m_node_map_a;
    node_map_lg = m_node_map_b;
    small_node_count = m_node_map_a->get_node_count();
    large_node_count = m_node_map_b->get_node_count();
    }
else{
    node_map_sm = m_node_map_b;
    node_map_lg = m_node_map_a;
    small_node_count = m_node_map_b->get_node_count();
    large_node_count = m_node_map_a->get_node_count();
    }
utility_vec.reserve(large_node_count);
AUTO_ASSERT(node_map_sm->get_node_count() <= node_map_lg->get_node_count());

/* reserve memory */
m_edge_map->reserve(small_node_count);

/* sort large points by distance to nearest small point.
find centroid of large point set */
typedef std::pair<np01_float64, np01_uint32> flt_idx_pair;
typedef std::vector<flt_idx_pair> flt_idx_pair_vec; 
typedef flt_idx_pair_vec::const_iterator flt_idx_pair_vec_citr;
typedef flt_idx_pair_vec::iterator flt_idx_pair_vec_itr;
np01_float64 near_d;
flt_idx_pair_vec near_d_idx_vec_lg;
near_d_idx_vec_lg.reserve(large_node_count);
np01_xy xy_sum_lg(0.0, 0.0);
for(i = 0; i < large_node_count; ++i){
    node_lg = node_map_lg->get_node_by_idx(i);
    near_d = 0.0;
    x = node_lg->get_x();
    y = node_lg->get_y();
    node_map_sm->get_near_node(np01_xy(x,y), &near_d, &utility_vec);
    near_d_idx_vec_lg.push_back(flt_idx_pair(near_d,i));
    xy_sum_lg.first += (node_lg->get_x());
    xy_sum_lg.second += (node_lg->get_y());
    }
std::sort(near_d_idx_vec_lg.begin(), near_d_idx_vec_lg.end());
np01_xy xy_ctr_lg(0.0,0.0); /* centroid */
if(large_node_count > 0){
    xy_ctr_lg.first =
        xy_sum_lg.first/static_cast<np01_float64>(large_node_count);
    xy_ctr_lg.second =
        xy_sum_lg.second/static_cast<np01_float64>(large_node_count);
    }

/* Using only the nearest <small_node_count> nodes of large node set,
initialize: usable large node table, large point-to-nearest small point
distance table, large point-to-centroid distance vector */
np01_float64_vec d_ctr_vec_lg; /* distances to centroid */
np01_float64 d_ctr_lg_sum = 0.0; /* sum of distances to centroid */
np01_float64_vec d_near_vec_lg; /* distances to nearest small node */
np01_float64 d_near_lg_sum = 0.0; /* sum of distances to nearest small node */
np01_xy_node_map::node_vec usable_node_lg_lookup_table;
d_ctr_vec_lg.resize(large_node_count, 0.0);
d_near_vec_lg.resize(large_node_count, 0.0);
usable_node_lg_lookup_table.reserve(small_node_count);
for(j = 0; j < small_node_count; ++j){
    const flt_idx_pair& di = near_d_idx_vec_lg.at(j);
    i = di.second;
    d_near_vec_lg.at(i) = di.first;
    d_near_lg_sum += di.first;
    node_lg = node_map_lg->get_node_by_idx(i);
    usable_node_lg_lookup_table.push_back(node_lg);
    x = node_lg->get_x();
    y = node_lg->get_y();
    dx = x - xy_ctr_lg.first;
    dy = y - xy_ctr_lg.second;
    dsq = (dx*dx) + (dy*dy);
    d = sqrt(dsq);
    d_ctr_vec_lg.at(i) = d;
    d_ctr_lg_sum += d;
    }
std::sort(usable_node_lg_lookup_table.begin(),
    usable_node_lg_lookup_table.end());
const np01_float64 d_ctr_lg_ave = (small_node_count > 0) ?
    d_ctr_lg_sum/static_cast<np01_float64>(small_node_count) : 0.0;
const np01_float64 d_near_lg_ave = (small_node_count > 0) ?
    d_near_lg_sum/static_cast<np01_float64>(small_node_count) : 0.0;

/* check that usable_node_lg_lookup_table contains unique, sorted values */
AUTO_ASSERT(usable_node_lg_lookup_table.size() == small_node_count);
if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1)){
    for(size_t ii = 1; ii < usable_node_lg_lookup_table.size(); ++ii){
        const np01_xy_node *p = usable_node_lg_lookup_table.at(ii-1);
        const np01_xy_node *n = usable_node_lg_lookup_table.at(ii);
        AA_XDBG_ASSERT(NULL != p, CF01_AA_DEBUG_LEVEL_1);
        AA_XDBG_ASSERT(NULL != n, CF01_AA_DEBUG_LEVEL_1);
        AA_XDBG_ASSERT(p->get_owner() == node_map_lg, CF01_AA_DEBUG_LEVEL_1);
        AA_XDBG_ASSERT(n->get_owner() == node_map_lg, CF01_AA_DEBUG_LEVEL_1);
        AA_XDBG_ASSERT(p < n, CF01_AA_DEBUG_LEVEL_1);
        }
    }

/* find centroid of small point set*/
np01_xy xy_sum_sm(0.0, 0.0);
for(i = 0; i < small_node_count; ++i){
    node_sm = node_map_sm->get_node_by_idx(i);
    xy_sum_sm.first += node_sm->get_x();
    xy_sum_sm.second += node_sm->get_y();
    }
np01_xy xy_ctr_sm(0.0, 0.0);
if(small_node_count > 0){
    xy_ctr_sm.first = xy_sum_sm.first/
        static_cast<np01_float64>(small_node_count);
    xy_ctr_sm.second = xy_sum_sm.second/
        static_cast<np01_float64>(small_node_count);
    }

/* initialize: small point-to-nearest large point distance table,
small point-to-centroid distance vector */
np01_float64_vec d_ctr_vec_sm; /* distances to centroid */
np01_float64 d_ctr_sm_sum = 0.0; /* sum of distances to centroid */
np01_float64_vec d_near_vec_sm; /* distances to nearest large node */
np01_float64 d_near_sm_sum = 0.0; /* sum of distances to nearest large node */
d_ctr_vec_sm.resize(small_node_count, 0.0);
d_near_vec_sm.resize(small_node_count, 0.0);
for(i = 0; i < small_node_count; ++i){
    node_sm = node_map_sm->get_node_by_idx(i);
    near_d = 0.0;
    x = node_sm->get_x();
    y = node_sm->get_y();
    node_map_lg->get_near_node(np01_xy(x,y), &near_d, &utility_vec);
    d_near_vec_sm.at(i) = near_d;
    d_near_sm_sum += near_d;
    dx = x - xy_ctr_sm.first;
    dy = y - xy_ctr_sm.second;
    dsq = (dx*dx) + (dy*dy);
    d = sqrt(dsq);
    d_ctr_vec_sm.at(i) = d;
    d_ctr_sm_sum += d;
    }
const np01_float64 d_ctr_sm_ave = (small_node_count > 0) ?
    d_ctr_sm_sum/static_cast<np01_float64>(small_node_count) : 0.0;
const np01_float64 d_near_sm_ave = (small_node_count > 0) ?
    d_near_sm_sum/static_cast<np01_float64>(small_node_count) : 0.0;

/* put both large and small nodes into a single sorted vector */
typedef std::pair<np01_float64,np01_xy_group_type> fv_gt_pair;
typedef std::pair<np01_uint32,np01_xy_node *> idx_node_pair;
typedef std::pair<fv_gt_pair, idx_node_pair> fv_gt_idx_node;
typedef std::vector<fv_gt_idx_node> fv_gt_idx_node_vec;
np01_float64 far_val;
fv_gt_idx_node_vec near_far_vec;
near_far_vec.reserve(2 * small_node_count);

/* put large nodes into near-far vector */
const np01_float64 d_ctr_lg_ave_pos=(d_ctr_lg_ave >0.0) ? d_ctr_lg_ave : 1.0;
const np01_float64 d_near_lg_ave_pos=(d_near_lg_ave>0.0) ? d_near_lg_ave : 1.0;
for(j = 0; j < small_node_count; ++j){
    const flt_idx_pair& di = near_d_idx_vec_lg.at(j);
    i = di.second;
    node_lg = node_map_lg->get_node_by_idx(i);
    d = d_ctr_vec_lg.at(i);
    near_d = d_near_vec_lg.at(i);
    far_val = (d/d_ctr_lg_ave_pos) + (near_d/d_near_lg_ave_pos);
    near_far_vec.push_back(fv_gt_idx_node(
        fv_gt_pair(far_val, node_lg->get_xy_group_type()),
        idx_node_pair(i, node_lg)));
    }

/* put small nodes into near-far vector */
const np01_float64 d_ctr_sm_ave_pos=(d_ctr_sm_ave >0.0) ? d_ctr_sm_ave : 1.0;
const np01_float64 d_near_sm_ave_pos=(d_near_lg_ave>0.0) ? d_near_sm_ave : 1.0;
for(i = 0; i < small_node_count; ++i){
    node_sm = node_map_sm->get_node_by_idx(i);
    d = d_ctr_vec_sm.at(i);
    near_d = d_near_vec_sm.at(i);
    far_val = (d/d_ctr_sm_ave_pos) + (near_d/d_near_sm_ave_pos);
    near_far_vec.push_back(fv_gt_idx_node(
        fv_gt_pair(far_val, node_sm->get_xy_group_type()),
        idx_node_pair(i, node_sm)));
    }

/* sort near-far vector */
AUTO_ASSERT(near_far_vec.size() == (2*small_node_count));
std::sort(near_far_vec.begin(), near_far_vec.end());

/* check that near_far_vec contains unique nodes */
AUTO_ASSERT(near_far_vec.size() == (2*small_node_count));
if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1)){
    std::set<const np01_xy_node*> node_set;
    size_t sm_count = 0;
    size_t lg_count = 0;
    for(size_t ii = 0; ii < near_far_vec.size(); ++ii){
        const fv_gt_idx_node& k = near_far_vec.at(ii);
        const np01_xy_node *nn = k.second.second;
        AA_ALWAYS_ASSERT(NULL != nn);
        if(nn->get_owner() == node_map_sm){ ++sm_count; }
        else if(nn->get_owner() == node_map_lg){ ++lg_count; } 
        else{AA_ALWAYS_ASSERT(false); }
        AA_XDBG_ASSERT(nn->get_edge() == NULL, CF01_AA_DEBUG_LEVEL_1);
        std::pair<std::set<const np01_xy_node*>::iterator,bool>
            insert_result = node_set.insert(nn);
        AA_XDBG_ASSERT(insert_result.second, CF01_AA_DEBUG_LEVEL_1);
        }
    AA_XDBG_ASSERT(sm_count == small_node_count, CF01_AA_DEBUG_LEVEL_1);
    AA_XDBG_ASSERT(lg_count == small_node_count, CF01_AA_DEBUG_LEVEL_1);
    AA_XDBG_ASSERT(node_set.size() == (2*small_node_count),
        CF01_AA_DEBUG_LEVEL_1);
    }

/* assign points far to near */
j = static_cast<np01_uint32>(near_far_vec.size());
while(j > 0){
    --j;
    const fv_gt_idx_node& fgin = near_far_vec.at(j);
    np01_xy_node *node = fgin.second.second;
    AA_ALWAYS_ASSERT(NULL != node);
    AUTO_ASSERT(node->get_idx() == fgin.second.first);
    AA_ALWAYS_ASSERT(node->get_xy_group_type() == fgin.first.second);
    if(node->get_edge() == NULL){
        const np01_xy xy(node->get_x(), node->get_y());
        if(node->get_owner() == node_map_sm){
            node_sm = node;
            node_lg = node_map_lg->get_near_unconnected_node(xy,
                &usable_node_lg_lookup_table, &utility_vec);
            AA_ALWAYS_ASSERT(NULL != node_lg);
            AA_ALWAYS_ASSERT(node_lg->get_owner() == node_map_lg);
            AA_ALWAYS_ASSERT(node_lg->get_edge() == NULL);
            }
        else{
            AA_ALWAYS_ASSERT(node->get_owner() == node_map_lg);
            node_lg = node;
            node_sm = node_map_sm->get_near_unconnected_node(xy, NULL,
                &utility_vec);
            AA_ALWAYS_ASSERT(NULL != node_sm);
            AA_ALWAYS_ASSERT(node_sm->get_owner() == node_map_sm);
            AA_ALWAYS_ASSERT(node_sm->get_edge() == NULL);
            }
        np01_xy_edge *edge = NULL;
        switch(node_sm->get_xy_group_type()){
            case NP01_XY_GROUP_TYPE_A:
                edge = m_edge_map->create_edge(node_sm, node_lg);
                break;
            case NP01_XY_GROUP_TYPE_B:
                edge = m_edge_map->create_edge(node_lg, node_sm);
                break;
            case NP01_XY_GROUP_TYPE_UNKNOWN:
            default:
                edge = NULL;
                AA_ALWAYS_ASSERT(false);
                break;
            }
        AA_ALWAYS_ASSERT(NULL != edge);
        AUTO_ASSERT(0 == edge->verify_data(AA_ERR_BUF(),
            AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
        }
    }

/* check that all nodes in near_far_vec were assigned */
AUTO_ASSERT(near_far_vec.size() == (2*small_node_count));
if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1)){
    for(size_t ii = 0; ii < near_far_vec.size(); ++ii){
        const fv_gt_idx_node& k = near_far_vec.at(ii);
        const np01_xy_node *nn = k.second.second;
        AA_ALWAYS_ASSERT(NULL != nn);
        AA_XDBG_ASSERT(nn->get_edge() != NULL, CF01_AA_DEBUG_LEVEL_1);
        }
    }

AUTO_ASSERT(m_edge_map->get_edge_count() == small_node_count);
AUTO_ASSERT(0 == verify_data(AA_ERR_BUF(), AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()));
AA_DECR_CALL_DEPTH();
}

int np01_xy_pair_map::verify_data( char *err_msg, const size_t err_msg_capacity,
    size_t *err_msg_pos ) const {
int err_cnt = 0;
if( NULL != m_node_map_a ) {
    err_cnt += m_node_map_a->verify_data(err_msg,err_msg_capacity,
        err_msg_pos);
    if( m_node_map_a->get_owner() != this ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_node_map_a->get_owner() != this\n"); }
    if( NP01_XY_GROUP_TYPE_A != m_node_map_a->get_xy_group_type() ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_node_map_a  group type != A\n"); }
    }

if( NULL != m_node_map_b ) {
    err_cnt += m_node_map_b->verify_data(err_msg,err_msg_capacity,
        err_msg_pos);
    if( m_node_map_b->get_owner() != this ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_node_map_b->get_owner() != this\n"); }
    if( NP01_XY_GROUP_TYPE_B != m_node_map_b->get_xy_group_type() ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_node_map_b  group type != B\n"); }
    }

if( NULL != m_edge_map ) {
    err_cnt += m_edge_map->verify_data(err_msg,err_msg_capacity,
        err_msg_pos);
    if( m_edge_map->get_owner() != this ) {
        ++err_cnt;
        np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
            "m_edge_map->get_owner() != this\n"); }
    }

if( (NULL != m_node_map_a) && (NULL != m_node_map_b) &&
    (m_node_map_a == m_node_map_b)){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "m_node_map_a == m_node_map_b\n"); 
    }
if( (NULL != m_node_map_a) && (NULL != m_node_map_b) &&
    (m_node_map_a->get_loc_grid_dim() != m_node_map_b->get_loc_grid_dim())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "m_node_map_a,m_node_map_b different locator grid dimensions\n");
    }
if( (NULL != m_node_map_a) && (NULL != m_edge_map) &&
    (m_node_map_a->get_loc_grid_dim() != m_edge_map->get_loc_grid_dim())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "m_node_map_a,m_edge_map different locator grid dimensions\n");
    }
if( (NULL != m_node_map_b) && (NULL != m_edge_map) &&
    (m_node_map_b->get_loc_grid_dim() != m_edge_map->get_loc_grid_dim())){
    ++err_cnt;
    np01_snprintf(err_msg, err_msg_capacity, err_msg_pos,
        "m_node_map_b,m_edge_map different locator grid dimensions\n");
    }
return err_cnt; 
}

std::ostream& np01_xy_pair_map::ostream_output(std::ostream& os) const{
os << "<xy_pair_map>\n";
if(NULL == m_node_map_a){
    os << "<node_map_a=NULL>\n";}
else{
    os << "<node_map_a>\n";
    m_node_map_a->ostream_output(os);
    os << "</node_map_a>\n";
    }
if(NULL == m_node_map_b){
    os << "<node_map_b=NULL>\n";}
else{
    os << "<node_map_b>\n";
    m_node_map_b->ostream_output(os);
    os << "</node_map_b>\n";
    }
if(NULL == m_node_map_a){
    os << "<edge_map=NULL>\n";}
else{
    m_edge_map->ostream_output(os);
    }
os << "</xy_pair_map>\n";

return os;
}


void np01_xy_pair_map::write_bmp_file(const char *file_name) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != m_node_map_a);
AA_ALWAYS_ASSERT(NULL != m_node_map_b);
AA_ALWAYS_ASSERT(NULL != m_edge_map);

np01_uint32 i,j;
np01_uint32 ii,jj,ii2,jj2;
np01_float64 x,xa,xb,y,ya,yb;
const np01_xy_node *node_a, *node_b;
const np01_xy_edge *edge;
np01_float64 x_shp_min = std::numeric_limits<np01_float64>::max();
np01_float64 x_shp_max = -std::numeric_limits<np01_float64>::max();
np01_float64 y_shp_min = std::numeric_limits<np01_float64>::max();
np01_float64 y_shp_max = -std::numeric_limits<np01_float64>::max();

for(i = 0; i < 3; ++i){
    np01_xy_loc_grid_dim loc_grid_dim;
    switch(i){
        case 0: loc_grid_dim = m_node_map_a->get_loc_grid_dim(); break;
        case 1: loc_grid_dim = m_node_map_b->get_loc_grid_dim(); break;
        case 2: loc_grid_dim = m_edge_map->get_loc_grid_dim(); break;
        default: break;
        }
    if((loc_grid_dim.get_w() > 0) && (loc_grid_dim.get_h() > 0)){
        x=loc_grid_dim.get_x_min() + (loc_grid_dim.get_sq_size() *
            static_cast<np01_float64>(loc_grid_dim.get_w()));
        y=loc_grid_dim.get_y_min() + (loc_grid_dim.get_sq_size() *
            static_cast<np01_float64>(loc_grid_dim.get_h()));
        if(loc_grid_dim.get_x_min()<x_shp_min){x_shp_min=loc_grid_dim.get_x_min();}
        if(loc_grid_dim.get_y_min()<y_shp_min){y_shp_min=loc_grid_dim.get_y_min();}
        if(x>x_shp_max){x_shp_max=x;}
        if(y>y_shp_max){y_shp_max=y;}
        }
    }

const np01_uint32 node_count_a = m_node_map_a->get_node_count();
const np01_uint32 node_count_b = m_node_map_b->get_node_count();
const np01_uint32 edge_count = m_edge_map->get_edge_count();
const np01_uint32 shape_count = node_count_a + node_count_b + edge_count;

for(i=0; i < node_count_a; ++i){
    node_a = m_node_map_a->get_node_by_idx(i);
    AA_ALWAYS_ASSERT(NULL != node_a);
    xa = node_a->get_x();
    ya = node_a->get_y();
    if(xa<x_shp_min){x_shp_min=xa;}
    if(xa>x_shp_max){x_shp_max=xa;}
    if(ya<y_shp_min){y_shp_min=ya;}
    if(ya>y_shp_max){y_shp_max=ya;}
    }

for(i=0; i < node_count_b; ++i){
    node_b = m_node_map_b->get_node_by_idx(i);
    AA_ALWAYS_ASSERT(NULL != node_b);
    xb = node_b->get_x();
    yb = node_b->get_y();
    if(xb<x_shp_min){x_shp_min=xb;}
    if(xb>x_shp_max){x_shp_max=xb;}
    if(yb<y_shp_min){y_shp_min=yb;}
    if(yb>y_shp_max){y_shp_max=yb;}
    }

np01_float64 dx_shp = (x_shp_max > x_shp_min) ? x_shp_max - x_shp_min : 0.0;
np01_float64 dy_shp = (y_shp_max > y_shp_min) ? y_shp_max - y_shp_min : 0.0;
np01_float64 area_shp = dx_shp * dy_shp;
np01_float64 scale_factor = 1.0;
if((area_shp > 0.0) && (shape_count > 0.0)){
    np01_float64 unit_area = area_shp/static_cast<np01_float64>(shape_count);
    np01_float64 unit_len = sqrt(unit_area);
    np01_float64 width_units = dx_shp/unit_len;
    np01_float64 height_units = dy_shp/unit_len;
    static const np01_float64 target_pixels_per_unit_len = 16.0;
    np01_float64 target_width_pixels=width_units*target_pixels_per_unit_len;
    np01_float64 target_height_pixels=height_units*target_pixels_per_unit_len;
    static const np01_float64 max_bmp_pixels_shp_w_h = 4000.0;
    if(target_width_pixels > max_bmp_pixels_shp_w_h){
        target_width_pixels = max_bmp_pixels_shp_w_h; }
    if(target_height_pixels > max_bmp_pixels_shp_w_h){
        target_height_pixels = max_bmp_pixels_shp_w_h; }
    np01_float64 scale_factor_w = target_width_pixels/dx_shp; 
    np01_float64 scale_factor_h = target_height_pixels/dy_shp;
    scale_factor = (scale_factor_w < scale_factor_h) ?
        scale_factor_h : scale_factor_h;
    }
AA_ALWAYS_ASSERT(scale_factor > 0.0);


static const int32_t pixel_border_width = 8;
np01_bmp_file_init_params bmp_init_params;
bmp_init_params.width_px = (2*pixel_border_width) + lrint(scale_factor*dx_shp);
bmp_init_params.height_px = (2*pixel_border_width) + lrint(scale_factor*dy_shp);
np01_bmp_file bmp_file(bmp_init_params);

const np01_xy_loc_grid_dim& loc_grid_dim_e = m_edge_map->get_loc_grid_dim();
for(i = 0; i <= loc_grid_dim_e.get_w(); ++i){
    x=loc_grid_dim_e.get_x_min() + 
        (static_cast<np01_float64>(i) * loc_grid_dim_e.get_sq_size());
    y=loc_grid_dim_e.get_y_min() + 
        (static_cast<np01_float64>(loc_grid_dim_e.get_h()) *
            loc_grid_dim_e.get_sq_size());
    ii = pixel_border_width +
        static_cast<np01_uint32>(floor((x-x_shp_min)*scale_factor));
    jj = pixel_border_width + static_cast<np01_uint32>(
        floor((loc_grid_dim_e.get_y_min()-y_shp_min)*scale_factor));
    jj2 = pixel_border_width +
        static_cast<np01_uint32>(floor((y-y_shp_min)*scale_factor));
    bmp_file.draw_line(ii,jj,ii,jj2,np01::np01_bmp_color(96,96,96));
    }
for(j = 0; j <= loc_grid_dim_e.get_h(); ++j){
    y=loc_grid_dim_e.get_y_min() + 
        (static_cast<np01_float64>(j) * loc_grid_dim_e.get_sq_size());
    x=loc_grid_dim_e.get_x_min() + 
        (static_cast<np01_float64>(loc_grid_dim_e.get_w()) *
            loc_grid_dim_e.get_sq_size());
    jj = pixel_border_width +
        static_cast<np01_uint32>(floor((y-y_shp_min)*scale_factor));
    ii = pixel_border_width + static_cast<np01_uint32>(
        floor((loc_grid_dim_e.get_x_min()-x_shp_min)*scale_factor));
    ii2 = pixel_border_width +
        static_cast<np01_uint32>(floor((x-x_shp_min)*scale_factor));
    bmp_file.draw_line(ii,jj,ii2,jj,np01::np01_bmp_color(96,96,96,false));
    }

for(i=0; i < node_count_a; ++i){
    node_a = m_node_map_a->get_node_by_idx(i);
    AA_ALWAYS_ASSERT(NULL != node_a);
    xa = node_a->get_x();
    ya = node_a->get_y();
    ii = pixel_border_width +
        static_cast<np01_uint32>(floor((xa-x_shp_min)*scale_factor));
    jj = pixel_border_width +
        static_cast<np01_uint32>(floor((ya-y_shp_min)*scale_factor));
    bmp_file.draw_diamond(ii,jj,7,np01::np01_bmp_color(0,255,0,true));
    }

for(i=0; i < node_count_b; ++i){
    node_b = m_node_map_b->get_node_by_idx(i);
    AA_ALWAYS_ASSERT(NULL != node_b);
    xb = node_b->get_x();
    yb = node_b->get_y();
    ii = pixel_border_width +
        static_cast<np01_uint32>(floor((xb-x_shp_min)*scale_factor));
    jj = pixel_border_width +
        static_cast<np01_uint32>(floor((yb-y_shp_min)*scale_factor));
    bmp_file.draw_box(ii-2,jj-2,ii+2,jj+2,np01::np01_bmp_color(255,0,0,true));
    }

for(i=0; i < edge_count;++i){
    edge=m_edge_map->get_edge_by_idx(i);
    AA_ALWAYS_ASSERT(NULL != edge);
    node_a = edge->get_node_a();
    node_b = edge->get_node_b();
    if((NULL!=node_a) && (NULL!=node_b)){
        xa = node_a->get_x();
        ya = node_a->get_y();
        xb = node_b->get_x();
        yb = node_b->get_y();
        ii = pixel_border_width + static_cast<np01_uint32>(floor((xa-x_shp_min)*scale_factor));
        jj = pixel_border_width + static_cast<np01_uint32>(floor((ya-y_shp_min)*scale_factor));
        ii2 = pixel_border_width + static_cast<np01_uint32>(floor((xb-x_shp_min)*scale_factor));
        jj2 = pixel_border_width + static_cast<np01_uint32>(floor((yb-y_shp_min)*scale_factor));
        bmp_file.draw_line(ii,jj,ii2,jj2,np01::np01_bmp_color(0,0,255,false));
        }
    }

bmp_file.write_file(file_name);
AA_DECR_CALL_DEPTH();
}


/* sort, remove duplicates */
void np01_xy_pair_sub_soln::sort_unique_ab_idx_pairs(){
std::sort(m_ab_idx_pair_vec.begin(), m_ab_idx_pair_vec.end());
np01_uint32_pair_vec_itr unique_end_itr = std::unique(
    m_ab_idx_pair_vec.begin(), m_ab_idx_pair_vec.end());
m_ab_idx_pair_vec.erase(unique_end_itr, m_ab_idx_pair_vec.end());
}

std::ostream& np01_xy_pair_sub_soln::ostream_output(std::ostream& os) const{
os << "<xy_pair_sub_soln=" << std::hex << this << std::dec << ">\n";
size_t i, i_end;
i_end = (m_ab_idx_pair_vec.size() > m_edge_len_vec.size()) ?
    m_ab_idx_pair_vec.size() : m_edge_len_vec.size();
for(i = 0; i < i_end; ++i){
    os << i << ",";
    if(i < m_ab_idx_pair_vec.size()){
        const np01_uint32_pair& ab_idx_pair = m_ab_idx_pair_vec.at(i);
        os << "[" << ab_idx_pair.first << "," << ab_idx_pair.second << "],";
        }
    else {
        os << "[,],";
        }
    if(i < m_edge_len_vec.size()){
        const np01_uint32_pair& ab_idx_pair = m_ab_idx_pair_vec.at(i);
        os << m_edge_len_vec.at(i) << ",\n";
        }
    else {
        os << ",\n";
        }
    }
os << "<cost></unp=" << m_cost.m_unpaired_count
   << "></vsq=" << m_cost.m_viol_edge_len_sq_sum
   << "></lensum=" << m_cost.m_edge_len_sum << ">\n";

os << "</xy_pair_sub_soln>\n";
return os;
}


void np01_xy_pair_genetic_wksp_init_params::set_default(
    np01_xy_pair_map *xy_pair_map){
m_xy_pair_map = xy_pair_map;
m_target_max_edge_len = 0.0;
m_max_iteration_count = 256;
m_no_change_max_iteration_count = 64;
m_inner_loop_max_iteration_count = 32;
}


np01_xy_pair_genetic_wksp::np01_xy_pair_genetic_wksp():
    m_xy_pair_map(NULL),
    m_target_max_edge_len(0.0),
    m_max_iteration_count(0),
    m_no_change_max_iteration_count(0),
    m_inner_loop_max_iteration_count(0),
    m_usage_count_a_vec(),
    m_usage_count_b_vec(),
    m_a_idx(0),
    m_b_idx(0),
    m_node_a_count(0),
    m_node_b_count(0),
    m_a_offset(0),
    m_b_offset(0),
    m_target_usage_count(0),
    m_soln_vec(),
    m_soln_free_chain(NULL),
    m_circle_a_nodes(),
    m_circle_b_nodes(),
    m_near_edges(),
    m_utility_vec(),
    m_unused_idx_set_a(),
    m_unused_idx_set_b(),
    m_changed(false),
    m_iteration_count(0),
    m_no_change_iteration_count(0),
    m_inner_loop_iteration_count(0),
    m_rand_uint32(0)
{}

np01_xy_pair_genetic_wksp::np01_xy_pair_genetic_wksp(
    const np01_xy_pair_genetic_wksp_init_params& init_params):
    m_xy_pair_map(NULL),
    m_target_max_edge_len(0.0),
    m_max_iteration_count(0),
    m_no_change_max_iteration_count(0),
    m_inner_loop_max_iteration_count(0),
    m_usage_count_a_vec(),
    m_usage_count_b_vec(),
    m_a_idx(0),
    m_b_idx(0),
    m_node_a_count(0),
    m_node_b_count(0),
    m_a_offset(0),
    m_b_offset(0),
    m_target_usage_count(0),
    m_circle_a_nodes(),
    m_circle_b_nodes(),
    m_near_edges(),
    m_utility_vec(),
    m_unused_idx_set_a(),
    m_unused_idx_set_b(),
    m_soln_vec(),
    m_soln_free_chain(NULL),
    m_changed(false),
    m_iteration_count(0),
    m_no_change_iteration_count(0),
    m_inner_loop_iteration_count(0),
    m_rand_uint32(0){
init(init_params);
}

np01_xy_pair_genetic_wksp::~np01_xy_pair_genetic_wksp(){
free_soln_vec();
free_free_chain();
}

void np01_xy_pair_genetic_wksp::init(
    const np01_xy_pair_genetic_wksp_init_params& init_params){
m_xy_pair_map = init_params.m_xy_pair_map;
m_target_max_edge_len = init_params.m_target_max_edge_len;
m_max_iteration_count = init_params.m_max_iteration_count;
m_no_change_max_iteration_count = init_params.m_no_change_max_iteration_count;
m_inner_loop_max_iteration_count=init_params.m_inner_loop_max_iteration_count;
m_usage_count_a_vec.clear(); /* count: node A[i] was processed */
m_usage_count_b_vec.clear(); /* count: node B[i] was processed */
m_a_idx = 0;
m_b_idx = 0;
m_node_a_count = 0;
m_node_b_count = 0;
m_a_offset = 0;
m_b_offset = 0;
m_target_usage_count = 0;
m_circle_a_nodes.clear();
m_circle_b_nodes.clear();
m_near_edges.clear();
m_utility_vec.clear();
m_unused_idx_set_a.clear();
m_unused_idx_set_b.clear();
free_soln_vec();
m_changed = false;
m_iteration_count = 0;
m_no_change_iteration_count = 0;
m_inner_loop_iteration_count = 0;
m_rand_uint32 = 4294967291ul;
}

void np01_xy_pair_genetic_wksp::run(){
AA_INCR_CALL_DEPTH();
time_t start_time, now, prev_update_time;
start_time = time(NULL);
prev_update_time = start_time;
const np01_xy_node_map *node_map_a = (NULL==m_xy_pair_map) ?
    NULL : m_xy_pair_map->get_node_map_a();
const np01_xy_node_map *node_map_b = (NULL==m_xy_pair_map) ?
    NULL : m_xy_pair_map->get_node_map_b();
np01_xy_edge_map *edge_map = (NULL==m_xy_pair_map) ? NULL :
    m_xy_pair_map->get_edge_map();
m_node_a_count = (NULL == node_map_a) ? 0 : node_map_a->get_node_count();
m_node_b_count = (NULL == node_map_b) ? 0 : node_map_b->get_node_count();
bool done = (0 == m_node_a_count) || (0 == m_node_b_count);
m_iteration_count = 0;
while(!done){
    m_usage_count_a_vec.clear();
    m_usage_count_a_vec.resize(m_node_a_count, 0);
    m_usage_count_b_vec.clear();
    m_usage_count_b_vec.resize(m_node_b_count, 0);
    m_a_idx=0;
    m_b_idx=0;
    advance_rand();
    m_a_offset = m_rand_uint32 % m_node_a_count;
    advance_rand();
    m_b_offset = m_rand_uint32 % m_node_b_count;
    m_changed = false;

    while((m_a_idx < m_node_a_count) || (m_b_idx < m_node_b_count)){
        advance_a_idx();
        advance_b_idx();
        if((m_a_idx < m_node_a_count) || (m_b_idx < m_node_b_count)){

            /* work on a subset */
            seed_soln_vec();
            for(m_inner_loop_iteration_count = 0; 
                m_inner_loop_iteration_count < m_inner_loop_max_iteration_count;
                ++m_inner_loop_iteration_count){
                mutate();
                gen_crossover();
                cull();
                }
            save_best_soln();
            free_soln_vec();
            }
        }

    ++m_iteration_count;
    if(m_iteration_count >= m_max_iteration_count){
        done = true; }

    if(m_changed){
        m_no_change_iteration_count = 0; }
    else{
        ++m_no_change_iteration_count; }
    if(m_no_change_iteration_count >= m_no_change_max_iteration_count){
        done = true; }

    /* report progress */
    now = time(NULL);
    const double time_elapsed = static_cast<double>(now - start_time);
    if( ((time_elapsed < 10.0) && (now >= (prev_update_time + 1))) ||
        (now >= (prev_update_time + 10))){
        const double progress_ratio = (m_max_iteration_count>0) ?
            static_cast<double>(m_iteration_count)/
            static_cast<double>(m_max_iteration_count) : 1.0;
        const double total_expected_time = (progress_ratio > 0.0) ?
            time_elapsed / progress_ratio : time_elapsed;
        const double expected_remaining_time=total_expected_time-time_elapsed;
        np01_xy_pair_genetic_wksp_result result = {};
        get_result(&result);
        const np01_uint32& crossovers = result.m_total_edge_crossover_count;
        const np01_xy_pair_sub_soln_cost& cost = result.m_solution.get_cost();
        bool cost_valid = (result.m_solution).is_cost_valid();
        AA_ALWAYS_ASSERT(cost_valid);
        const np01_float64 edge_count_f = (NULL == edge_map) ?
            0.0 : edge_map->get_edge_count();
        np01_float64 ave_len = (edge_count_f > 0.0) ?
            cost.m_edge_len_sum/ edge_count_f : 0.0;
        if(prev_update_time == start_time){
            std::cout <<
                "[iteration count/max iteration count]\n"
                "|        progress \n"
                "|        |       time remaining\n"
                "|        |       |            sum of sq of violation lengths\n"
                "|        |       |            |           average length\n"
                "|        |       |            |           |              crossovers\n";
            }
        std::cout
           << "[" <<  m_iteration_count << "/" << m_max_iteration_count << "] "
           << std::fixed << std::setprecision(2) 
           << (progress_ratio * 100.0) << "%  "
           << static_cast<int>(expected_remaining_time) << " s left "
           << std::setprecision(3) 
           << "  viol=" << cost.m_viol_edge_len_sq_sum
           << "  ave_len=" << ave_len
           << "  X=" << crossovers << "\n";
        prev_update_time = now;
        }
    }

AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::get_result(
    np01_xy_pair_genetic_wksp_result *result) const{
AA_INCR_CALL_DEPTH();
if(NULL == result){AA_ALWAYS_ASSERT(false);}
(result->m_solution).reset();
(result->m_edge_crossover_count_vec).clear();
result->m_total_edge_crossover_count = 0;
result->m_len_viol_count = 0;
result->m_max_result_length = 0.0;
if(NULL != m_xy_pair_map){
    np01_xy_edge_map *edge_map = m_xy_pair_map->get_edge_map();
    AA_ALWAYS_ASSERT(NULL != edge_map);
    if( NULL != edge_map){
        const np01_uint32 edge_count = edge_map->get_edge_count();
        for(np01_uint32 i = 0; i < edge_count; ++i ){
            const np01_xy_edge *edge = edge_map->get_edge_by_idx(i);
            AA_ALWAYS_ASSERT(NULL != edge);
            const np01_xy_node *node_a = edge->get_node_a();
            const np01_xy_node *node_b = edge->get_node_b();
            const np01_uint32 idx_a = (NULL == node_a) ?
                np01_xy_node::get_invalid_idx() : node_a->get_idx();
            const np01_uint32 idx_b = (NULL == node_b) ?
                np01_xy_node::get_invalid_idx() : node_b ->get_idx();
            (result->m_solution).add_ab_idx_pair(
                np01_uint32_pair(idx_a, idx_b));
            }            
        (result->m_solution).sort_unique_ab_idx_pairs();
        edge_map->get_edge_crossover_count(
            &(result->m_edge_crossover_count_vec),
            &(result->m_total_edge_crossover_count));
        compute_len_cost(&(result->m_solution));
        for(np01_uint32 i = 0; i < edge_count; ++i ){
            const np01_float64& len = (result->m_solution).get_edge_len(i);
            if(len > result->m_max_result_length){
                result->m_max_result_length = len; }
            if((m_target_max_edge_len > 0.0) && (len > m_target_max_edge_len)){
                ++(result->m_len_viol_count);
                }
            }
        }    
    }

AA_DECR_CALL_DEPTH();
}

/* advance m_a_idx until unused spot is found or m_node_a_count is reached */
void np01_xy_pair_genetic_wksp::advance_a_idx(){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(m_usage_count_a_vec.size() == m_node_a_count);
bool done = (m_a_idx >= m_node_a_count);
while (!done) {
    const np01_uint32 i = (m_a_idx + m_a_offset) % m_node_a_count;
    if(m_usage_count_a_vec.at(i) == 0){
        done = true; } /* unused */
    else if(m_a_idx >= m_node_a_count){
        done = true; } /* reached end */
    else{
        ++m_a_idx; }
    }
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::advance_b_idx(){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(m_usage_count_b_vec.size() == m_node_b_count);
bool done = (m_b_idx >= m_node_b_count);
while (!done) {
    const np01_uint32 i = (m_b_idx + m_b_offset) % m_node_b_count;
    if(m_usage_count_b_vec.at(i) == 0){
        done = true; } /* unused */
    else if(m_b_idx >= m_node_b_count){
        done = true; } /* reached end */
    else{
        ++m_b_idx; }
    }
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::seed_soln_vec(){
AA_INCR_CALL_DEPTH();
np01_xy xy_ctr(0.0, 0.0);
np01_xy_edge *connected_edge = NULL;
if((m_node_a_count - m_a_idx) > (m_node_b_count - m_b_idx)){
    const np01_uint32 ia = (m_a_idx + m_a_offset) % m_node_a_count;
    np01_xy_node *node_a=m_xy_pair_map->get_node_map_a()->get_node_by_idx(ia);
    xy_ctr.first  = node_a->get_x();
    xy_ctr.second = node_a->get_y();
    connected_edge = node_a->get_edge();

    /* make sure a_idx advances by ensuring the usage count is incremented */
    ++m_usage_count_a_vec.at(ia);
    }
else{
    const np01_uint32 ib = (m_b_idx + m_b_offset) % m_node_b_count;
    np01_xy_node *node_b=m_xy_pair_map->get_node_map_b()->get_node_by_idx(ib);
    xy_ctr.first  = node_b->get_x();
    xy_ctr.second = node_b->get_y();
    connected_edge = node_b->get_edge();

    /* make sure b_idx advances by ensuring the usage count is incremented */
    ++m_usage_count_b_vec.at(ib);
    }

advance_rand();
np01_float64 r_ratio = 1.0;
switch(get_rand_int() % 16){
    case 0: r_ratio = 0.1; break;
    case 1: r_ratio = 0.2; break;
    case 2: r_ratio = 0.25; break;
    case 3: r_ratio = 0.3; break;
    case 4: r_ratio = 0.35; break;
    case 5: r_ratio = 0.4; break;
    case 6: r_ratio = 0.45; break;
    case 7: r_ratio = 0.5; break;
    case 8: r_ratio = 0.55; break;
    case 9: r_ratio = 0.6; break;
    case 10: r_ratio = 0.7; break;
    case 11: r_ratio = 0.75; break;
    case 12: r_ratio = 0.8; break;
    case 13: r_ratio = 0.9; break;
    case 14: r_ratio = 1.0; break;
    case 15: r_ratio = 1.5; break;
    default:
        break;
    }

np01_float64 grid_sq_sz =
    (m_xy_pair_map->get_edge_map())->get_loc_grid_dim().get_sq_size();
np01_float64 search_r = r_ratio * grid_sq_sz;

AUTO_ASSERT(m_circle_a_nodes.empty());
AUTO_ASSERT(m_circle_b_nodes.empty());
AUTO_ASSERT(m_near_edges.empty());
AUTO_ASSERT(m_utility_vec.empty());

m_xy_pair_map->get_node_map_a()->get_nodes_in_circle(xy_ctr, search_r,
    &m_circle_a_nodes);
m_xy_pair_map->get_node_map_b()->get_nodes_in_circle(xy_ctr, search_r,
    &m_circle_b_nodes);
advance_rand();
if(((get_rand_int() % 2) == 0) || (NULL==connected_edge)){
    m_xy_pair_map->get_edge_map()->get_edges_in_circle(xy_ctr, search_r,
        &m_near_edges, &m_utility_vec);
    }
else{
    m_xy_pair_map->get_edge_map()->get_intersecting_edges(connected_edge,
        &m_near_edges, &m_utility_vec);
    }

np01_xy_pair_sub_soln *first_soln = alloc_sub_soln();
m_soln_vec.push_back(first_soln);

const np01_xy_node *node_a;
const np01_xy_node *node_b;
const np01_xy_edge *edge;
np01_uint32 node_a_idx;
np01_uint32 node_b_idx;
np01_xy_node_map::node_vec_citr node_a_itr = m_circle_a_nodes.begin();
for(; node_a_itr != m_circle_a_nodes.end(); ++node_a_itr){
    node_a = *node_a_itr;
    AA_ALWAYS_ASSERT(NULL != node_a);
    node_a_idx = node_a->get_idx();
    edge = node_a->get_edge();
    node_b = (NULL == edge) ? NULL : edge->get_node_b();
    node_b_idx = (NULL == node_b) ?
        np01_xy_node::get_invalid_idx() : node_b->get_idx();
    first_soln->add_ab_idx_pair(np01_uint32_pair(node_a_idx, node_b_idx));
    }

np01_xy_node_map::node_vec_citr node_b_itr = m_circle_b_nodes.begin();
for(; node_b_itr != m_circle_b_nodes.end(); ++node_b_itr){
    node_b = *node_b_itr;
    AA_ALWAYS_ASSERT(NULL != node_b);
    node_b_idx = node_b->get_idx();
    edge = node_b->get_edge();
    node_a = (NULL == edge) ? NULL : edge->get_node_a();
    node_a_idx = (NULL == node_a) ?
        np01_xy_node::get_invalid_idx() : node_a->get_idx();
    first_soln->add_ab_idx_pair(np01_uint32_pair(node_a_idx, node_b_idx));
    }

np01_xy_edge_map::edge_vec_citr edge_itr = m_near_edges.begin();
for(; edge_itr != m_near_edges.end(); ++edge_itr){
    edge = *edge_itr;
    AA_ALWAYS_ASSERT(NULL != edge);
    node_a = edge->get_node_a();
    AA_ALWAYS_ASSERT(NULL != node_a);
    node_a_idx = node_a->get_idx();
    node_b = edge->get_node_b();
    AA_ALWAYS_ASSERT(NULL != node_b);
    node_b_idx = node_b->get_idx();
    first_soln->add_ab_idx_pair(np01_uint32_pair(node_a_idx, node_b_idx));
    }

first_soln->sort_unique_ab_idx_pairs();
compute_len_cost(first_soln);

/* mark nodes as used */
const size_t idx_pair_count = first_soln->get_ab_idx_pair_count();
for(size_t ii = 0; ii < idx_pair_count; ++ii ){
    const np01_uint32_pair& idx_pair = first_soln->get_ab_idx_pair(ii);
    node_a_idx = idx_pair.first;
    node_b_idx = idx_pair.second;
    if(np01_xy_node::get_invalid_idx() != node_a_idx){
        ++(m_usage_count_a_vec.at(node_a_idx)); }
    if(np01_xy_node::get_invalid_idx() != node_b_idx){
        ++(m_usage_count_b_vec.at(node_b_idx)); }
    }

m_circle_a_nodes.clear();
m_circle_b_nodes.clear();
m_near_edges.clear();
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::mutate(){
AA_INCR_CALL_DEPTH();
/* pick random existing solution*/
AA_ALWAYS_ASSERT(!m_soln_vec.empty());
advance_rand();
size_t base_soln_idx = m_rand_uint32 % (m_soln_vec.size());
np01_xy_pair_sub_soln *base_soln = m_soln_vec.at(base_soln_idx);
AA_ALWAYS_ASSERT(NULL != base_soln);

size_t ab_idx_pair_count = base_soln->get_ab_idx_pair_count();
if( ab_idx_pair_count > 1){
    /* copy base solution */
    np01_xy_pair_sub_soln *s = alloc_sub_soln();
    m_soln_vec.push_back(s);
    s->reserve(ab_idx_pair_count);
    for(size_t i = 0; i < ab_idx_pair_count; ++i){
        const np01_uint32_pair& ab_idx_pair = base_soln->get_ab_idx_pair(i);
        s->add_ab_idx_pair(ab_idx_pair);
        }

    /* swap around 2, 3 or 4 pairs */
    advance_rand();
    switch(ab_idx_pair_count){
        case 0: break;
        case 1: break;
        case 2:
            mutate2(s);
            break;
        case 3:
            switch(m_rand_uint32 % 2){
                case 0: mutate2(s); break;
                case 1: default: mutate3(s); break;
                }
            break;
        case 4:
        default: 
            switch(m_rand_uint32 % 3){
                case 0: mutate2(s); break;
                case 1: mutate3(s); break;
                case 2: default: mutate4(s); break;
                }
            break;
        }
    compute_len_cost(s);
    }
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::mutate2(np01_xy_pair_sub_soln *s){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != s);

/* pick two random indices */
size_t ab_idx_pair_count = s->get_ab_idx_pair_count();
AA_ALWAYS_ASSERT(ab_idx_pair_count >= 2);
advance_rand();
size_t ii = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t jj = m_rand_uint32 % ab_idx_pair_count;
if(ii == jj){ jj = (ii+1) % ab_idx_pair_count; }
AA_ALWAYS_ASSERT(ii!=jj);

/* swap B indices */
np01_uint32_pair p_ii = s->get_ab_idx_pair(ii);
np01_uint32_pair p_jj = s->get_ab_idx_pair(jj);
np01_uint32 temp = p_ii.second;
p_ii.second = p_jj.second;
p_jj.second = temp;
s->set_ab_idx_pair(ii, p_ii);
s->set_ab_idx_pair(jj, p_jj);

AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::mutate3(np01_xy_pair_sub_soln *s){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != s);

/* pick three random indices */
size_t ab_idx_pair_count = s->get_ab_idx_pair_count();
AA_ALWAYS_ASSERT(ab_idx_pair_count >= 3);
advance_rand();
size_t ii = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t jj = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t kk = m_rand_uint32 % ab_idx_pair_count;
if((ii == jj) || (jj == kk) || (ii == kk)){
    jj = (ii+1) % ab_idx_pair_count; 
    kk = (ii+2) % ab_idx_pair_count; 
    }
AA_ALWAYS_ASSERT(ii!=jj);
AA_ALWAYS_ASSERT(ii!=kk);
AA_ALWAYS_ASSERT(jj!=kk);

/* swap B indices */
np01_uint32_pair p_ii = s->get_ab_idx_pair(ii);
np01_uint32_pair p_jj = s->get_ab_idx_pair(jj);
np01_uint32_pair p_kk = s->get_ab_idx_pair(kk);
np01_uint32 temp = p_ii.second;
p_ii.second = p_jj.second;
p_jj.second = p_kk.second;
p_kk.second = temp;
s->set_ab_idx_pair(ii, p_ii);
s->set_ab_idx_pair(jj, p_jj);
s->set_ab_idx_pair(kk, p_kk);

AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::mutate4(np01_xy_pair_sub_soln *s){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != s);

/* pick four random indices */
size_t ab_idx_pair_count = s->get_ab_idx_pair_count();
AA_ALWAYS_ASSERT(ab_idx_pair_count >= 4);
advance_rand();
size_t ii = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t jj = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t kk = m_rand_uint32 % ab_idx_pair_count;
advance_rand();
size_t ll = m_rand_uint32 % ab_idx_pair_count;
if((ii==jj) || (ii==kk) || (ii==ll) || (jj==kk) || (jj==ll) || (kk==ll)){
    jj = (ii+1) % ab_idx_pair_count; 
    kk = (ii+2) % ab_idx_pair_count;  
    ll = (ii+3) % ab_idx_pair_count; 
    }
AA_ALWAYS_ASSERT(ii!=jj);
AA_ALWAYS_ASSERT(ii!=kk);
AA_ALWAYS_ASSERT(ii!=ll);
AA_ALWAYS_ASSERT(jj!=kk);
AA_ALWAYS_ASSERT(jj!=ll);
AA_ALWAYS_ASSERT(kk!=ll);

/* swap B indices */
np01_uint32_pair p_ii = s->get_ab_idx_pair(ii);
np01_uint32_pair p_jj = s->get_ab_idx_pair(jj);
np01_uint32_pair p_kk = s->get_ab_idx_pair(kk);
np01_uint32_pair p_ll = s->get_ab_idx_pair(ll);
np01_uint32 temp = p_ii.second;
p_ii.second = p_jj.second;
p_jj.second = p_kk.second;
p_kk.second = p_ll.second;
p_ll.second = temp;
s->set_ab_idx_pair(ii, p_ii);
s->set_ab_idx_pair(jj, p_jj);
s->set_ab_idx_pair(kk, p_kk);
s->set_ab_idx_pair(ll, p_ll);

AA_DECR_CALL_DEPTH();
}

/* perform genetic crossovers only occasionally */
void np01_xy_pair_genetic_wksp::gen_crossover(){
size_t i=0;
AA_INCR_CALL_DEPTH();
/* perform genetic crossovers only occasionally */
advance_rand();
if( m_soln_vec.size() >= (NP01_GEN_MAX_POPULATION_SIZE/2) &&
    ((m_rand_uint32 % 32) < 4)){
    /* pick two existing solutions as parents */
    AA_ALWAYS_ASSERT(m_soln_vec.size() >= 2);
    advance_rand();
    size_t parent_idx_0 = get_rand_int() % (m_soln_vec.size());
    advance_rand();
    size_t parent_idx_1 = get_rand_int() % (m_soln_vec.size());
    if(parent_idx_0 == parent_idx_1){
        parent_idx_1 = (parent_idx_0 + 1) % (m_soln_vec.size());
        }
    AA_ALWAYS_ASSERT(parent_idx_0 != parent_idx_1);
    const np01_xy_pair_sub_soln *parent_0 = m_soln_vec.at(parent_idx_0);
    const np01_xy_pair_sub_soln *parent_1 = m_soln_vec.at(parent_idx_1);
    AA_ALWAYS_ASSERT(NULL != parent_0);
    AA_ALWAYS_ASSERT(NULL != parent_1);

    /* create empty child solution */
    const size_t ab_idx_pair_count = parent_0->get_ab_idx_pair_count();
    AUTO_ASSERT(parent_1->get_ab_idx_pair_count() == ab_idx_pair_count);
    np01_xy_pair_sub_soln *s = alloc_sub_soln();
    m_soln_vec.push_back(s);
    s->reserve(ab_idx_pair_count);

    /* initialize child solution by copying elements that parents
    have in common */
    AUTO_ASSERT(m_unused_idx_set_a.empty());
    AUTO_ASSERT(m_unused_idx_set_b.empty());
    for(i = 0; i < ab_idx_pair_count; ++i){
        const np01_uint32_pair& ab_idx_pair_0 = parent_0->get_ab_idx_pair(i);
        const np01_uint32_pair& ab_idx_pair_1 = parent_1->get_ab_idx_pair(i);
        if(ab_idx_pair_0 == ab_idx_pair_1){
            AUTO_ASSERT(m_unused_idx_set_a.find(ab_idx_pair_0.first)
                    ==m_unused_idx_set_a.end());
            AUTO_ASSERT(m_unused_idx_set_b.find(ab_idx_pair_0.second)
                    ==m_unused_idx_set_b.end());
            s->add_ab_idx_pair(ab_idx_pair_0);
            }
        else{
            /* save unused indices */
            if(np01_xy_node::get_invalid_idx() != ab_idx_pair_0.first){ 
                m_unused_idx_set_a.insert(ab_idx_pair_0.first);}
            if(np01_xy_node::get_invalid_idx() != ab_idx_pair_1.first){ 
                m_unused_idx_set_a.insert(ab_idx_pair_1.first);}
            if(np01_xy_node::get_invalid_idx() != ab_idx_pair_0.second){ 
                m_unused_idx_set_b.insert(ab_idx_pair_0.second);}
            if(np01_xy_node::get_invalid_idx() != ab_idx_pair_1.second){ 
                m_unused_idx_set_b.insert(ab_idx_pair_1.second);}

            /* store dummy value */
            s->add_ab_idx_pair(np01_uint32_pair(
                np01_xy_node::get_invalid_idx(),
                np01_xy_node::get_invalid_idx()));
            }        
        }
    AUTO_ASSERT(s->get_ab_idx_pair_count() == ab_idx_pair_count);

    /* initialize child solution by copying elements from one or
    the other parent. */
    for(i = 0; i < ab_idx_pair_count; ++i){
        const np01_uint32_pair& ab_idx_pair_0 = parent_0->get_ab_idx_pair(i);
        const np01_uint32_pair& ab_idx_pair_1 = parent_1->get_ab_idx_pair(i);
        if(ab_idx_pair_0 == ab_idx_pair_1){
            /* parent data matches */
            AUTO_ASSERT(m_unused_idx_set_a.find(ab_idx_pair_0.first)
                    ==m_unused_idx_set_a.end());
            AUTO_ASSERT(m_unused_idx_set_b.find(ab_idx_pair_0.second)
                    ==m_unused_idx_set_b.end());
            }
        else{
            /* pick data from one or the other parent, choose randomly */
            advance_rand();
            np01_uint32_pair ab_idx_pair = ((get_rand_int() % 2) == 0) ?
                ab_idx_pair_0 : ab_idx_pair_1;

            /* modify as needed to avoid using any index more than once.*/
            if(np01_xy_node::get_invalid_idx() != ab_idx_pair.first){
                if(m_unused_idx_set_a.empty()){
                    /* no index A available*/
                    ab_idx_pair.first = np01_xy_node::get_invalid_idx();
                    }
                else{
                    np01_uint32_set_itr itr_a =
                        m_unused_idx_set_a.find(ab_idx_pair.first);
                    if(itr_a == m_unused_idx_set_a.end()){
                        /* idx_a already used.  pick a different index */
                        itr_a = m_unused_idx_set_a.begin();
                        AA_ALWAYS_ASSERT(itr_a != m_unused_idx_set_a.end());
                        ab_idx_pair.first = *itr_a;
                        }
                    AUTO_ASSERT(*itr_a == ab_idx_pair.first);
                    m_unused_idx_set_a.erase(itr_a);
                    }
                }

            if(np01_xy_node::get_invalid_idx() != ab_idx_pair.second){
                if(m_unused_idx_set_b.empty()){
                    /* no index B available*/
                    ab_idx_pair.second = np01_xy_node::get_invalid_idx();
                    }
                else{
                    np01_uint32_set_itr itr_b =
                        m_unused_idx_set_b.find(ab_idx_pair.second);
                    if(itr_b == m_unused_idx_set_b.end()){
                        /* idx_b already used.  pick a different index */
                        itr_b = m_unused_idx_set_b.begin();
                        AA_ALWAYS_ASSERT(itr_b != m_unused_idx_set_b.end());
                        ab_idx_pair.second = *itr_b;
                        }
                    AUTO_ASSERT(*itr_b == ab_idx_pair.second);
                    m_unused_idx_set_b.erase(itr_b);
                    }
                }

            /* add index pair to child solution */
            s->set_ab_idx_pair(i, ab_idx_pair);
            }
        }

    /* use any unused indices by replacing invalid indices */
    for(i = 0; (i < ab_idx_pair_count) && 
       ( (!m_unused_idx_set_a.empty()) ||
         (!m_unused_idx_set_b.empty()) ); ++i){
        np01_uint32_pair ab_idx_pair = s->get_ab_idx_pair(i);
        if((np01_xy_node::get_invalid_idx() == ab_idx_pair.first) ||
           (np01_xy_node::get_invalid_idx() == ab_idx_pair.second)){
            if((np01_xy_node::get_invalid_idx() == ab_idx_pair.first)&&
                (!m_unused_idx_set_a.empty())){
                ab_idx_pair.first = *(m_unused_idx_set_a.begin());
                m_unused_idx_set_a.erase(m_unused_idx_set_a.begin());
                }
            if((np01_xy_node::get_invalid_idx() == ab_idx_pair.second)&&
                (!m_unused_idx_set_b.empty())){
                ab_idx_pair.second = *(m_unused_idx_set_b.begin());
                m_unused_idx_set_b.erase(m_unused_idx_set_b.begin());
                }
            s->set_ab_idx_pair(i,ab_idx_pair);
            }
        }

    s->sort_unique_ab_idx_pairs();
    compute_len_cost(s);

    AUTO_ASSERT(s->get_ab_idx_pair_count() == ab_idx_pair_count);
    AUTO_ASSERT(m_unused_idx_set_a.empty());
    AUTO_ASSERT(m_unused_idx_set_b.empty());


    if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_2)){
        np01_uint32_set parent_0_a_idx_set;
        np01_uint32_set parent_0_b_idx_set;
        np01_uint32_set parent_1_a_idx_set;
        np01_uint32_set parent_1_b_idx_set;
        np01_uint32_set s_a_idx_set;
        np01_uint32_set s_b_idx_set;
        for(size_t ii = 0; ii < ab_idx_pair_count; ++ii){
            const np01_uint32_pair& ab_idx_pair_0=parent_0->get_ab_idx_pair(ii);
            const np01_uint32_pair& ab_idx_pair_1=parent_1->get_ab_idx_pair(ii);
            const np01_uint32_pair& ab_idx_pair_s=s->get_ab_idx_pair(ii);
            parent_0_a_idx_set.insert(ab_idx_pair_0.first);
            parent_0_b_idx_set.insert(ab_idx_pair_0.second);
            parent_1_a_idx_set.insert(ab_idx_pair_1.first);  
            parent_1_b_idx_set.insert(ab_idx_pair_1.second); 
            s_a_idx_set.insert(ab_idx_pair_s.first);  
            s_b_idx_set.insert(ab_idx_pair_s.second); 
            }
        AA_XDBG_ASSERT(parent_0_a_idx_set == parent_1_a_idx_set,
            CF01_AA_DEBUG_LEVEL_2);
        AA_XDBG_ASSERT(parent_1_a_idx_set==s_a_idx_set, CF01_AA_DEBUG_LEVEL_2);
        AA_XDBG_ASSERT(s_a_idx_set==parent_0_a_idx_set, CF01_AA_DEBUG_LEVEL_2);
        AA_XDBG_ASSERT(parent_0_b_idx_set == parent_1_b_idx_set,
            CF01_AA_DEBUG_LEVEL_2);
        AA_XDBG_ASSERT(parent_1_b_idx_set==s_b_idx_set, CF01_AA_DEBUG_LEVEL_2);
        AA_XDBG_ASSERT(s_b_idx_set==parent_0_b_idx_set, CF01_AA_DEBUG_LEVEL_2);
        }
    }
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::cull(){
AA_INCR_CALL_DEPTH();
while(m_soln_vec.size() > NP01_GEN_MAX_POPULATION_SIZE){
    /* find highest cost solution */
    soln_ptr_vec_itr cull_itr = m_soln_vec.begin();
    np01_xy_pair_sub_soln *cull_soln = *cull_itr;
    AA_ALWAYS_ASSERT(NULL != cull_soln);
    AUTO_ASSERT(cull_soln->is_cost_valid());
    soln_ptr_vec_itr search_itr = m_soln_vec.begin();
    for(++search_itr; search_itr != m_soln_vec.end(); ++search_itr){
        const np01_xy_pair_sub_soln *s = *search_itr;
        AA_ALWAYS_ASSERT(NULL != s);
        AUTO_ASSERT(s->is_cost_valid());
        if(cull_soln->get_cost() < s->get_cost()){
            /* found solution with higher cost*/
            cull_itr = search_itr;
            cull_soln = *cull_itr;
            }
        }
    /* remove highest cost solution */
    free_sub_soln(cull_soln);
    m_soln_vec.erase(cull_itr);
    }
AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::save_best_soln(){
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != m_xy_pair_map);
AA_XDBG_ASSERT(0 == m_xy_pair_map->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AA_ALWAYS_ASSERT(!m_soln_vec.empty());

/* find best solution*/
soln_ptr_vec_itr best_itr = m_soln_vec.begin();
AA_ALWAYS_ASSERT(NULL != *best_itr);
AUTO_ASSERT((*best_itr)->is_cost_valid());
soln_ptr_vec_itr search_itr = m_soln_vec.begin();
for(++search_itr; search_itr != m_soln_vec.end(); ++search_itr){
    np01_xy_pair_sub_soln *s = *search_itr;
    AA_ALWAYS_ASSERT(NULL != s);
    AUTO_ASSERT(s->is_cost_valid());
    if(s->get_cost() < (*best_itr)->get_cost()){
        /* found solution with lower cost*/
        best_itr = search_itr;
        }
    }
const np01_xy_pair_sub_soln *best_soln = *best_itr;

/* copy back to main data set */
size_t i;
np01_xy_node *node_a, *node_b;
np01_xy_node *node_a_old_b, *node_b_old_a;
np01_xy_edge *node_a_edge, *node_b_edge;

np01_xy_node_map *node_map_a = m_xy_pair_map->get_node_map_a();
np01_xy_node_map *node_map_b = m_xy_pair_map->get_node_map_b();
np01_xy_edge_map *edge_map = m_xy_pair_map->get_edge_map();
AA_ALWAYS_ASSERT(NULL != node_map_a);
AA_ALWAYS_ASSERT(NULL != node_map_b);
AA_ALWAYS_ASSERT(NULL != edge_map);
const size_t ab_idx_pair_count = best_soln->get_ab_idx_pair_count();
np01_xy_group_type xy_small_group_type = 
    (node_map_a->get_node_count() < node_map_b->get_node_count() ) ?
    NP01_XY_GROUP_TYPE_A : NP01_XY_GROUP_TYPE_B;
bool done = false;

/* disconnect edges */
for(i = 0; i < ab_idx_pair_count; ++i){
    const np01_uint32_pair& ab_idx_pair = best_soln->get_ab_idx_pair(i);
    node_a=node_map_a->get_node_by_idx(ab_idx_pair.first);
    node_a_edge=(NULL == node_a)?NULL:node_a->get_edge();
    node_b=node_map_b->get_node_by_idx(ab_idx_pair.second);
    node_b_edge=(NULL == node_b)?NULL:node_b->get_edge();
    switch(xy_small_group_type){
        case NP01_XY_GROUP_TYPE_A:
            /* disconnect existing B node from edge. leave A with edge */
            if((NULL != node_a) && (NULL != node_b) &&
               (NULL != node_a_edge) && (node_a_edge != node_b_edge)){
                node_a_old_b = node_a_edge->get_node_b();
                AA_ALWAYS_ASSERT(NULL != node_a_old_b);
                AUTO_ASSERT(node_b != node_a_old_b);
                node_a_old_b->disconnect_edge();
                edge_map->remove_edge_from_loc_grid(node_a_edge);
                m_changed = true;
                }
            break;
        case NP01_XY_GROUP_TYPE_B:
            /* disconnect existing A node from edge. leave B with edge */
            if((NULL != node_a) && (NULL != node_b) &&
               (NULL != node_b_edge) && (node_a_edge != node_b_edge)){
                node_b_old_a = node_b_edge->get_node_a();
                AA_ALWAYS_ASSERT(NULL != node_b_old_a);
                AUTO_ASSERT(node_a != node_b_old_a);
                node_b_old_a->disconnect_edge();
                edge_map->remove_edge_from_loc_grid(node_b_edge);
                m_changed = true;
                }
            break;
        case NP01_XY_GROUP_TYPE_UNKNOWN:
        default:
            AA_ALWAYS_ASSERT(false);
            break;
        }
    }

/* connect edges */
for(i = 0; i < ab_idx_pair_count; ++i){
    const np01_uint32_pair& ab_idx_pair = best_soln->get_ab_idx_pair(i);
    node_a=node_map_a->get_node_by_idx(ab_idx_pair.first);
    node_a_edge=(NULL == node_a)?NULL:node_a->get_edge();
    node_b=node_map_b->get_node_by_idx(ab_idx_pair.second);
    node_b_edge=(NULL == node_b)?NULL:node_b->get_edge();
    switch(xy_small_group_type){
        case NP01_XY_GROUP_TYPE_A:
            /* Connect B to A's edge */
            if((NULL != node_a) && (NULL != node_b) &&
               (NULL != node_a_edge) && (NULL == node_b_edge)){
                edge_map->update_edge_ab(node_a_edge, node_a, node_b); 
                AUTO_ASSERT(0 == node_a_edge->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
                AUTO_ASSERT(m_changed);
                }
            break;
        case NP01_XY_GROUP_TYPE_B:
            /* Connect A to B's edge */
            if((NULL != node_a) && (NULL != node_b) &&
               (NULL != node_b_edge) && (NULL == node_a_edge)){
                edge_map->update_edge_ab(node_b_edge, node_a, node_b); 
                AUTO_ASSERT(0 == node_b_edge->verify_data(AA_ERR_BUF(),
                    AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
                AUTO_ASSERT(m_changed);
                }
            break;
        case NP01_XY_GROUP_TYPE_UNKNOWN:
        default:
            AA_ALWAYS_ASSERT(false);
            break;
        }
    }

if(AA_SHOULD_RUN_XDBG(CF01_AA_DEBUG_LEVEL_1)){
    for(i = 0; i < ab_idx_pair_count; ++i){
        const np01_uint32_pair& ab_idx_pair = best_soln->get_ab_idx_pair(i);
        node_a=node_map_a->get_node_by_idx(ab_idx_pair.first);
        node_a_edge=(NULL == node_a)?NULL:node_a->get_edge();
        node_b=node_map_b->get_node_by_idx(ab_idx_pair.second);
        node_b_edge=(NULL == node_b)?NULL:node_b->get_edge();
        if((NULL != node_a) && (NULL != node_b)){
            AA_XDBG_ASSERT(node_a_edge != NULL, CF01_AA_DEBUG_LEVEL_1);
            AA_XDBG_ASSERT(node_b_edge != NULL, CF01_AA_DEBUG_LEVEL_1);
            AA_XDBG_ASSERT(node_a_edge == node_a_edge, CF01_AA_DEBUG_LEVEL_1);
            }
        }
    }

AA_XDBG_ASSERT(0 == m_xy_pair_map->verify_data(AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(),
    AA_ERR_BUF_POS_PTR()), CF01_AA_DEBUG_LEVEL_3);
AUTO_ASSERT(0 == edge_map->verify_edges_fully_connected(
    AA_ERR_BUF(),AA_ERR_BUF_CAPACITY(), AA_ERR_BUF_POS_PTR()));
AA_DECR_CALL_DEPTH();
}

np01_xy_pair_sub_soln *np01_xy_pair_genetic_wksp::alloc_sub_soln(){
AA_INCR_CALL_DEPTH();
np01_xy_pair_sub_soln *s = m_soln_free_chain;
if(NULL == s){
    s = new np01_xy_pair_sub_soln();
    }
else{
    m_soln_free_chain = s->get_free_chain_next();
    s->reset();
    }
AUTO_ASSERT(s->get_free_chain_next() == NULL);
AA_DECR_CALL_DEPTH();
return s;
}

void np01_xy_pair_genetic_wksp::free_sub_soln(np01_xy_pair_sub_soln *s){
s->set_free_chain_next(m_soln_free_chain);
m_soln_free_chain = s;
}

void np01_xy_pair_genetic_wksp::free_soln_vec(){
soln_ptr_vec_citr s_itr = m_soln_vec.begin();
for(; s_itr != m_soln_vec.end(); ++s_itr){
    np01_xy_pair_sub_soln *s = *s_itr;
    free_sub_soln(s);
    }
m_soln_vec.clear();
}

void np01_xy_pair_genetic_wksp::free_free_chain(){
while(NULL != m_soln_free_chain){
    np01_xy_pair_sub_soln *next = m_soln_free_chain->get_free_chain_next();
    delete m_soln_free_chain;
    m_soln_free_chain = next;
    }
}

/* compute lengths and cost */
void np01_xy_pair_genetic_wksp::compute_len_cost(np01_xy_pair_sub_soln *s) const{
AA_INCR_CALL_DEPTH();
AA_ALWAYS_ASSERT(NULL != s);
AA_ALWAYS_ASSERT(NULL != m_xy_pair_map);
np01_xy_pair_sub_soln_cost cost;

np01_xy_node_map *node_map_a = m_xy_pair_map->get_node_map_a();
np01_xy_node_map *node_map_b = m_xy_pair_map->get_node_map_b();
AA_ALWAYS_ASSERT(NULL != node_map_a);
AA_ALWAYS_ASSERT(NULL != node_map_b);
np01_xy_group_type xy_small_group_type = 
    (node_map_a->get_node_count() < node_map_b->get_node_count() ) ?
    NP01_XY_GROUP_TYPE_A : NP01_XY_GROUP_TYPE_B;

cost.m_unpaired_count = 0; /* unpaired A+B. zero if only A or only B */
cost.m_viol_edge_len_sq_sum = 0.0; /* sum(len^2) for violating edges */
cost.m_edge_len_sum = 0.0; /* sum of edge lengths */
size_t i;
const size_t idx_pair_count = s->get_ab_idx_pair_count();
for(i = 0; i < idx_pair_count; ++i){
    const np01_uint32_pair& ab_idx_pair = s->get_ab_idx_pair(i);
    np01_float64 edge_len = 0.0;
    const np01_xy_node *node_a = 
        m_xy_pair_map->get_node_map_a()->get_node_by_idx(ab_idx_pair.first);
    const np01_xy_node *node_b = 
        m_xy_pair_map->get_node_map_b()->get_node_by_idx(ab_idx_pair.second);
    if(NULL != node_a){
        if(NULL != node_b){
            /* A not NULL, B not NULL*/
            const np01_float64 dx = (node_b->get_x()) - (node_a->get_x());
            const np01_float64 dy = (node_b->get_y()) - (node_a->get_y());
            const np01_float64 edge_len_sq = (dx*dx)+(dy*dy);
            edge_len = sqrt(edge_len_sq);
            if((m_target_max_edge_len > 0.0) && 
                (edge_len > m_target_max_edge_len)){
                cost.m_viol_edge_len_sq_sum += edge_len_sq;
                }
            cost.m_edge_len_sum += edge_len;
            }
        else if(NP01_XY_GROUP_TYPE_A == xy_small_group_type){
            /* A not NULL, B NULL*/
            AA_ALWAYS_ASSERT(NULL != node_a);
            AA_ALWAYS_ASSERT(NULL == node_b);
            ++cost.m_unpaired_count;
            }
        }
    else if(NULL != node_b && (NP01_XY_GROUP_TYPE_B == xy_small_group_type)){
        /* A NULL, B not NULL*/
        AA_ALWAYS_ASSERT(NULL == node_a);
        AA_ALWAYS_ASSERT(NULL != node_b);
        ++cost.m_unpaired_count;
        }
    s->set_len(static_cast<np01_uint32>(i),edge_len);
    }
s->set_cost(cost);
s->set_cost_valid(true);

AA_DECR_CALL_DEPTH();
}

void np01_xy_pair_genetic_wksp::advance_rand(){
/* advance random number, pg 359, The Standard C Library (c)1991,
P. J. Plauger, ISBN 0-13-131509-9 */
m_rand_uint32= (m_rand_uint32 * 1103515245) + 12345;
}


const char *np01_main::m_prog_name = "Newport Algorithm 1";
const char *np01_main::m_version_str = "0.1.1 " __DATE__ " " __TIME__;

int np01_main::run(int argc, char *argv[]){
np01_main main(argc, argv);
return main.execute();
}

np01_main::np01_main(int argc, char *argv[]):m_argc(argc),
    m_argv(const_cast<const char **>(&(argv[0]))),
    m_test_option(false),
    m_test_number(0),
    m_pair_points_option(true), /*default algorithm*/
    m_xy_a_file_option(false),
    m_xy_a_filename(),
    m_xy_b_file_option(false),
    m_xy_b_filename(),
    m_out_ab_file_option(false),
    m_out_ab_filename(),
    m_max_len_option(false),
    m_max_len(0.0),
    m_iterations_option(false),
    m_iterations(0),
    m_test_iterations_option(false),
    m_test_iterations(0),
    m_test_rand_seed_option(false),
    m_test_rand_seed(0)
{
parse_cmd_line();
}

np01_main::~np01_main()
{}

int np01_main::execute(){
int error_code = 0;
if(m_test_option){
    std::cout << m_prog_name << " v" << m_version_str << "\n";
    switch( m_test_number){
        case 0: execute_test_0(); break;
        case 1: execute_test_1(); break;
        case 2: execute_test_2(); break;
        case 3: execute_test_3(); break;
        case 4: execute_test_4(); break;
        case 5: execute_test_5(); break;
        default:
            std::cout << "attempt to run test " << m_test_number 
                << " fails.  no such test\n";
            break;
        }
    }
else if(m_pair_points_option){
    std::cout << m_prog_name << " v" << m_version_str << "\n";
    execute_pp();
    }

return error_code;
}

void np01_main::parse_cmd_line(){
/* reset options */
m_test_option = false;
m_test_number = 0;
m_pair_points_option = true; /*default algorithm*/
m_xy_a_file_option = false;
m_xy_a_filename.erase();
m_xy_b_file_option = false;
m_xy_b_filename.erase();
m_out_ab_file_option = false;
m_out_ab_filename.erase();
m_max_len_option = false;
m_max_len = 0.0;
m_iterations_option = false;
m_iterations = 0;
m_test_iterations_option = false;
m_test_iterations = 0;
m_test_rand_seed_option = false;
m_test_rand_seed = 0;

int i = 1;
bool should_print_help = (m_argc < 2);
while (i < m_argc ){
    if(strncmp(m_argv[i], "--help", strlen("--help"))==0){
        m_test_option=false;
        m_pair_points_option=false;
        should_print_help = true;
        ++i;
        }
    else if(strncmp(m_argv[i], "--test=", strlen("--test="))==0){
        m_test_option=true;
        m_pair_points_option=false;
        m_test_number = atoi( (m_argv[i]+strlen("--test=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "-t=", strlen("-t="))==0){
        m_test_option=true;
        m_pair_points_option=false;
        m_test_number = atoi( (m_argv[i]+strlen("-t=")) );
        ++i;
        }
    else if(strcmp(m_argv[i], "--pair-points")==0){
        m_pair_points_option=true;
        m_test_option=false;
        ++i;
        }
    else if(strcmp(m_argv[i], "-pp")==0){
        m_pair_points_option=true;
        ++i;
        }
    else if(strncmp(m_argv[i], "--xy-a-file=", strlen("--xy-a-file="))==0){
        m_xy_a_file_option=true;
        m_xy_a_filename.assign( (m_argv[i]+strlen("--xy-a-file=")) );
        ++i;
        if(m_xy_a_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "-xyaf=", strlen("-xyaf="))==0){
        m_xy_a_file_option=true;
        m_xy_a_filename.assign( (m_argv[i]+strlen("-xyaf=")) );
        ++i;
        if(m_xy_a_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "--xy-b-file=", strlen("--xy-b-file="))==0){
        m_xy_b_file_option=true;
        m_xy_b_filename.assign( (m_argv[i]+strlen("--xy-b-file=")) );
        ++i;
        if(m_xy_b_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "-xybf=", strlen("-xybf="))==0){
        m_xy_b_file_option=true;
        m_xy_b_filename.assign( (m_argv[i]+strlen("-xybf=")) );
        ++i;
        if(m_xy_b_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "--out-ab-file=", strlen("--out-ab-file="))==0){
        m_out_ab_file_option=true;
        m_out_ab_filename.assign( (m_argv[i]+strlen("--out-ab-file=")) );
        ++i;
        if(m_out_ab_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "--oabf=", strlen("--oabf="))==0){
        m_out_ab_file_option=true;
        m_out_ab_filename.assign( (m_argv[i]+strlen("--oabf=")) );
        ++i;
        if(m_out_ab_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "--out-ab-file=", strlen("--out-ab-file="))==0){
        m_out_ab_file_option=true;
        m_out_ab_filename.assign( (m_argv[i]+strlen("--out-ab-file=")) );
        ++i;
        if(m_out_ab_filename.empty()){ should_print_help = true; }
        }
    else if(strncmp(m_argv[i], "--max-len=", strlen("--max-len="))==0){
        m_max_len_option=true;
        m_max_len = atof( (m_argv[i]+strlen("--max-len=")) );
        if(m_max_len < 0.0){
            m_max_len = 0.0; }
        ++i;
        }
    else if(strncmp(m_argv[i], "--iterations=", strlen("--iterations="))==0){
        m_iterations_option=true;
        m_iterations = atoi( (m_argv[i]+strlen("--iterations=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--iter=", strlen("--iter="))==0){
        m_iterations_option=true;
        m_iterations = atoi( (m_argv[i]+strlen("--iter=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--test-iterations=", strlen("--test-iterations="))==0){
        m_test_iterations_option=true;
        m_test_iterations = atoi( (m_argv[i]+strlen("--test-iterations=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--ti=", strlen("--ti="))==0){
        m_test_iterations_option=true;
        m_test_iterations = atoi( (m_argv[i]+strlen("--ti=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--test-rand-seed=", strlen("--test-rand-seed="))==0){
        m_test_rand_seed_option=true;
        m_test_rand_seed = atoi( (m_argv[i]+strlen("--test-rand-seed=")) );
        ++i;
        }
    else if(strncmp(m_argv[i], "--seed=", strlen("--seed="))==0){
        m_test_rand_seed_option=true;
        m_test_rand_seed = atoi( (m_argv[i]+strlen("--seed=")) );
        ++i;
        }
    else {
        std::cout << "unknown option: " << m_argv[i] << "\n";
        should_print_help = true;
        ++i;
        }
    }
if(should_print_help){
    std::cout << m_prog_name << " v" << m_version_str << "\n\n";
    std::cout << 
        "This program connects two groups of points, A and B, which may be of different\n"
        "size.  Each point in the smaller group is assigned to a unique point in the\n"
        "larger group.  The assignment is repeatedly optimized towards two criteria.\n"
        "One, minimize average length, thus minimizing crossovers.  Two, try to keep\n"
        "each length below an optional maximum value.\n\n"
        "There are two input files and two output files.  Each input file has two\n"
        "columns: x and y.  The first output file is a csv file with two columns: A\n"
        "and B.  Each column has the zero-based index of the x,y point in the\n"
        "corresponding input file.  The second output file is a bitmap image of the\n"
        "assignment.    Author: Mac Stevens   stevensm@earthlink.net\n\n";
    std::cout << "usage: " << m_argv[0] << " --pair-points \\\n"
        "       --xy-a-file=<input x,y csv file A> \\           default=a.csv\n"
        "       --xy-b-file=<input x,y csv file B> \\           default=b.csv\n"
        "       --out-ab-file=<output (a idx,b idx) csv file>\\ default=ab.csv\n"
        "       --max-len=<target max length>\\                 default=0.0 =>no max\n"
        "       --iterations=<optimization iteration count>    default=256\n\n";
    std::cout << "Regression Tests\n";
    std::cout << "usage: " << m_argv[0] << " --test=<test_number>\n";
    std::cout << "usage: " << m_argv[0] << " --test=0 \\\n"
        "                  --test-rand-seed=<random num gen seed> \\\n"
        "                  --test-iterations=<test iteration count>\n";
    m_test_option=false;
    m_pair_points_option=false;
    }

/* post-process options*/
if(m_test_option){
    }
else if (m_pair_points_option){
    if(m_xy_a_file_option || m_xy_a_filename.empty()){
        if(m_xy_a_filename.empty()){
            m_xy_a_filename="a.csv"; }
        std::ifstream temp_is(m_xy_a_filename.c_str());
        if(!temp_is.good()){
            std::cout << "can't open xy-a-file " << m_xy_a_filename << "\n"; }
        temp_is.close();
        }
    if(m_xy_b_file_option || m_xy_b_filename.empty()){
        if(m_xy_b_filename.empty()){
            m_xy_b_filename="b.csv"; }
        std::ifstream temp_is(m_xy_b_filename.c_str());
        if(!temp_is.good()){
            std::cout << "can't open xy-b-file " << m_xy_b_filename << "\n"; }
        temp_is.close();
        }
    if(m_out_ab_file_option || m_out_ab_filename.empty()){
        if(m_out_ab_filename.empty()){
            m_out_ab_filename="ab.csv"; }
        std::ofstream temp_os(m_out_ab_filename.c_str(),
            std::ofstream::out | std::ofstream::app );
        if(!temp_os.good()){
            std::cout<<"can't open out-ab-file "<<m_out_ab_filename<<"\n"; }
        temp_os.close();
        }
    }
}

void np01_main::execute_pp(){
const time_t start_time = time(NULL);
time_t finish_time;
std::string line;
std::string x_num_str;
std::string y_num_str;
size_t x_num_pos;
size_t x_end_delim_pos;
size_t y_num_pos;
size_t y_end_delim_pos;
np01::np01_xy xy;
std::string bmp_filename;

np01_xy_pair_map_init_ab_params init_params = {};
#define NP01_DEFAULT_LOC_GRID_DENSITY (16)
#define NP01_DEFAULT_MAX_LOC_GRID_SQ_COUNT (16384)
init_params.loc_grid_density = NP01_DEFAULT_LOC_GRID_DENSITY;
init_params.max_loc_grid_sq_count = NP01_DEFAULT_MAX_LOC_GRID_SQ_COUNT;

/* parse input file A */
std::ifstream ifs_a(m_xy_a_filename.c_str());
while(ifs_a.good()){
    getline(ifs_a, line);
    x_num_pos = line.find_first_of("+-0123456789.", 0);
    x_end_delim_pos = line.find_first_not_of("+-0123456789.eE", x_num_pos);
    y_num_pos = line.find_first_of("+-0123456789.", x_end_delim_pos);
    y_end_delim_pos = line.find_first_not_of("+-0123456789.eE", y_num_pos);
    if(x_num_pos < line.length()){
        x_num_str = line.substr(x_num_pos, x_end_delim_pos); }
    else{
        x_num_str.erase(); }
    if(y_num_pos < line.length()){
        y_num_str = line.substr(y_num_pos, y_end_delim_pos); }
    else{
        y_num_str.erase(); }
    if(x_num_str.empty() || y_num_str.empty()) {
        /* error */
        }
    else{
        xy.first = atof(x_num_str.c_str());
        xy.second = atof(y_num_str.c_str());
        init_params.xy_vec_a.push_back(xy);
        }

    }
ifs_a.close();
std::cout << "Input A:" << m_xy_a_filename << "  (x,y) count=" 
    << init_params.xy_vec_a.size() << "\n";

/* parse input file B */
std::ifstream ifs_b(m_xy_b_filename.c_str());
while(ifs_b.good()){
    getline(ifs_b, line);
    x_num_pos = line.find_first_of("+-0123456789.", 0);
    x_end_delim_pos = line.find_first_not_of("+-0123456789.eE", x_num_pos);
    y_num_pos = line.find_first_of("+-0123456789.", x_end_delim_pos);
    y_end_delim_pos = line.find_first_not_of("+-0123456789.eE", y_num_pos);
    if(x_num_pos < line.length()){
        x_num_str = line.substr(x_num_pos, x_end_delim_pos); }
    else{
        x_num_str.erase(); }
    if(y_num_pos < line.length()){
        y_num_str = line.substr(y_num_pos, y_end_delim_pos); }
    else{
        y_num_str.erase(); }
    if(x_num_str.empty() || y_num_str.empty()) {
        /* error */
        }
    else{
        xy.first = atof(x_num_str.c_str());
        xy.second = atof(y_num_str.c_str());
        init_params.xy_vec_b.push_back(xy);
        }
    }
ifs_b.close();
std::cout << "Input B:" << m_xy_b_filename << "  (x,y) count=" 
    << init_params.xy_vec_b.size() << "\n";
if( m_max_len_option && ( m_max_len > 0.0 ) ){
    std::cout << "target max length=" << m_max_len << "\n";
    }
else{
    std::cout << "no maximum length\n";
    }

/* Create nodes.  Initialize locator grids for nodes and edges. */
np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);

/* initial assignment */
xy_pair_map.init_edges_far_near();

/* run genetic optimization algorithm */
np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
wksp_init_params.m_target_max_edge_len = (m_max_len_option) ?
    m_max_len : 0.0;
if(m_iterations_option && (m_iterations > 0)){
    wksp_init_params.m_max_iteration_count = m_iterations;
    }
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result = {};
w.get_result(&w_result);
const np01_xy_pair_sub_soln& solution =  w_result.m_solution;

/* write index pair CSV file*/
std::ofstream ofs_ab(m_out_ab_filename.c_str());
ofs_ab << "a,b\n";
const size_t ab_idx_pair_count = solution.get_ab_idx_pair_count();
for(size_t abi = 0; abi < ab_idx_pair_count; ++abi){
    const np01_uint32_pair& ab_idx_pair = solution.get_ab_idx_pair(abi);
    ofs_ab << ab_idx_pair.first << "," << ab_idx_pair.second << "\n";
    }
ofs_ab.close();
std::cout << "Output AB:" << m_out_ab_filename << "  (a,b) count=" 
    << ab_idx_pair_count << "\n";

/* write bitmap file */
std::string out_ab_filename_no_suffix = m_out_ab_filename;
size_t out_ab_filename_suffix_pos = m_out_ab_filename.rfind(".csv");
if((out_ab_filename_suffix_pos != std::string::npos) &&
    (out_ab_filename_suffix_pos > 0)){
    out_ab_filename_no_suffix.resize(out_ab_filename_suffix_pos); }
bmp_filename = out_ab_filename_no_suffix + ".bmp";
xy_pair_map.write_bmp_file(bmp_filename.c_str());
std::cout << "Output bitmap:" << bmp_filename << "\n";

/* output summary */
finish_time = time(NULL);
std::cout << "time elapsed=" << (finish_time - start_time) << " seconds\n";
const np01_xy_pair_sub_soln_cost& cost = solution.get_cost();
const np01_float64 ave_len = (ab_idx_pair_count > 0 ) ?
    (cost.m_edge_len_sum)/static_cast<np01_float64>(ab_idx_pair_count) :
    0.0;
std::cout << "average length = " << ave_len << "\n";
std::cout << "max length = " << (w_result.m_max_result_length) << "\n";
std::cout << "length violation count = " << (w_result.m_len_viol_count) << "\n";
std::cout << "crossover count = " 
    << w_result.m_total_edge_crossover_count << "\n";
}

/* Test 0 */
void np01_main::execute_test_0(){

#define TEST_0_DEFAULT_MAX_ITERATION_COUNT (16)
#define TEST_0_DEFAULT_RAND_SEED (1)
const int max_iteration_count = (m_test_iterations_option) ?
    m_test_iterations : TEST_0_DEFAULT_MAX_ITERATION_COUNT;
const int test_rand_seed = (m_test_rand_seed_option) ?
    m_test_rand_seed : TEST_0_DEFAULT_RAND_SEED;
int error_count = 0;
const time_t start_time = time(NULL);

std::cout << "test 0\n";

np01_uint32 rand_uint32 = test_rand_seed;

for(int iteration=0; iteration < max_iteration_count; ++iteration){
    std::cout << "test 0 iteration=" << iteration << "/" << max_iteration_count << "\n";
    /* choose point-to-point spacing */
    np01_float64 typical_point_spacing;
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 7){
        case 0: typical_point_spacing = 0.001; break;
        case 1: typical_point_spacing = 0.01; break;
        case 2: typical_point_spacing = 0.1; break;
        case 3: typical_point_spacing = 1.0; break;
        case 4: typical_point_spacing = 1.0; break;
        case 5: typical_point_spacing = 1.0; break;
        case 6: default: typical_point_spacing = 10.0; break;
        }

    /* choose size of groups A and B */
    np01_uint32 approx_point_count_a;
    np01_uint32 approx_point_count_b;
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 14){
        case 0: approx_point_count_a = 10; break;
        case 1: approx_point_count_a = 20; break;
        case 2: approx_point_count_a = 50; break;
        case 3: approx_point_count_a = 100; break;
        case 4: approx_point_count_a = 200; break;
        case 5: approx_point_count_a = 500; break;
        case 6: approx_point_count_a = 1000; break;
        case 7: approx_point_count_a = 2000; break;
        case 8: approx_point_count_a = 3000; break;
        case 9: approx_point_count_a = 4000; break;
        case 10: approx_point_count_a = 5000; break;
        case 11: approx_point_count_a = 6000; break;
        case 12: approx_point_count_a = 8000; break;
        case 13: default: approx_point_count_a = 10000; break;
        }
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 5){
        case 0: approx_point_count_b = (approx_point_count_a * 9) / 10; break;
        case 1: approx_point_count_b = approx_point_count_a; break;
        case 2: approx_point_count_b = (approx_point_count_a * 11) / 10; break;
        case 4:
            approx_point_count_b = approx_point_count_a;
            approx_point_count_a = (approx_point_count_b * 9) / 10;
            break;
        case 5: default:
            approx_point_count_b = approx_point_count_a;
            approx_point_count_a = (approx_point_count_b * 11) / 10;
            break;
        }

    /* choose points */
    np01_xy_pair_map_init_ab_params init_params = {};
    init_params.loc_grid_density = NP01_DEFAULT_LOC_GRID_DENSITY;
    init_params.max_loc_grid_sq_count = NP01_DEFAULT_MAX_LOC_GRID_SQ_COUNT;

    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 3){
        case 0: 
            test_init_xy_vec_method_0(approx_point_count_a,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_a));
            break;
        case 1: 
            test_init_xy_vec_method_1(approx_point_count_a,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_a));
            break;
        default: 
            test_init_xy_vec_method_2(approx_point_count_a,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_a));
            break;
        }

    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 3){
        case 0: 
            test_init_xy_vec_method_0(approx_point_count_b,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_b));
            break;
        case 1: 
            test_init_xy_vec_method_1(approx_point_count_b,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_b));
            break;
        default: 
            test_init_xy_vec_method_2(approx_point_count_b,
                typical_point_spacing, (rand_uint32>>16) | (rand_uint32<<16),
                &(init_params.xy_vec_b));
            break;
        }
    const np01_uint32 xy_count_a =
        static_cast<np01_uint32>(init_params.xy_vec_a.size());
    const np01_uint32 xy_count_b =
        static_cast<np01_uint32>(init_params.xy_vec_b.size());
    const np01_uint32 expected_edge_count = static_cast<np01_uint32>(
        (xy_count_a < xy_count_b) ? xy_count_a : xy_count_b);

    /* Create nodes.  Initialize locator grids for nodes and edges. */
    np01::np01_xy_pair_map xy_pair_map; 
    xy_pair_map.init_ab(&init_params);

    /* initial assignment */
    xy_pair_map.init_edges_far_near();

    /* run genetic optimization algorithm */
    np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
    wksp_init_params.set_default(&xy_pair_map);
    wksp_init_params.m_target_max_edge_len = 0.0;
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    switch((rand_uint32 >> 16) % 3){
        case 0: wksp_init_params.m_target_max_edge_len = 0.0;
            break;
        case 1:
            wksp_init_params.m_target_max_edge_len =
                typical_point_spacing * 5.0; 
            break;
        case 2:
        default:
            wksp_init_params.m_target_max_edge_len =
                typical_point_spacing * 10.0; 
            break;
        }
    if((approx_point_count_a > 5000) || (approx_point_count_b > 5000)){
        wksp_init_params.m_max_iteration_count = 32;
        }
    np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
    w.run();
    np01::np01_xy_pair_genetic_wksp_result w_result = {};
    w.get_result(&w_result);
    const np01_xy_pair_sub_soln& solution =  w_result.m_solution;

    /* write bitmap file */
    if( (iteration < 10) ||
        ((iteration < 100) && ((iteration%10) == 0)) ||
        ((iteration < 1000) && ((iteration%100) == 0)) ||
        ((iteration < 10000) && ((iteration%1000) == 0)) ){
        char bmp_file_name[128];
        sprintf(bmp_file_name, "test0_%i.bmp", iteration);
        xy_pair_map.write_bmp_file(bmp_file_name);
        std::cout << "Output bitmap:" << bmp_file_name << "\n";
        }

    /* check */
    const np01_xy_pair_sub_soln_cost& cost = solution.get_cost();
    const np01_float64 ave_len = (expected_edge_count > 0 ) ?
        (cost.m_edge_len_sum)/static_cast<np01_float64>(expected_edge_count) :
        0.0;
    if(xy_pair_map.get_node_map_a()->get_node_count() != xy_count_a){
        ++error_count;
        std::cout << "error: A count mismatch: "
            << (xy_pair_map.get_node_map_a()->get_node_count()) << "/"
            << xy_count_a << "\n";
        }
    if(xy_pair_map.get_node_map_b()->get_node_count() != xy_count_b){
        ++error_count;
        std::cout << "error: B count mismatch: "
            << (xy_pair_map.get_node_map_b()->get_node_count()) << "/"
            << xy_count_b << "\n";
        }
    if(xy_pair_map.get_edge_map()->get_edge_count() != expected_edge_count){
        ++error_count;
        std::cout << "error: edge count mismatch: "
            << (xy_pair_map.get_edge_map()->get_edge_count()) << "/"
            << expected_edge_count << "\n";
        }
    if(wksp_init_params.m_target_max_edge_len > 0.0){
        if((wksp_init_params.m_max_iteration_count > 200) &&
            ( w_result.m_max_result_length > (ave_len * 4.0)) &&
            (w_result.m_len_viol_count >  (expected_edge_count/20)) ) {
            ++error_count;
            std::cout << "error: excessive length violations\n";
            }
        }
    else{
        if((wksp_init_params.m_max_iteration_count > 200) && 
            (w_result.m_total_edge_crossover_count>  (expected_edge_count/10))) {
            ++error_count;
            std::cout << "error: excessive crossovers\n";
            }
        }

    /* output summary */
    std::cout << "size a= " << xy_count_a << "    size b= " << xy_count_b << "\n";
    std::cout << "average length = " << ave_len << "\n";
    std::cout << "max length = " << (w_result.m_max_result_length) << "\n";
    std::cout << "length violation count = " << (w_result.m_len_viol_count) << "\n";
    std::cout << "crossover count = " 
        << w_result.m_total_edge_crossover_count << "\n";
    std::cout << "\n";
    }

const time_t done_time = time(NULL);
std::cout << "test 0 ";
//#if(!defined NDEBUG)
//std::cout << "(debug build) ";
//#endif
std::cout << "done  " << ctime(&done_time);
std::cout << "iterations=" << max_iteration_count << "\n";
std::cout << "rand_seed=" << test_rand_seed << "\n";
std::cout << "elapsed_time=" << (done_time-start_time) << " s\n";
std::cout << "error_count=" << error_count << "\n";
}


/* Test 1 */
void np01_main::execute_test_1(){
std::cout << "test 1\n";

np01::np01_xy_pair_map_init_ab_params init_params = {};
init_params.xy_vec_a.push_back( np01::np01_xy(-1.0,-1.0) );
init_params.xy_vec_a.push_back( np01::np01_xy(-1.0,1.0) );
init_params.xy_vec_a.push_back( np01::np01_xy(1.0,-1.0) );
init_params.xy_vec_a.push_back( np01::np01_xy(1.0,1.0) );
init_params.xy_vec_b.push_back( np01::np01_xy(1.0,0.0) );
init_params.xy_vec_b.push_back( np01::np01_xy(0.0,1.0) );
init_params.xy_vec_b.push_back( np01::np01_xy(0.0,-1.0) );
init_params.loc_grid_density = 2;
init_params.max_loc_grid_sq_count = 12;


np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);
xy_pair_map.init_edges_far_near();
np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result = {};
w.get_result(&w_result);

np01::np01_xy_loc_grid_dim d = xy_pair_map.get_node_map_a()->get_loc_grid_dim();
std::cout << "node_count_a = " << (xy_pair_map.get_node_map_a()->get_node_count()) << "\n";
std::cout << "node_count_b = " << (xy_pair_map.get_node_map_b()->get_node_count()) << "\n";
std::cout << "edge_count = " << (xy_pair_map.get_edge_map()->get_edge_count()) << "\n";
std::cout << "loc grid   w=" << d.get_w() << "  h=" << d.get_h() 
    << "  x_min=" << d.get_x_min() << "  y_min=" << d.get_y_min() 
    << "  sq_sz=" << d.get_sq_size() << "\n";
np01::np01_uint32_vec edge_crossover_count_vec;
np01::np01_uint32 total_edge_crossover_count;
xy_pair_map.get_edge_map()->get_edge_crossover_count(
    &edge_crossover_count_vec, &total_edge_crossover_count);
std::cout << "edge crossover count = " << total_edge_crossover_count << "\n";

if( xy_pair_map.get_edge_map()->get_edge_count() < 100){
    std::cout << "crossovers by edge\n";
    for(size_t ci=0; ci < edge_crossover_count_vec.size(); ++ci){
        std::cout << "["<<ci<<"]" << edge_crossover_count_vec.at(ci) << "\n";
        }
    }
xy_pair_map.write_bmp_file("test1.bmp");
CF01_AA_WKSP_OUTPUT(std::cout);
}

void np01_main::execute_test_2(){
std::cout << "test 2\n";

np01::np01_xy_pair_map_init_ab_params init_params = {};
for(int ii = -4; ii <=4; ++ii ){
    for(int jj = -4; jj <=4; ++jj ){
        init_params.xy_vec_a.push_back( np01::np01_xy(ii, jj) );
        }
    }
for(int rr = 1; rr <=5; rr+=2 ){
    for(int th = 0; th < 360; th += 15 ){
        np01::np01_xy xy;
        xy.first = rr * cos(th * 3.14159/180.0);
        xy.second = rr * sin(th * 3.14159/180.0);
        init_params.xy_vec_b.push_back( xy );
        }
    }
init_params.loc_grid_density = 10;
init_params.max_loc_grid_sq_count = 1000;


np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);
xy_pair_map.init_edges_far_near();
np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
wksp_init_params.m_target_max_edge_len = 3.0;
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result = {};
w.get_result(&w_result);

np01::np01_xy_loc_grid_dim d = xy_pair_map.get_node_map_a()->get_loc_grid_dim();
std::cout << "node_count_a = " << (xy_pair_map.get_node_map_a()->get_node_count()) << "\n";
std::cout << "node_count_b = " << (xy_pair_map.get_node_map_b()->get_node_count()) << "\n";
std::cout << "edge_count = " << (xy_pair_map.get_edge_map()->get_edge_count()) << "\n";
std::cout << "loc grid   w=" << d.get_w() << "  h=" << d.get_h() 
    << "  x_min=" << d.get_x_min() << "  y_min=" << d.get_y_min() 
    << "  sq_sz=" << d.get_sq_size() << "\n";
np01::np01_uint32_vec edge_crossover_count_vec;
np01::np01_uint32 total_edge_crossover_count;
xy_pair_map.get_edge_map()->get_edge_crossover_count(
    &edge_crossover_count_vec, &total_edge_crossover_count);
std::cout << "edge crossover count = " << total_edge_crossover_count << "\n";

if( xy_pair_map.get_edge_map()->get_edge_count() < 100){
    std::cout << "crossovers by edge\n";
    for(size_t ci=0; ci < edge_crossover_count_vec.size(); ++ci){
        std::cout << "["<<ci<<"]" << edge_crossover_count_vec.at(ci) << "\n";
        }
    }

xy_pair_map.write_bmp_file("test2.bmp");

CF01_AA_WKSP_OUTPUT(std::cout);
}

void np01_main::execute_test_3(){
std::cout << "test 3\n";


np01::np01_xy_pair_map_init_ab_params init_params = {};
for(int ii = -20; ii <=20; ++ii ){
    for(int jj = -20; jj <=20; ++jj ){
        init_params.xy_vec_a.push_back( np01::np01_xy(ii, jj) );
        }
    }
for(int rr = 1; rr <=25; rr+=2 ){
    for(int th = 0; th < 360; th += 5 ){
        np01::np01_xy xy;
        xy.first = rr * cos(th * 3.14159/180.0);
        xy.second = rr * sin(th * 3.14159/180.0);
        init_params.xy_vec_b.push_back( xy );
        }
    }
init_params.loc_grid_density = 10;
init_params.max_loc_grid_sq_count = 1000;


np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);
xy_pair_map.init_edges_far_near();
np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
wksp_init_params.m_target_max_edge_len = 5.0;
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result = {};
w.get_result(&w_result);

np01::np01_xy_loc_grid_dim d = xy_pair_map.get_node_map_a()->get_loc_grid_dim();
std::cout << "node_count_a = " << (xy_pair_map.get_node_map_a()->get_node_count()) << "\n";
std::cout << "node_count_b = " << (xy_pair_map.get_node_map_b()->get_node_count()) << "\n";
std::cout << "edge_count = " << (xy_pair_map.get_edge_map()->get_edge_count()) << "\n";
std::cout << "loc grid   w=" << d.get_w() << "  h=" << d.get_h() 
    << "  x_min=" << d.get_x_min() << "  y_min=" << d.get_y_min() 
    << "  sq_sz=" << d.get_sq_size() << "\n";
np01::np01_uint32_vec edge_crossover_count_vec;
np01::np01_uint32 total_edge_crossover_count;
xy_pair_map.get_edge_map()->get_edge_crossover_count(
    &edge_crossover_count_vec, &total_edge_crossover_count);
std::cout << "edge crossover count = " << total_edge_crossover_count << "\n";

if( xy_pair_map.get_edge_map()->get_edge_count() < 100){
    std::cout << "crossovers by edge\n";
    for(size_t ci=0; ci < edge_crossover_count_vec.size(); ++ci){
        std::cout << "["<<ci<<"]" << edge_crossover_count_vec.at(ci) << "\n";
        }
    }

xy_pair_map.write_bmp_file("test3.bmp");

CF01_AA_WKSP_OUTPUT(std::cout);
}

void np01_main::execute_test_4(){
std::cout << "test 4\n";

np01::np01_xy_pair_map_init_ab_params init_params = {};
for(int ii = -50; ii <=50; ++ii ){
    for(int jj = -50; jj <=50; ++jj ){
        init_params.xy_vec_a.push_back( np01::np01_xy(ii, jj) );
        }
    }
for(double rr = 1.0; rr <=49.0; rr+=1.0 ){
    for(double th = 0.0; th < 360.0; th += (50.0/rr)  ){
        np01::np01_xy xy;
        xy.first = rr * cos(th * 3.14159/180.0);
        xy.second = rr * sin(th * 3.14159/180.0);
        init_params.xy_vec_b.push_back( xy );
        }
    }
init_params.loc_grid_density = 16;
init_params.max_loc_grid_sq_count = 4000;


np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);
xy_pair_map.init_edges_far_near();
np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
wksp_init_params.m_target_max_edge_len = 50.0;
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result;
w.get_result(&w_result);

np01::np01_xy_loc_grid_dim d = xy_pair_map.get_node_map_a()->get_loc_grid_dim();
std::cout << "node_count_a = " << (xy_pair_map.get_node_map_a()->get_node_count()) << "\n";
std::cout << "node_count_b = " << (xy_pair_map.get_node_map_b()->get_node_count()) << "\n";
std::cout << "edge_count = " << (xy_pair_map.get_edge_map()->get_edge_count()) << "\n";
std::cout << "loc grid   w=" << d.get_w() << "  h=" << d.get_h() 
    << "  x_min=" << d.get_x_min() << "  y_min=" << d.get_y_min() 
    << "  sq_sz=" << d.get_sq_size() << "\n";
np01::np01_uint32_vec edge_crossover_count_vec;
np01::np01_uint32 total_edge_crossover_count;
xy_pair_map.get_edge_map()->get_edge_crossover_count(
    &edge_crossover_count_vec, &total_edge_crossover_count);
std::cout << "edge crossover count = " << total_edge_crossover_count << "\n";

if( xy_pair_map.get_edge_map()->get_edge_count() < 100){
    std::cout << "crossovers by edge\n";
    for(size_t ci=0; ci < edge_crossover_count_vec.size(); ++ci){
        std::cout << "["<<ci<<"]" << edge_crossover_count_vec.at(ci) << "\n";
        }
    }

xy_pair_map.write_bmp_file("test4.bmp");

CF01_AA_WKSP_OUTPUT(std::cout);
}


void np01_main::execute_test_5(){
std::cout << "test 5\n";

np01::np01_xy_pair_map_init_ab_params init_params = {};
for(int ii = -100; ii <=100; ++ii ){
    for(int jj = -100; jj <=100; ++jj ){
        init_params.xy_vec_a.push_back( np01::np01_xy(ii, jj) );
        }
    }
for(double rr = 1.0; rr <=99.0; rr+=1.0 ){
    for(double th = 0.0; th < 360.0; th += (50.0/rr)  ){
        np01::np01_xy xy;
        xy.first = rr * cos(th * 3.14159/180.0);
        xy.second = rr * sin(th * 3.14159/180.0);
        init_params.xy_vec_b.push_back( xy );
        }
    }
init_params.loc_grid_density = 20;
init_params.max_loc_grid_sq_count = 4000;


np01::np01_xy_pair_map xy_pair_map; 
xy_pair_map.init_ab(&init_params);

#if 0
const np01_xy_loc_grid_dim& loc_grid_dim = 
    xy_pair_map.get_edge_map()->get_loc_grid_dim();
const np01_float64 dsq_noise_ampl = 0.125 *
    static_cast<np01_float64>(loc_grid_dim.get_w()) *
    static_cast<np01_float64>(loc_grid_dim.get_h()) *
    loc_grid_dim.get_sq_size() * loc_grid_dim.get_sq_size();
xy_pair_map.init_edges_ctr_out(dsq_noise_ampl);
std::cout << "dsq_noise_ampl=" << dsq_noise_ampl << "\n";
#endif
xy_pair_map.init_edges_far_near();

np01::np01_xy_pair_genetic_wksp_init_params wksp_init_params;
wksp_init_params.set_default(&xy_pair_map);
wksp_init_params.m_target_max_edge_len = 50.0;
np01::np01_xy_pair_genetic_wksp w(wksp_init_params);
w.run();
np01::np01_xy_pair_genetic_wksp_result w_result = {};
w.get_result(&w_result);

np01::np01_xy_loc_grid_dim d = xy_pair_map.get_node_map_a()->get_loc_grid_dim();
std::cout << "node_count_a = " << (xy_pair_map.get_node_map_a()->get_node_count()) << "\n";
std::cout << "node_count_b = " << (xy_pair_map.get_node_map_b()->get_node_count()) << "\n";
std::cout << "edge_count = " << (xy_pair_map.get_edge_map()->get_edge_count()) << "\n";
std::cout << "loc grid   w=" << d.get_w() << "  h=" << d.get_h() 
    << "  x_min=" << d.get_x_min() << "  y_min=" << d.get_y_min() 
    << "  sq_sz=" << d.get_sq_size() << "\n";
np01::np01_uint32_vec edge_crossover_count_vec;
np01::np01_uint32 total_edge_crossover_count;
xy_pair_map.get_edge_map()->get_edge_crossover_count(
    &edge_crossover_count_vec, &total_edge_crossover_count);
std::cout << "edge crossover count = " << total_edge_crossover_count << "\n";

if( xy_pair_map.get_edge_map()->get_edge_count() < 100){
    std::cout << "crossovers by edge\n";
    for(size_t ci=0; ci < edge_crossover_count_vec.size(); ++ci){
        std::cout << "["<<ci<<"]" << edge_crossover_count_vec.at(ci) << "\n";
        }
    }

xy_pair_map.write_bmp_file("test5.bmp");

CF01_AA_WKSP_OUTPUT(std::cout);
}


void np01_main::test_init_xy_vec_method_0(const np01_uint32& approx_point_count,
    const np01_float64& typical_point_spacing, 
    const np01_uint32& rand_seed, np01_xy_vec *xy_vec){

np01_uint32 rand_uint32 = rand_seed;


const np01_uint32 sqrt_approx_pt_count = static_cast<np01_uint32>(
    sqrt(static_cast<np01_float64>(approx_point_count)));
np01_uint32 w_pt_count=0, h_pt_count=0;

rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 4){
    case 0:
        w_pt_count = sqrt_approx_pt_count;
        h_pt_count = sqrt_approx_pt_count;
        break;
    case 1:
        w_pt_count = (sqrt_approx_pt_count * 5) / 4;
        h_pt_count = approx_point_count/w_pt_count;
        break;
    case 2:
        w_pt_count = (sqrt_approx_pt_count * 3) / 2;
        h_pt_count = approx_point_count/w_pt_count;
        break;
    case 3:
    default:
        w_pt_count = sqrt_approx_pt_count * 2;
        h_pt_count = approx_point_count/w_pt_count;
        break;
    }
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
if(((rand_uint32 >> 16) % 2) == 0){
    np01_uint32 t = w_pt_count;
    w_pt_count = h_pt_count;
    h_pt_count = t;
    }
if(w_pt_count < 1){ w_pt_count = 1; }
if(h_pt_count < 1){ h_pt_count = 1; }


np01_float64 point_spacing;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 4){
    case 0: point_spacing = typical_point_spacing * 0.9; break;
    case 1: point_spacing = typical_point_spacing * 0.95; break;
    case 2: point_spacing = typical_point_spacing; break;
    case 3:
    default:
        point_spacing = typical_point_spacing * 1.05;
        break;
    }

/* center offset  */
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 rand_x_ctr = typical_point_spacing *
    (static_cast<np01_float64>((rand_uint32 >> 16) % 21)-10.0);
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 rand_y_ctr = typical_point_spacing *
    (static_cast<np01_float64>((rand_uint32 >> 16) % 21)-10.0);
np01_float64 x_ctr = 0.0;
np01_float64 y_ctr = 0.0;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 16){
    case 0: x_ctr = rand_x_ctr; y_ctr = 0.0; break;
    case 1: x_ctr = 0.0; y_ctr = rand_y_ctr; break;
    case 2: x_ctr = rand_x_ctr; y_ctr = rand_y_ctr; break;
    case 3:  case 4:  case 5:  case 6:  case 7:
    case 8:  case 9:  case 10:  case 11:  case 12:
    case 13:  case 14:  case 15:  default:
        x_ctr = 0.0; y_ctr = 0.0; break;
        break;
    }

const np01_float64 x_min = x_ctr -
    (static_cast<np01_float64>(w_pt_count-1)*point_spacing/2.0);
const np01_float64 y_min = y_ctr -
    (static_cast<np01_float64>(h_pt_count-1)*point_spacing/2.0);
np01_float64 i,j;
for(i = 0.0; i < static_cast<np01_float64>(w_pt_count); ++i){
    const np01_float64 x = x_min + (i * point_spacing);
    for(j = 0.0; j < static_cast<np01_float64>(h_pt_count); ++j){
        const np01_float64 y = y_min + (j * point_spacing);
        xy_vec->push_back(np01_xy(x,y));
        }
    }
}

void np01_main::test_init_xy_vec_method_1(const np01_uint32& approx_point_count,
    const np01_float64& typical_point_spacing, 
    const np01_uint32& rand_seed, np01_xy_vec *xy_vec){
np01_float64 point_spacing;
np01_uint32 rand_uint32 = rand_seed;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 4){
    case 0: point_spacing = typical_point_spacing * 0.9; break;
    case 1: point_spacing = typical_point_spacing * 0.95; break;
    case 2: point_spacing = typical_point_spacing; break;
    case 3:
    default:
        point_spacing = typical_point_spacing * 1.05;
        break;
    }

/* pi*r^2 = point_count * spacing * spacing 
r = sqrt(point_count * spacing * spacing / pi)
*/
const np01_float64 r_max_0 = sqrt(approx_point_count *
    point_spacing * point_spacing / 3.1415926); ;


np01_float64 r_min = 0.0;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 5){
    case 0: r_min = point_spacing * 0.5; break;
    case 1: point_spacing = typical_point_spacing * 2.5; break;
    case 2: point_spacing = typical_point_spacing * 5.5; break;
    case 3:
    default:
        rand_uint32=(rand_uint32*1103515245)+12345;/* advance rand number */
        r_min=(r_max_0/16.0)*static_cast<np01_float64>((rand_uint32>>16)%24);
        point_spacing = typical_point_spacing * 1.05;
        break;
    }

/* pi*(r_max^2-r_min^2) = point_count * spacing * spacing 
r_max = sqrt(r_min^2 + (point_count*spacing*spacing/pi))
*/
const np01_float64 r_max = sqrt((r_min*r_min)+(approx_point_count *
    point_spacing * point_spacing / 3.1415926)); ;

/* center offset  */
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 rand_x_ctr = typical_point_spacing *
    (static_cast<np01_float64>((rand_uint32 >> 16) % 21)-10.0);
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 rand_y_ctr = typical_point_spacing *
    (static_cast<np01_float64>((rand_uint32 >> 16) % 21)-10.0);
np01_float64 x_ctr = 0.0;
np01_float64 y_ctr = 0.0;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 16){
    case 0: x_ctr = rand_x_ctr; y_ctr = 0.0; break;
    case 1: x_ctr = 0.0; y_ctr = rand_y_ctr; break;
    case 2: x_ctr = rand_x_ctr; y_ctr = rand_y_ctr; break;
    case 3:  case 4:  case 5:  case 6:  case 7:
    case 8:  case 9:  case 10:  case 11:  case 12:
    case 13:  case 14:  case 15:  default:
        x_ctr = 0.0; y_ctr = 0.0; break;
        break;
    }

np01_float64 r,th,th_incr,x,y;
th = 0.0;
for(r = r_min; ((r <= r_max) && (xy_vec->size() < approx_point_count )); 
  r += point_spacing){
    th_incr = (r > 0.0) ? (point_spacing/r) : 6.283185307179;
    for(; th < 6.283185307179; th += th_incr){
        x = x_ctr + (r * cos(th));
        y = y_ctr + (r * sin(th));
        xy_vec->push_back(np01_xy(x,y));
        }
    th -= 6.283185307179;
    }
}


void np01_main::test_init_xy_vec_method_2(const np01_uint32& approx_point_count,
    const np01_float64& typical_point_spacing, 
    const np01_uint32& rand_seed, np01_xy_vec *xy_vec){
np01_float64 point_spacing;
np01_uint32 rand_uint32 = rand_seed;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
switch((rand_uint32 >> 16) % 4){
    case 0: point_spacing = typical_point_spacing * 0.9; break;
    case 1: point_spacing = typical_point_spacing * 0.95; break;
    case 2: point_spacing = typical_point_spacing; break;
    case 3:
    default:
        point_spacing = typical_point_spacing * 1.05;
        break;
    }

np01_float64 x,y;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_uint32 site_point_approx_count = ((rand_uint32 >> 16) % 50) + 10;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_uint32 site_point_w = ((rand_uint32 >> 16) % 11) + 2;
const np01_uint32 site_point_h = site_point_approx_count / site_point_w;
const np01_uint32 site_point_count = site_point_w * site_point_h;
const np01_float64 site_width = point_spacing *
    (static_cast<np01_float64>(site_point_w) - 1.0);
const np01_float64 site_height = point_spacing *
    (static_cast<np01_float64>(site_point_h) - 1.0);
const np01_float64 site_x_min = -site_width/2.0;
const np01_float64 site_y_min = -site_height/2.0;
np01_xy_vec site_xy_vec;
while(site_xy_vec.size() < site_point_count){
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    x = site_x_min + (point_spacing * (1.0/8.0) *
        static_cast<np01_float64>((rand_uint32>>16) % (8 * site_point_w)));
    rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
    y = site_y_min + (point_spacing * (1.0/8.0) *
        static_cast<np01_float64>((rand_uint32>>16) % (8 * site_point_h)));
    site_xy_vec.push_back(np01_xy(x,y));
    }

const np01_float64 site_count_f=static_cast<np01_float64>(approx_point_count) /
    static_cast<np01_float64>(site_point_count);
const np01_float64 sqrt_site_count_f = sqrt(site_count_f);
const np01_uint32 sqrt_site_count_ii =
    static_cast<np01_uint32>(sqrt_site_count_f);
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
np01_uint32 site_ww = sqrt_site_count_ii + ((rand_uint32 >> 16) % 2);
if(site_ww < 1){ site_ww = 1; }
np01_uint32 site_hh = static_cast<np01_uint32>(sqrt_site_count_f) / site_ww;
if(site_hh < 1){ site_hh = 1; }
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
if(((rand_uint32 >> 16) % 2)==0){
    np01_uint32 t = site_ww;
    site_ww = site_hh;
    site_hh = t;
    }

rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 site_spacing_x = site_width + 
    static_cast<np01_float64>((rand_uint32 >> 16) % 4) * point_spacing;
rand_uint32 = (rand_uint32 * 1103515245) + 12345; /* advance rand number */
const np01_float64 site_spacing_y = site_height + 
    static_cast<np01_float64>((rand_uint32 >> 16) % 4) * point_spacing;

const np01_float64 site_min_x =
    -static_cast<np01_float64>(site_ww - 1) * site_spacing_x / 2.0;
const np01_float64 site_min_y =
    -static_cast<np01_float64>(site_hh - 1) * site_spacing_y / 2.0;

np01_uint32 site_i, site_j;
for(site_i = 0; site_i < site_ww; ++site_i){
    const np01_float64 site_ctr_x = site_min_x +
        (static_cast<np01_float64>(site_i) * site_spacing_x);
    for(site_j = 0; site_j < site_hh; ++site_j){
        const np01_float64 site_ctr_y = site_min_y +
            (static_cast<np01_float64>(site_j) * site_spacing_y);
        np01_xy_vec_citr site_xy_itr = site_xy_vec.begin();
        for(; site_xy_itr != site_xy_vec.end(); ++site_xy_itr){
            x = site_ctr_x + site_xy_itr->first;
            y = site_ctr_y + site_xy_itr->second;
            xy_vec->push_back(np01_xy(x,y));
            }
        }
    }
}


/* print to buffer
@param buf (output) character buffer
@param buf_capacity (input) size of character buffer
@param buf_pos (input & output)  input: starting position
  output: one character past last non-0x0 character stored
@param fmt C string that contains the text to be written 
*/
void np01_snprintf( char *buf, const size_t buf_capacity, size_t *buf_pos,
               const char *fmt, ... )
{
if((NULL!=buf) &&  (buf_capacity > 0) && (NULL!=buf_pos) &&
    (*buf_pos < buf_capacity) && (NULL != fmt))
    {
    va_list ap;
    size_t n;
    size_t vn;
    va_start(ap, fmt);
    n = buf_capacity - *buf_pos;
    vn = static_cast<size_t>(vsnprintf(buf+*buf_pos, n, fmt, ap));
    *buf_pos = (vn > n) ? buf_capacity : (*buf_pos + vn);
    va_end(ap);
    }
}

} /* namespace np01 */
