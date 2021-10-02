/*
 * functions.hpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include "user_types.hpp"

void init_adj_mat(bool** adj_mat, int size);
void init_weight_mat(float** weight_mat, int size);
void init_weight_mat(float** weight_mat, float** weight_mat_ref, int size);
void init_coordinates(euc_c* coordinates, euc_c* coordinates_ref, int size);
void print_adj_mat(bool** adj_mat, int size);
void print_weight_mat(float** weight_mat, int size);
void populate_adj_and_weight(bool** adj_mat, euc_c* coordinates, float** weight_mat, int size_graph, float density);
void make_edge_set(std::vector <edge>& edge_set, int** adj_mat, float** weight_mat, int size);
bool** bool2D(const int size);
float** float2D(const int size);
void delete_bool2D(bool **p, int size);
void delete_float2D(float **p, int size);
void init_node_arr(std::vector <edge>& edge_set, node* node_arr, euc_c* coordinates, int size);
void set_heap(node* A, node* B, int heap_size);

int** int2D(const int size);
void fib_heap_insert(FibHeap* H, node* x);
void print_root_circle(node* z);
void make_child_of(FibHeap* H, node* y, node* x);
void consolidate(FibHeap* H);
void print_child_circle(node* child);
void print_circle(node* z);
bool numbers_children_match(node* z, int& num_nodes);
fib_props numbers_match(node* z);
bool is_fib_heap_children(node* z);
void nullify_children_parent_node(node* z);
bool is_fib_heap(node* z);
node* fib_heap_extract_min(FibHeap* H);
void cut(FibHeap* H, node* x, node* y);
void cascading_cut(FibHeap* H, node* y);
void fib_heap_decrease_key(FibHeap* H, node* x, float k);
void set_index_map(int size_graph, int* index_map, int s);
void populate_adj_and_weight_hr(int* index_map, int** adj_mat, float** weight_mat, int size_graph, std::vector<edge>& edges);
void prim(FibHeap* H, float** w, node** v_ref);
float weight_mst(int size_heap, node** v_ref);
void print_mst(int size_heap, node** node_arr);
mst_props mst(int n, std::vector<edge>& edges, int s);
#endif /* FUNCTIONS_HPP_ */
