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
void make_edge_set(std::vector <edge>& edge_set, bool** adj_mat, float** weight_mat, int size);
bool** bool2D(const int size);
float** float2D(const int size);
void delete_bool2D(bool **p, int size);
void delete_float2D(float **p, int size);
void init_node_arr(std::vector <edge>& edge_set, node* node_arr, euc_c* coordinates, int size);
void set_heap(node* A, node* B, int heap_size);

#endif /* FUNCTIONS_HPP_ */
