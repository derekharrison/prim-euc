/*
 * functions.hpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include "user_types.hpp"

void print_mst(int size_heap, node** node_arr);
mst_props mst(int n, std::vector<edge>& edges, int s);
void delete_node_ref(node **p, int size);

#endif /* FUNCTIONS_HPP_ */
