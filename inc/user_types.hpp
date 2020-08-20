/*
 * user_types.hpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

#include <vector>

const int inf = 3e+8;

typedef struct Euclidead_coordinate {
	float x;
	float y;
} euc_c;

typedef struct Edge {
    int start_vertex;
    int end_vertex;
    float weight;
} edge;

typedef struct Node {
    float key;
    Node* pi;
    euc_c coordinates;
    int index;
    int index_a;
    int parent_index;
    bool in_q;
    std::vector <int> adj_nodes;
} node;

#endif /* USER_TYPES_HPP_ */
