/*
 * prim.hpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#ifndef PRIM_HPP_
#define RPIM_HPP_

#include "user_types.hpp"

class Prim {
private:
    int heap_size;
    int length;
    node* heap;
    node* node_array;
    euc_c* coordinates;
    float** weight_mat;
    node* min_node_arr;
    std::vector <edge> edge_set;

    int parent(int i);
    int left(int i);
    int right(int i);
    void min_heapify(node A[], int i);
    void build_min_heap();
    node heap_extract_min();

public:
    Prim(int** adj_mat, float** weight_mat, euc_c* coordinates, int size);
    ~Prim();

    void prim_algo();
    bool min_heap_verify();
    void print_heap();
    void print_mst();
};

#endif /* PRIM_HPP_ */
