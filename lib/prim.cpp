/*
 * prim.cpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#include <iostream>

#include "../inc/functions.hpp"
#include "../inc/prim.hpp"
#include "../inc/user_types.hpp"

const char* file_name = "mst_data.txt";

Prim::Prim(bool** adj_mat, float** weight_mat, euc_c* coordinates, int size) {
    this->heap_size = size;
    this->length = size;
    this->heap = new node[size+1];
    this->heap[0].key = inf;
    this->node_array = new node[size];
    this->min_node_arr = new node[size+1];
    this->weight_mat = float2D(size);
    init_weight_mat(this->weight_mat, weight_mat, size);
    make_edge_set(this->edge_set, adj_mat, weight_mat, size);
    init_node_arr(this->edge_set, this->node_array, coordinates, size);
    set_heap(this->heap, this->node_array, size);
    Prim::build_min_heap();
}

Prim::~Prim() {
    delete [] this->heap;
    delete [] this->node_array;
    delete [] this->min_node_arr;
    delete_float2D(this->weight_mat, this->length);
}

int Prim::parent(int i) {
    return i/2;
}

int Prim::left(int i) {
    return 2*i;
}

int Prim::right(int i) {
    return 2*i + 1;
}

void Prim::min_heapify(node A[], int i) {
    int l, r, smallest;
    l = Prim::left(i);
    r = Prim::right(i);
    if(l < this->heap_size + 1 && A[l].key < A[i].key) {
        smallest = l;
    }
    else {
        smallest = i;
    }
    if(r < this->heap_size + 1 && A[r].key < A[smallest].key) {
        smallest = r;
    }
    if(smallest != i) {
        node dummy;
        dummy = A[i];
        A[i] = A[smallest];
        A[smallest] = dummy;
        Prim::min_heapify(A, smallest);
    }
}

void Prim::build_min_heap() {
    for(int i = this->heap_size/2; i > 0; --i) {
        Prim::min_heapify(this->heap, i);
    }
}

bool Prim::min_heap_verify() {
    bool is_min_heap = true;
    for(int i = (this->heap_size - 1)/2; i > 0; --i) {
        int l, r;
        l = Prim::left(i);
        r = Prim::right(i);
        if(this->heap[i].key > this->heap[l].key || this->heap[i].key > this->heap[r].key) {
            is_min_heap = false;
        }
    }

    return is_min_heap;
}

void Prim::print_heap() {
    for(int i = this->heap_size/2; i > 0; --i) {
        int l, r;
        l = Prim::left(i);
        r = Prim::right(i);
        if(l < this->heap_size + 1 && r < this->heap_size + 1) {
            printf("node: %i, key: %.2f, key left child: %.2f, key right child: %.2f\n", i, this->heap[i].key,  this->heap[l].key,  this->heap[r].key);
        }
    }
}

void Prim::print_mst() {
    double total_weight_mst = 0.0;
    for(int i = 1; i < this->length; ++i) {
        int parent_index = this->node_array[i].parent_index;
        int current_index = this->node_array[i].index;
        float weight = this->weight_mat[parent_index][current_index];
        total_weight_mst += weight;
        printf("went from %i to %i totaling: %.4f\n", parent_index, current_index, weight);
    }

    printf("total weight mst: %.8f\n", total_weight_mst);
}

node Prim::heap_extract_min() {
    if(this->heap_size < 1) {
        printf("heap size is less than 1\n");
    }
    node min = this->heap[1];
    this->node_array[this->heap[1].index].in_q = false;
    this->heap[1] = this->heap[this->heap_size];
    this->heap_size = this->heap_size - 1;
    Prim::min_heapify(this->heap, 1);

    return min;
}

void Prim::prim_algo() {
    int* index_map = new int[this->length+1];
    index_map[0] = 0;
    FILE *file_ptr = fopen(file_name, "w");
    /* Main iterations of Prim algorithm */
    for(int it = 0; it < this->length; ++it) {
        node min_node = Prim::heap_extract_min();
        this->min_node_arr[it] = min_node;
        for(int i = 1; i < this->heap_size + 1; ++i) {
            index_map[this->heap[i].index] = i;
        }

        /* Export data to file */
        int parent_index = node_array[min_node.index].parent_index;
        int index = min_node.index;
        fprintf(file_ptr, "%i %f %f %i %f %f %f\n", parent_index,
                                                    node_array[parent_index].coordinates.x,
                                                    node_array[parent_index].coordinates.y,
                                                    node_array[index].index,
                                                    node_array[index].coordinates.x,
                                                    node_array[index].coordinates.y,
                                                    weight_mat[parent_index][index]);

        /* Update index info and keys of neighboring nodes */
        for(unsigned int i = 0; i < min_node.adj_nodes.size(); ++i) {
            int start_vertex = min_node.index;
            int end_vertex = min_node.adj_nodes[i];
            int index_a = index_map[end_vertex];
            node v = node_array[end_vertex];
            if(v.in_q && this->weight_mat[start_vertex][end_vertex] < v.key) {
                node_array[end_vertex].parent_index = min_node.index;
                node_array[end_vertex].index = end_vertex;
                node_array[end_vertex].pi = &min_node_arr[it];
                node_array[end_vertex].key = weight_mat[start_vertex][end_vertex];
                heap[index_a].key = weight_mat[start_vertex][end_vertex];
            }
        }

        Prim::build_min_heap();
    }

    fclose(file_ptr);
}
