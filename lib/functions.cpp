/*
 * functions.cpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>

#include "../inc/user_types.hpp"

void init_adj_mat(bool** adj_mat, int size) {
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j) {
            adj_mat[i][j] = false;
        }
}

void init_weight_mat(float** weight_mat, int size) {
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j) {
            weight_mat[i][j] = 0.0;
        }
}

void init_weight_mat(float** weight_mat, float** weight_mat_ref, int size) {
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j) {
            weight_mat[i][j] = weight_mat_ref[i][j];
        }
}

void print_adj_mat(bool** adj_mat, int size) {
    printf("adjancy matrix\n");
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            printf("(%i,%i): %d ", i, j, adj_mat[i][j]);
        }
        printf("\n");
    }
}

void print_weight_mat(float** weight_mat, int size) {
    printf("weight matrix\n");
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j) {
            printf("(%i,%i): %.2f ", i, j, weight_mat[i][j]);
        }
        printf("\n");
    }
}

void populate_adj_and_weight(bool** adj_mat, euc_c* coordinates, float** weight_mat, int size_graph, float density) {
    init_adj_mat(adj_mat, size_graph);
    init_weight_mat(weight_mat, size_graph);

    srand(time(NULL));
    float max_coordinates = 10.0;
    for(int i = 0; i < size_graph; ++i) {
        float rand_num_x = (float) rand() / RAND_MAX;
        float rand_num_y = (float) rand() / RAND_MAX;
        coordinates[i].x = rand_num_x * max_coordinates;
        coordinates[i].y = rand_num_y * max_coordinates;
    }

    for(int i = 0; i < size_graph; ++i)
        for(int j = i; j < size_graph; ++j) {
            float rand_num = (float) rand() / RAND_MAX;
            if(i != j) {
                adj_mat[j][i] = adj_mat[i][j] = rand_num > (1 - density);
                if(adj_mat[i][j]) {
                    weight_mat[j][i] = weight_mat[i][j] = sqrt((coordinates[i].x - coordinates[j].x)*(coordinates[i].x - coordinates[j].x) +
                                                               (coordinates[i].y - coordinates[j].y)*(coordinates[i].y - coordinates[j].y));
                }
            }
        }
}

void make_edge_set(std::vector <edge>& edge_set, bool** adj_mat, float** weight_mat, int size) {
for(int i = 0; i < size; ++i)
    for(int j = 0; j < size; ++j)
        if(adj_mat[i][j]) {
            edge_set.push_back({i, j, weight_mat[i][j]});
        }
}

bool** bool2D(const int size) {
    bool** p = new bool*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new bool[size];

    return p;
}

float** float2D(const int size) {
    float** p = new float*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new float[size];

    return p;
}

void delete_bool2D(bool **p, int size) {
    for(int i = 0; i < size; ++i)
            delete [] p[i];

    delete [] p;
}

void delete_float2D(float **p, int size) {
    for(int i = 0; i < size; ++i)
            delete [] p[i];

    delete [] p;
}

void init_node_arr(std::vector <edge>& edge_set, node* node_arr, euc_c* coordinates, int size) {
    /* populate node array with index information */
    for(unsigned int i = 0; i < edge_set.size(); ++i) {
        int start_vertex = edge_set[i].start_vertex;
        int end_vertex = edge_set[i].end_vertex;
        node_arr[start_vertex].index = start_vertex;
        node_arr[start_vertex].adj_nodes.push_back(end_vertex);
    }

    /* Initialize node array */
    node_arr[0].key = 0;
    node_arr[0].parent_index = 0;
    node_arr[0].pi = NULL;
    node_arr[0].in_q = true;
    node_arr[0].coordinates = coordinates[0];
    for(int i = 1; i < size; ++i) {
        node_arr[i].key = inf;
        node_arr[i].pi = NULL;
        node_arr[i].in_q = true;
        node_arr[i].coordinates = coordinates[i];
    }
}

void set_heap(node* A, node* B, int heap_size) {
    for(int i = 1; i < heap_size + 1; ++i) {
        A[i] = B[i-1];
    }
}
