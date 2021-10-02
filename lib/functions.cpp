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

int** int2D(const int size) {
    int** p = new int*[size];

    for(int i = 0; i < size; ++i)
        p[i] = new int[size];

    return p;
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

void init_coordinates(euc_c* coordinates, euc_c* coordinates_ref, int size) {
    for(int i = 0; i < size; ++i) {
        coordinates[i] = coordinates_ref[i];
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

void make_edge_set(std::vector <edge>& edge_set, int** adj_mat, float** weight_mat, int size) {
for(int i = 0; i < size; ++i)
    for(int j = 0; j < size; ++j)
        if(adj_mat[i][j] == SETVAR) {
            edge_set.push_back({i, j, weight_mat[i][j]});
        }
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

void fib_heap_insert(FibHeap* H, node* x) {
    x->degree = 0;
    x->p = NULL;
    x->pi = NULL;
    x->child = NULL;
    x->mark = false;

    if(H->min == NULL) {
        x->left = x;
        x->right = x;
        H->min = x;
        H->n = 0;
    }
    else {
        x->left = H->min;
        x->right = H->min->right;
        H->min->right->left = x;
        H->min->right = x;
        if(x->key < H->min->key) {
            H->min = x;
        }
    }

    H->n = H->n + 1;
}

void print_root_circle(node* z) {
    node* xt = z;
    if(xt != NULL) {
        if(xt->right != z) {
            while(xt->right != z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                xt = xt->right;
            }
            if(xt->right == z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
            }
        }
        else {
            std::cout << "X == X->RIGHT" << std::endl;
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
        }
    }
}

void make_child_of(FibHeap* H, node* y, node* x) {

    //Remove node from root list
    y->left->right = y->right;
    y->right->left = y->left;

    if(x->child == NULL) {
        x->child = y;
        y->p = x;
        y->left = y;
        y->right = y;
    }
    else {
        y->left = x->child;
        y->right = x->child->right;
        y->p = x;
        x->child->right->left = y;
        x->child->right = y;
    }

    //Set mark
    y->mark = false;

    x->degree = x->degree + 1;
}

void consolidate(FibHeap* H) {

    double golden = (1.0 + sqrt(5.0)) / 2.0;
    double f = log(H->n) / log(golden);
    int D = floor(f + 0.01) + 1;

    node** A = new node*[D + 2];
    for(int i = 0; i < D + 2; ++i) {
        A[i] = NULL;
    }

    node* x = H->min;
    if(x != NULL) {
        if(x->right != H->min) {

            //Ensure all root nodes have unique degrees
            bool there_is_dup = true;
            while(there_is_dup) {
                there_is_dup = false;
                x = H->min;
                while(x->right != H->min) {
                    int d = x->degree;
                    if(A[d] != NULL && A[d] != x) {
                        there_is_dup = true;
                        node* y = A[d];
                        if(y->key > x->key) {
                            //Make y child of x;
                             make_child_of(H, y, x);

                             A[d] = NULL;
                             A[d+1] = x;

                            if(y == H->min) {
                                H->min = x;
                            }
                        }
                        else {
                            //Make x child of y;
                            make_child_of(H, x, y);

                            A[d] = NULL;
                            A[d+1] = y;

                            if(x == H->min) {
                                H->min = y;
                            }

                            x = y;
                        }
                    }
                    else {
                        A[d] = x;
                    }
                    x = x->right;
                }

                if(x->right == H->min) {
                    int d = x->degree;
                    if(A[d] != NULL && A[d] != x) {
                        there_is_dup = true;
                        node* y = A[d];
                        if(y->key > x->key) {
                            //Make y child of x;
                            make_child_of(H, y, x);

                            A[d] = NULL;
                            A[d+1] = x;

                            if(y == H->min) {
                                H->min = x;
                            }
                        }
                        else {
                            //Make x child of y;
                            make_child_of(H, x, y);

                            A[d] = NULL;
                            A[d+1] = y;

                            if(x == H->min) {
                                H->min = y;
                            }

                            x = y;
                        }
                    }
                    else {
                        A[d] = x;
                    }
                }
            }
        }
        else {
            int d = x->degree;
            A[d] = x;
        }
    }

    //Reconstruct root list
    H->min = NULL;
    for(int i = 0; i < D + 2; ++i) {
        if(A[i] != NULL) {
            if(H->min == NULL) {
                A[i]->left = A[i];
                A[i]->right = A[i];
                A[i]->p = NULL;
                H->min = A[i];
            }
            else {
                A[i]->left = H->min;
                A[i]->right = H->min->right;
                H->min->right->left = A[i];
                H->min->right = A[i];
                A[i]->p = NULL;
                if(A[i]->key < H->min->key) {
                    H->min = A[i];
                }
            }
        }
    }
}

void print_child_circle(node* child) {
    node* xt = child;
    if(xt != NULL) {
        if(xt->right != child) {
            while(xt->right != child) {
                std::cout << "xt->child->key: " << xt->key;
                std::cout << ", xt->child->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
                }
                xt = xt->right;
            }
            if(xt->right == child) {
                std::cout << "xt->child->key: " << xt->key;
                std::cout << ", xt->child->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    if(xt->child != NULL) {
                        std::cout << "xt->child->child->key: " << xt->child->key << std::endl;
                    }
                }
            }
        }
        else {
            std::cout << "X->CHILD == X->CHILD->RIGHT" << std::endl;
            std::cout << "xt->child->key: " << xt->key;
            std::cout << ", xt->child->degree: " << xt->degree << std::endl;
        }
    }
}

void print_circle(node* z) {
    node* xt = z;
    if(xt != NULL) {
        if(xt->right != z) {
            while(xt->right != z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    print_child_circle(xt->child);
                }
                xt = xt->right;
            }
            if(xt->right == z) {
                std::cout << "xt->key: " << xt->key;
                std::cout << ", xt->degree: " << xt->degree << std::endl;
                if(xt->child != NULL) {
                    if(xt->child != NULL) {
                        print_child_circle(xt->child);
                    }
                }
            }
        }
        else {
            std::cout << "X == X->RIGHT" << std::endl;
            std::cout << "xt->key: " << xt->key;
            std::cout << ", xt->degree: " << xt->degree << std::endl;
            if(xt->child != NULL) {
                print_child_circle(xt->child);
            }
        }
    }
}

bool numbers_children_match(node* z, int& num_nodes) {
    bool nums_match = true;
    int num_of_nodes = 0;

    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
            xt = xt->right;
        }
        if(xt->right == z->child) {
            num_of_nodes++;
            if(xt->child != NULL) {
                nums_match = numbers_children_match(xt, num_nodes);
                if(!nums_match) { return false; }
            }
        }

        num_nodes = num_nodes + num_of_nodes;

        if(num_of_nodes == z->degree) { nums_match = true; }
        else { nums_match = false; }
    }

    return nums_match;
}

fib_props numbers_match(node* z) {

    bool nums_match = true;
    int num_nodes = 0;
    fib_props fib_heap_props = { nums_match, num_nodes };

    node* xt = z;
    if(xt != NULL) {
        while(xt->right != z) {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
            xt = xt->right;
        }
        if(xt->right == z) {
            num_nodes++;
            nums_match = numbers_children_match(xt, num_nodes);
            fib_heap_props.deg_is_num_child = nums_match;
            fib_heap_props.num_nodes = num_nodes;
            if(!nums_match) { return fib_heap_props; }
        }
    }

    fib_heap_props.deg_is_num_child = nums_match;
    fib_heap_props.num_nodes = num_nodes;

    return fib_heap_props;
}

bool is_fib_heap_children(node* z) {
    bool is_fibheap = true;

    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
            xt = xt->right;
        }
        if(xt->right == z->child) {
            if(xt->p->key > xt->key) {
                return is_fibheap = false;
            }
            if(xt->child != NULL) {
                is_fibheap = is_fib_heap_children(xt);
                if(!is_fibheap) { return false; }
            }
        }
    }

    return is_fibheap;
}

void nullify_children_parent_node(node* z) {
    node* xt = z->child;
    if(xt != NULL) {
        while(xt->right != z->child) {
            xt->p = NULL;
            xt = xt->right;
        }
        if(xt->right == z->child) {
            xt->p = NULL;
        }
    }
}

bool is_fib_heap(node* z) {
    bool is_fibheap = true;

    node* xt = z;
    if(xt != NULL) {
        while(xt->right != z) {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
            xt = xt->right;
        }
        if(xt->right == z) {
            is_fibheap = is_fib_heap_children(xt);
            if(!is_fibheap) { return false; }
        }
    }

    return is_fibheap;
}

node* fib_heap_extract_min(FibHeap* H) {

    node* z = H->min;

    if(z != NULL) {

        //Add each child of z to root list
        node* y = z->child;
        if(y != NULL) {
            //Set children's parent node to NULL
            nullify_children_parent_node(z);

            y->left->right = z->right;
            z->right->left = y->left;
            y->left = z;
            z->right = y;
            z->degree = 0;

            z->child = NULL;
        }

        //Remove z from root list
        z->left->right = z->right;
        z->right->left = z->left;

        if(z == z->right) {
            H->min = NULL;
        }
        else {

            H->min = z->right;
            consolidate(H);
        }

        H->n = H->n - 1;

    }

    return z;

}

void cut(FibHeap* H, node* x, node* y) {

    //If x is only child set child of parent to null
    if(x == x->right) {
        y->child = NULL;
        y->degree = 0;
    }
    else {
        y->child = x->right;
        y->degree = y->degree - 1;
    }

    //Remove x from child list of y and add x to root list of H
    x->left->right = x->right;
    x->right->left = x->left;

    x->right = H->min->right;
    x->left = H->min;

    H->min->right->left = x;
    H->min->right = x;

    x->p = NULL;
    x->mark = false;
}

void cascading_cut(FibHeap* H, node* y) {
    node* z = y->p;
    if(z != NULL) {
        if(y->mark == false) {
            y->mark = true;
        }
        else {
            cut(H, y, z);
            cascading_cut(H, z);
        }
    }
}

void fib_heap_decrease_key(FibHeap* H, node* x, float k) {
    if(k > x->key) {
        const char* s = "new key is greater than current key";
        std::cout << s << std::endl;
        throw s;
    }

    x->key = k;
    node* y = x->p;
    if(y != NULL && x->key < y->key) {
        cut(H, x, y);
        cascading_cut(H, y);
    }

    if(x->key < H->min->key) {
        H->min = x;
    }
}

void set_index_map(int size_graph, int* index_map, int s) {

    int index_track = 0;
    for(int i = s; i < size_graph; ++i) {
        index_map[i] = index_track;
        index_track++;
    }
    for(int i = 0; i < s; ++i) {
        index_map[i] = index_track;
        index_track++;
    }
}

void create_edges_and_coordinates(bool** adj_mat,
		                          euc_c* coordinates,
		                          int* index_map,
		                          float** weight_mat,
		                          int size_graph,
		                          float density,
		                          std::vector<edge>& edges) {

    srand(time(NULL));
    float max_coordinates = 10.0;
    for(int i = 0; i < size_graph; ++i) {
        float rand_num_x = (float) rand() / RAND_MAX;
        float rand_num_y = (float) rand() / RAND_MAX;
        coordinates[i].x = rand_num_x * max_coordinates;
        coordinates[i].y = rand_num_y * max_coordinates;
    }

    int** elem_is_set = int2D(size_graph);

    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start = index_map[edges[i].start_vertex - 1];
        int end = index_map[edges[i].end_vertex - 1];
        float weight = edges[i].weight;
        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
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

void populate_adj_and_weight_hr(int* index_map, int** adj_mat, float** weight_mat, int size_graph, std::vector<edge>& edges) {

    int** elem_is_set = int2D(size_graph);

    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start = index_map[edges[i].start_vertex - 0];
        int end = index_map[edges[i].end_vertex - 0];
        float weight = edges[i].weight;
        if(elem_is_set[start][end] != SETVAR) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
            elem_is_set[start][end] = elem_is_set[end][start] = SETVAR;
        }
        else if(elem_is_set[start][end] == SETVAR && weight_mat[start][end] >= weight) {
            weight_mat[start][end] = weight_mat[end][start] = weight;
        }
        adj_mat[start][end] = adj_mat[end][start] = SETVAR;
    }
}

bool check_fib_heap(FibHeap* H) {
    /*This is the general test for the fibonacci heap.
      The function returns true if the heap satisfies
      the fibonacci heap properties
     */

    //Compute heap properties
    fib_props fh_props = numbers_match(H->min);
    bool heap_is_fibheap = is_fib_heap(H->min);

    //Check if number of children equal node degrees
    bool deg_is_num_childs = fh_props.deg_is_num_child;

    //Check if number of nodes counted in heap equals H.n
    int num_nodes = fh_props.num_nodes;
    bool num_nodes_match = (num_nodes == H->n);

    //Check to see if heap is properly structured
    bool heap_is_ok = num_nodes_match && deg_is_num_childs && heap_is_fibheap;

    return heap_is_ok;
}


void prim(FibHeap* H, float** w, node** v_ref) {

    //Perform Prim's algorithm
    while(H->n > 0) {
        node* u = fib_heap_extract_min(H);

        //Set u presence in set Q to false
        u->in_q = false;

        int num_adj_nodes = (int) u->adj_nodes.size();
        for(int i = 0; i < num_adj_nodes; ++i) {
            int index_ref = u->adj_nodes[i];
            node* v = v_ref[index_ref];

            if(v->in_q && v->key > w[u->index][v->index]) {
                float weight = w[u->index][v->index];
                v->pi = u;
                v->key = weight;
                fib_heap_decrease_key(H, v, weight);
            }

        }
    }
}

float weight_mst(int size_heap, node** v_ref) {
    float total_weight_mst = 0.0;
    for(int i = 0; i < size_heap; ++i) {
        if(v_ref[i]->pi != NULL) {
            float weight = v_ref[i]->key;
            total_weight_mst += weight;
        }
    }

    return total_weight_mst;
}

void print_mst(int size_heap, node** node_arr) {
    for(int i = 0; i < size_heap; ++i) {
        if(node_arr[i]->pi != NULL) {
            int parent_index = node_arr[i]->pi->index;
            int current_index = node_arr[i]->index;
            float weight = node_arr[i]->key;
            std::cout << "node index: " << current_index << ", ";
            std::cout << "parent index: " << parent_index << ", ";
            std::cout << "weight: " << weight;
            std::cout << std::endl;
        }
    }
}

mst_props mst(int n, std::vector<edge>& edges, int s) {

    //Declarations
    float mst_weight = 0.0;
    FibHeap H;
    const float inf = 3e+8;

    //Set index map
//    s = s - 1;
    int* index_map = new int[n];
    set_index_map(n, index_map, s);

    //Initialize heap
    int num_nodes = n;
    node** v_ref = new node*[num_nodes];
    for(int i = 0; i < num_nodes; ++i) {
        v_ref[i] = new node;
        v_ref[i]->key = inf;
        v_ref[i]->index = i;
        v_ref[i]->in_q = true;
        v_ref[i]->pi = NULL;
        if(i == 0) {
            v_ref[i]->key = 0.0;
        }
        fib_heap_insert(&H, v_ref[i]);
    }

    //Add references to adjacent nodes
    int num_edges = (int) edges.size();
    for(int i = 0; i < num_edges; ++i) {
        int start_index = edges[i].start_vertex;
        int end_index = edges[i].end_vertex;

        int start = index_map[start_index - 0];
        int end = index_map[end_index - 0];
        v_ref[start]->adj_nodes.push_back(end);
        v_ref[end]->adj_nodes.push_back(start);
    }

    //Initialize weight and adjacency matrices
    int** adj_mat = int2D(n);
    float** weight_mat = float2D(n);

    populate_adj_and_weight_hr(index_map, adj_mat, weight_mat, n, edges);

    //Perform Prim's algorithm
    prim(&H, weight_mat, v_ref);

    //Compute MST weight
    mst_weight = weight_mst(n, v_ref);

    //Store MST properties
    mst_props min_span_props;
    min_span_props.mst_weight = mst_weight;
    min_span_props.node_arr = v_ref;

    return min_span_props;
}
