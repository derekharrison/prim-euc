/*
 * main.cpp
 *
 *  Created on: Aug 9, 2020
 *      Author: d-w-h
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>

#include "../inc/functions.hpp"
#include "../inc/prim.hpp"

int main(int argc, char* argv[])
{
    //Declarations
    int s = 0; //Start vertex
    int n = 10; //Number of vertices
    int num_edges = 15; //Number of edges
    euc_c* r_vec = new euc_c[n];
    std::vector< edge > edges;

    //Create coordinates
    float max_r_vec = 10.0;
    srand(time(NULL));
    for(int i = 0; i < n; ++i) {
        float num_x = (float) rand() / RAND_MAX;
        float num_y = (float) rand() / RAND_MAX;
        r_vec[i].x = num_x * max_r_vec;
        r_vec[i].y = num_y * max_r_vec;
    }

    //Create edges
    for(int i = 0; i < num_edges; ++i) {
        int start = rand() % n + 0;
        int end = rand() % n + 0;
        float weight = sqrt((r_vec[start].x - r_vec[end].x)*(r_vec[start].x - r_vec[end].x) +
                            (r_vec[start].y - r_vec[end].y)*(r_vec[start].y - r_vec[end].y));

        edge edge_elem;
        edge_elem.start_vertex = start;
        edge_elem.end_vertex = end;
        edge_elem.weight = weight;
        edges.push_back(edge_elem);
    }

    //Compute minimum spanning tree
    mst_props min_span_props = mst(n, edges, s);

    //Print results
    print_mst(n, min_span_props.node_arr);
    std::cout << "size of minimum spanning tree: " << min_span_props.mst_weight << std::endl;
    std::cout << "done" << std::endl;


    return 0;
}
