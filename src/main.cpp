/*
 * main.cpp
 *
 *  Created on: Aug 9, 2020
 *      Author: d-w-h
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h>

#include "../inc/functions.hpp"
#include "../inc/prim.hpp"

int main(int argc, char* argv[])
{
    int size = 1000;
    float density = 0.5;
    bool** adj_mat = bool2D(size);
    float** weight_mat = float2D(size);
    euc_c* coordinates = new euc_c[size];

    /* Populate adjancy and weight matrices and coordinates with random data */
    populate_adj_and_weight(adj_mat, coordinates, weight_mat, size, density);

    /* Create Prim object */
    Prim myg(adj_mat, weight_mat, coordinates, size);

    /* Execute Prim's algorithm */
    myg.prim_algo();

    /* Print minimum spanning tree */
    myg.print_mst();

    /* Free allocated data */
    delete_bool2D(adj_mat, size);
    delete_float2D(weight_mat, size);
    delete [] coordinates;

    return 0;
}
