/*
 * File: BC.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:16:20 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 7:54:28 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#include <headers.h>

#include <globals.h>
#include <classes.h>
#include <functions.h>

VectorXd delta_r(int BC_option, VectorXd r_i, VectorXd r_j, VectorXd L)
{
    switch (BC_option)
    {
    case 1:
        break;

    case 2:
        return periodic_BC(r_i, r_j, L);
        break;

    default:
        return periodic_BC(r_i, r_j, L);
        break;
    }
}

VectorXd r_BC(int BC_option, VectorXd r, VectorXd L)
{
    switch (BC_option)
    {
    case 1:
        break;

    case 2:
        return periodic_BC(r, L);
        break;

    default:
        return periodic_BC(r, L);
        break;
    }
}

VectorXd periodic_BC(const VectorXd &pos_i, const VectorXd &pos_j, const VectorXd &L)
{
    VectorXd mod_pos = (pos_i - pos_j).array() - L.array() * round((pos_i - pos_j).array() / L.array());

    return mod_pos;
}

VectorXd periodic_BC(const VectorXd &pos, const VectorXd &L)
{
    VectorXd mod_pos = pos.array() - L.array() * floor((pos).array() / L.array());

    return mod_pos;
}

/* VectorXi periodic_BC(VectorXi index, VectorXi N_grid)
{
    VectorXi mod_index = index.array() - N_grid.array() * floor(index.array() / N_grid.array());

    return mod_index;
} */

/* VectorXd index_box(VectorXd index, VectorXd N_grid)
{
    VectorXd mod_index = 0.5*((index + N_grid) - (index - N_grid).cwiseAbs());

    return mod_index;
} */

void index_correction(VectorXd &index_1, VectorXd &index_2, const VectorXd &N_grid)
{
    for (int k = 0; k < d; ++k)
    {
        if ((index_1[k] - index_2[k] + 1) >= N_grid[k])
        {
            index_1[k] = N_grid[k] - 1;
            index_2[k] = 0;
        }
    }
}

/* VectorXd stationary_BC(VectorXd pos_i, VectorXd pos_j, VectorXd L)
{
     VectorXd mod_pos = (pos_i - pos_j);

     return mod_pos;
} */
