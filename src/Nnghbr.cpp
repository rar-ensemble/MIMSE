/*
 * File: Nnghbr.cpp
 * Project: Quenching
 * File Created: Monday, 2nd July 2020 10:17:15 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:32:41 pm
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

void N_neighbor(const double &nl_skin, const int &nl_dest_flag, const double &r_cut_off)
{
    if (bubbles[0].nl_flag > nl_dest_flag)
    {
        if (bubbles[0].nl_flag == 3)
        {
            double L_dim = bubbles[0].radius;//1; //max_radius();
            VectorXd N_grid = (floor(L_box.array() / L_dim));

            vector<vector<vector<vector<int> > > > grid_list;

            if (d == 3)
            {
                grid_list.assign(N_grid[0], vector<vector<vector<int> > >(N_grid[1], vector<vector<int> >(N_grid[2], vector<int>(0))));
            }
            else if (d == 2)
            {
                grid_list.assign(N_grid[0], vector<vector<vector<int> > >(N_grid[1], vector<vector<int> >(1, vector<int>(0))));
            }

            VectorXd L_grid = (L_box.array() / N_grid.array());
            VectorXd grid_index;

            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    bubbles[i].grid_index.setConstant(d, -1);
                    grid_index = floor(bubbles[i].position.array() / L_grid.array());
                    grid_index = r_BC(type_BC, grid_index, N_grid);

                    if (d == 3)
                    {
                        grid_list[grid_index[0]][grid_index[1]][grid_index[2]].push_back(i);
                    }
                    else if (d == 2)
                    {
                        grid_list[grid_index[0]][grid_index[1]][0].push_back(i);
                    }
                    bubbles[i].grid_index = grid_index;

                    bubbles[i].non_red_neighbor_list.clear();
                    bubbles[i].non_red_com_neighbor_list.clear();

                    bubbles[i].com_neighbor_list.clear();
                    bubbles[i].neighbor_list.clear();

                    bubbles[i].skin_distance.setZero(d);
                    bubbles[i].skin_length = 0.0;
                    bubbles[i].skin_radius_length = 0.0;

                    //bubbles[i].force.setZero(d);
                    //bubbles[i].Q_flux = 0.0;
                }
            }

            VectorXd mod_index(d);
            vector<VectorXd> check_index;
            VectorXd check_delta;
            VectorXd check_delta_1;
            VectorXd check_delta_2;

            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    double sigma_ij = param_sigma_ij(type_sigma, i);
                    check_delta = ((r_cut_off * sigma_ij + sigma_tol + nl_skin) * L_grid.cwiseInverse());
                    check_delta = ceil(check_delta.array());
                    check_delta_1 = check_delta + bubbles[i].grid_index;
                    check_delta_2 = -check_delta + bubbles[i].grid_index;
                    index_correction(check_delta_1, check_delta_2, N_grid);

                    for (int l = check_delta_2[0]; l <= check_delta_1[0]; ++l)
                    {
                        for (int m = check_delta_2[1]; m <= check_delta_1[1]; ++m)
                        {
                            if (d == 3)
                            {
                                for (int n = check_delta_2[2]; n <= check_delta_1[2]; ++n)
                                {
                                    mod_index << l, m, n;
                                    mod_index = r_BC(type_BC, mod_index, N_grid);

                                    for (unsigned int j = 0; j < grid_list[mod_index[0]][mod_index[1]][mod_index[2]].size(); ++j)
                                    {
                                        neighbor_criteria(i, grid_list[mod_index[0]][mod_index[1]][mod_index[2]][j], nl_skin, r_cut_off);
                                    }
                                }
                            }
                            else if (d == 2)
                            {
                                mod_index << l, m;
                                mod_index = r_BC(type_BC, mod_index, N_grid);

                                for (unsigned int j = 0; j < grid_list[mod_index[0]][mod_index[1]][0].size(); ++j)
                                {
                                    neighbor_criteria(i, grid_list[mod_index[0]][mod_index[1]][0][j], nl_skin, r_cut_off);
                                }
                            }
                        }
                    }
                    check_index.clear();
                }
            }
            if (nl_skin == 0.0)
            {
                bubbles[0].nl_flag = 0;
            }
            else
            {
                bubbles[0].nl_flag = 2;
            }
        }

        if ((bubbles[0].nl_flag <= 2) && (nl_dest_flag < bubbles[0].nl_flag))
        {
            exact_neighbor_list();
        }
    }
}