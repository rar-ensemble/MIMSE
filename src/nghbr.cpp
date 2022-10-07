/*
 * File: nghbr.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:16:20 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:30:59 pm
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

void neighbor_list(const int &neighbor_list_option, const int &nl_dest_flag, const double &nl_skin)
{
    vector<double> nl_length(bubbles[0].N);

    if (bubbles[0].nl_flag != 3)
    {
        if (nl_skin != 0)
        {
            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    bubbles[i].skin_length = bubbles[i].skin_distance.norm();
                    nl_length[bubbles[i].real_index[0] - 1] = bubbles[i].skin_length + bubbles[i].skin_radius_length;
                }
            }
            partial_sort(nl_length.begin(), nl_length.begin() + 2, nl_length.end(), greater<double>());
        }

        if (((nl_length[0] + nl_length[1]) >= nl_skin) || ((nl_skin == 0) && (bubbles[0].nl_flag > 0)))
        {
            bubbles[0].nl_flag = 3;
        }
    }

    switch (neighbor_list_option)
    {
    // case 1:
    //     N2_neighbor(nl_skin, nl_dest_flag);
    //     break;

    case 2:
        N_neighbor(nl_skin, nl_dest_flag, r_force_cut_off);
        break;
    }
}

void exact_neighbor_list()
{
    if (bubbles[0].nl_flag == 2)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                bubbles[i].neighbor_list.clear();
                bubbles[i].non_red_neighbor_list.clear();
            }
        }

        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                for (unsigned int k = 0; k < bubbles[i].non_red_com_neighbor_list.size(); ++k)
                {
                    int j = bubbles[i].non_red_com_neighbor_list[k];

                    if (bubbles[j].flag == 1)
                    {
                        neighbor_criteria(i, j, -1, r_force_cut_off);
                    }
                }
            }
        }

        bubbles[0].nl_flag = 0;
    }
    else if (bubbles[0].nl_flag == 1)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                for (unsigned int j = 0; j < bubbles[i].non_red_neighbor_list.size(); ++j)
                {
                    if (bubbles[bubbles[i].non_red_neighbor_list[j]].flag != 1)
                    {
                        bubbles[i].non_red_neighbor_list.erase(bubbles[i].non_red_neighbor_list.begin() + j);
                    }
                }

                for (unsigned int j = 0; j < bubbles[i].neighbor_list.size(); ++j)
                {
                    if (bubbles[bubbles[i].neighbor_list[j]].flag != 1)
                    {
                        bubbles[i].neighbor_list.erase(bubbles[i].neighbor_list.begin() + j);
                    }
                }
            }
        }

        bubbles[0].nl_flag = 0;
    }

    /* for (int i = 1; i <= N_init; ++i)
    {   if (bubbles[i].flag == 1)
    {
        //sort(bubbles[i].non_red_com_neighbor_list.begin(), bubbles[i].non_red_com_neighbor_list.end());
        sort(bubbles[i].non_red_neighbor_list.begin(), bubbles[i].non_red_neighbor_list.end());
        sort(bubbles[i].neighbor_list.begin(), bubbles[i].neighbor_list.end());}
    } */
}