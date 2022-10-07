/*
 * File: init.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:16:20 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:08:13 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#include <globals.h>
#include <functions.h>

void initialization()
{
    input();

    bubbles.clear();

    counter_system = 0;
    time_system = 0.0;

    counter_phyproc = 0;
    time_phyproc = 0.0;

    bubbles.resize(N_init + 1);
    bubbles[0].N = N_init;

    // Hessian_system.resize(d * N_init, d * N_init);
    // Hessianeigenvalues_system.setZero(d * N_init);
    // Hessianeigenvectors_system.setZero(d * N_init, d * N_init);
    h_flag = 3;

    r_ij_system.resize(d * N_init, d * N_init);
    r_flag = 2;

    indexing();

    generator_1.seed(random_seed_size_dist);

    size_distribution(type_distribution);

    mass_initial(type_mass);

    generator_2.seed(random_seed_pos_init);

    for (int j = 1; j <= random_seed_pos_init_discard; ++j)
    {
        position_initialization(1);
    }

    // position_initialization(1);    
    position_initialization(type_initialization);

    generator_3.seed(random_seed_random_number);

    generator_4.seed(random_seed_random_number);

    init_run();

    output();
}

void size_distribution(const int &distribution_option)
{   
    switch (distribution_option)
    {
    case 1:
        uniform_real_dist(para_1, para_2);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 2:
        real_gaussian_dist(para_1, para_2);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 5:
        real_weibull_dist(para_1, para_2);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 7:
        bidisperse_dist(para_1, para_2);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 8:
        monodisperse_dist(para_1);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 10:
        swap_power_law(para_1, para_2);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 11:
        binary_dist(para_1, para_2, para_3);
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].z_flag = 1;
        break;

    case 12:
        file_init();
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].com_flag = 1;
        bubbles[0].z_flag = 1;
        break;
    }

    bubbles[0].radius = mean_radius();
}

void position_initialization(const int &initialization_option)
{   
    switch (initialization_option)
    {
    case 1:
        random_position_init();
        bubbles[0].nl_flag = 3;
        if (bubbles[0].bl_flag > -1)
        {
            bubbles[0].bl_flag = 3;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].com_flag = 1;
        bubbles[0].z_flag = 1;
        break;

    case 2:
        break;
    }
}