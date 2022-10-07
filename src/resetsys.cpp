/*
 * File: resetsys.cpp
 * Project: Quenching
 * File Created: Friday, 28th June 2020 4:44:42 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 11:13:01 pm
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

void start_system()
{
    for (int i = 0; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1) || (i == 0))
        {
            bubbles[i].s = 0.0;
            bubbles[i].contour = 0.0;
            bubbles[i].displacement.setZero(d);

            bubbles[i].ds = 0.0;
            bubbles[i].dcontour = 0.0;
            bubbles[i].ddisplacement.setZero(d);

            bubbles[i].deltadisplacement.setZero(d);
        }
    }

    bubbles_athermal = bubbles;
    bubbles_thermal = bubbles;

    time_system += random_seed_pos_init_discard * dt_reset_system;
}

void reset_system(const int &reset_system_option, const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt, const int &distribution_option, const int &initialization_option)
{
    if ((reset_system_option == 1) || (reset_system_option == 2))
    {
        bubbles_thermal[0].state_flag = bubbles[0].state_flag;
        bubbles = bubbles_thermal;
    }
    else if ((reset_system_option >= 1) && (bias_counter[1] <= -1))
    {
        for (int i = 0; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1) || (i == 0))
            {
                bubbles[i].position = bubbles_thermal[i].position;

                if (bias_counter[1] <= -2)
                {
                    bubbles[i].bias_U_position_list = bubbles_thermal[i].bias_U_position_list;

                    if (bias_counter[1] <= -3)
                    {   
                        if (bubbles[i].bias_U_position_list.size() > 0)
                        {
                            bubbles[i].bias_U_position_list.pop_back();
                        }
                    }
                }
            }
        }

        if (bias_counter[1] <= -2)
        {
            bias_counter[0] = bubbles[0].bias_U_position_list.size();
        }
        
        bubbles[0].nl_flag = 3;
        if (bias_counter[0] > 0)
        {
            bubbles[0].bl_flag = 3;
        }
        else
        {
            bubbles[0].bl_flag = -1;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].com_flag = 1;
        bubbles[0].z_flag = 1;
    }

    update_force_U(type_force, type_bias_force);

    counter_phyproc = -1;
    time_phyproc = -1;

    for (int i = 0; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1) || (i == 0))
        {
            // bubbles[i].s = 0.0;
            // bubbles[i].s_non_rattlers = 0.0;
            // bubbles[i].contour = 0.0;
            // bubbles[i].displacement.setZero(d);

            bubbles[i].ds = 0.0;
            bubbles[i].ds_non_rattlers = 0.0;
            bubbles[i].dcontour = 0.0;
            bubbles[i].ddisplacement.setZero(d);

            bubbles[i].U_non_bias_min = INFINITY;
            bubbles[i].U_non_bias_max = -INFINITY;
            bubbles[i].F_max = 0.0;

            if (i > 0)
            {
                if (bubbles[i].Z[0] >= d + 1)
                {
                    bubbles[i].z_flag = 1;
                }
                else
                {
                    bubbles[i].z_flag = 0;
                }
            }
        }
    }

    switch (reset_system_option)
    {
    case 1:
        position_reset_system(initialization_option);
        break;

    case 2:
        relaxation_reset_system(relaxation_option, dynamics_option, MD_MC_option, dt);
        break;

    case 3:
        metadynamics_reset_system(bias_Z_cut_off);
        break;

    // case 4:
    //     dispalcement_reset_system(bias_U_sigma, bias_U_0, bias_Z_cut_off);
    }

    if (bubbles[0].state_flag <= 0.5)
    {
        bubbles[0].state_flag = 4;
    }
    else if (bubbles[0].state_flag <= 1.5)
    {
        bubbles[0].state_flag = 3;
    }
    else
    {
        bubbles[0].state_flag = 4;
    }

    // update_center_of_mass();

    update_force_U(type_force, type_bias_force);

    bubbles[0].s_displacement = 0.0;
    double ds = 0.0;
    double ds_nr = 0.0;

    if (reset_system_option != 2)
    {
        update_force_U(type_force, type_bias_force);

        update_Z();
        for (int i = 1; i <= N_init; ++i)
        {
            // bubbles[i].displacement = bubbles_thermal[i].displacement + bubbles[i].ddisplacement;

            bubbles[i].ds = bubbles[i].ddisplacement.norm();
            bubbles[i].s = bubbles_thermal[i].s + bubbles[i].ds;

            bubbles[i].dcontour = bubbles[i].ds;
            bubbles[i].contour = bubbles_thermal[i].contour + bubbles[i].dcontour;

            bubbles[0].s_displacement += bubbles[i].displacement.squaredNorm();

            ds += gsl_pow_2(bubbles[i].ds);

            if (bubbles[i].Z[0] >= d + 1)
            {
                ds_nr += gsl_pow_2(bubbles[i].ds);
            }

            if (bias_counter[0] > 0)
            {
                bubbles[i].U -= bubbles[i].bias_U;
            }
            bubbles[i].dU = bubbles_thermal[i].U - bubbles[i].U;
            bubbles[i].dU_U = (1.0 - bubbles[i].U / bubbles_thermal[i].U);
        }

        bubbles[0].s_displacement = sqrt(bubbles[0].s_displacement);
        ds = sqrt(ds);
        ds_nr = sqrt(ds_nr);

        bubbles[0].ds = ds;
        bubbles[0].s = bubbles_thermal[0].s + ds;

        bubbles[0].dcontour = bubbles[0].ds;
        bubbles[0].contour = bubbles_thermal[0].contour + bubbles[0].dcontour;

        bubbles[0].ds_non_rattlers = ds_nr;
        bubbles[0].s_non_rattlers = bubbles_thermal[0].s_non_rattlers + ds_nr;

        if (bias_counter[0] > 0)
        {
            bubbles[0].U -= bubbles[0].bias_U;
        }
        bubbles[0].dU = bubbles_thermal[0].U - bubbles[0].U;
        bubbles[0].dU_U = (1.0 - bubbles[0].U / bubbles_thermal[0].U);
    }

    if (reset_system_option >= 3)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                bubbles[i].old_position = bubbles_thermal[i].position;
            }
        }

        if ((unit_vector_flag >= 0) && (bias_counter[0] > 0))
        {
            output_displacement_vector(displacement_vector_filename, counter_system, time_system, counter_phyproc, time_phyproc, bias_counter[0] - 1, 0);
            unit_vector_flag = 1;
        }
    }

    bubbles_thermal = bubbles;

    if (reset_system_option >= 3)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                bubbles_thermal[i].position = bubbles_thermal[i].old_position;
            }
        }
    }

    if (bias_counter[0] > 0)
    {
        for (int i = 0; i <= N_init; ++i)
        {
            bubbles[i].U += bubbles[i].bias_U;
        }
    }

    if (bubbles[0].old_position[0] == -1)
    {
        bubbles[0].s_non_rattlers = ds_nr;
        bubbles[0].ds_non_rattlers = ds_nr;
    }
}

void position_reset_system(const int &initialization_option)
{
    position_initialization(initialization_option);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].ddisplacement = periodic_BC(bubbles[i].position, bubbles_thermal[i].position, L_box);
        }
    }

    bubbles[0].com_flag = 0;

    //start_system();
}