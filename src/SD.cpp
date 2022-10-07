/*
 * File: SD.cpp
 * Project: Quenching
 * File Created: Wednesday, 7th November20205:39:09 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 11:44:16 pm
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

void SD(const int &dynamics_option, const int &MD_MC_option, const double &dt)
{
    double U = 0.0;
    double U_old = 0.0;
    double U_conv = 0.0;

    int force_flag = 0;

    int count = 0;
    double count_SUCCESSFUL = 0;

    double gamma = dt;

    double tSD = 0.0;
    double dtSD = b * gamma; //dt_minimizer;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].velocity.setZero(d);
        }
    }

    update_force_U(type_force, type_bias_force);
    U = bubbles[0].U;

    counter_phyproc = count;
    time_phyproc = tSD;

    output_analysis(3, counter_system, time_system, counter_phyproc, time_phyproc);

    do
    {
        U_old = U;
        U = 0.0;

        ++count;
        tSD += dtSD;

        counter_phyproc = count;
        time_phyproc = tSD;

        count_SUCCESSFUL += update_pos_vel(dynamics_option, MD_MC_option, dtSD);
        U = bubbles[0].U;

        output_analysis(3, counter_system, time_system, counter_phyproc, time_phyproc);

        force_flag = 0;
        for (int i = 1; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1) && (bubbles[i].f_flag == 0))
            {
                force_flag = 1;
                break;
            }
        }

        if (force_flag == 1)
        {
            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    bubbles[i].f_flag = -1;
                }
            }
            bubbles[0].state_flag = 1.5;
            break;
        }

        U_conv = abs(U / U_old - 1); //abs((U - U_old) / eff_U_0) / bubbles[0].N;

        if (((U_conv < U_tol) && (count > 0)) || ((abs(U) / eff_U_0 / bubbles[0].N < U_tol) && (type_force == 1)))
        {
            force_flag = 1;
            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    if (bubbles[i].force.squaredNorm() > gsl_pow_2(F_tol * eff_F_0))
                    {
                        force_flag = 0;
                        break;
                    }
                }
            }

            if (force_flag == 1)
            {
                bubbles[0].state_flag = 1;
                break;
            }
        }

    } while (count <= SD_STEP_MAX);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].velocity.setZero(d);
        }
    }

    min_counter_SUCCESSFUL = count_SUCCESSFUL;
    min_counter = count;
}