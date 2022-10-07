/*
 * File: FIRE.cpp
 * Project: Quenching
 * File Created: Thursday, 28th June 2020 11:33:18 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:06:49 pm
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

void FIRE(const int &dynamics_option, const int &MD_MC_option, const double &dt)
{
    double U = 0.0;
    double U_old = 0.0;
    double U_conv = 0.0;

    int force_flag = 0;

    double f_norm = 0.0;
    double v_norm = 0.0;

    double P = 0.0;

    double tfire = 0.0;
    double dtfire = dt;
    double dtfiremax = 10 * dtfire;
    int Nmin = 5;
    double finc = 1.1;
    double fdec = 0.5;
    double alphastart = 0.1;
    double falpha = 0.99;
    int ncount = 0;
    int count = 0;
    double count_SUCCESSFUL = 0;

    double alpha = alphastart;

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
    time_phyproc = tfire;

    output_analysis(3, counter_system, time_system, counter_phyproc, time_phyproc);

    do
    {
        U_old = U;
        // U = 0.0;

        if (count > 0)
        {
            tfire += dtfire;

            counter_phyproc = count;
            time_phyproc = tfire;

            f_norm = 0.0;
            v_norm = 0.0;

            P = 0.0;

            count_SUCCESSFUL += update_pos_vel(dynamics_option, MD_MC_option, dtfire, f_norm, v_norm, P);
            U = bubbles[0].U;

            output_analysis(3, counter_system, time_system, counter_phyproc, time_phyproc);
        }

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
            // break;
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
                // cout << setprecision(numeric_limits<double>::digits10) << fixed << U_conv << " " << U << " " << U_old << endl;
                bubbles[0].state_flag = 1;
                break;
            }
        }

        if (count > 0)
        {
            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    bubbles[i].velocity = ((1.0 - alpha) * bubbles[i].velocity) + (alpha * (bubbles[i].force / f_norm) * v_norm);
                }
            }

            if (P > 0.0)
            {
                ++ncount;

                if (ncount > Nmin)
                {
                    dtfire = min(dtfire * finc, dtfiremax);
                    alpha = alpha * falpha;

                    // ncount = 0; // extra addition for careful descent
                }
            }
            else
            {
                ncount = 0;
                dtfire = dtfire * fdec;
                for (int i = 1; i <= N_init; ++i)
                {
                    if (bubbles[i].flag == 1)
                    {
                        bubbles[i].velocity.setZero(d);
                    }
                }

                alpha = alphastart;
            }
        }

    } while (++count <= FIRE_STEP_MAX);

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