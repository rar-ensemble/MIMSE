/*
 * File: MD.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:18:14 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 9:05:17 pm
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

double MD(const int &MD_option, const double &dt, const int &option, double &f_norm, double &v_norm, double &P)
{
    double status = 0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            status += MD(MD_option, i, dt, option, f_norm, v_norm, P);

            for (int k = 0; k < d; ++k)
            {
                if (abs(bubbles[i].force[k] / (eff_F_0 * sqrt_d)) > F_tol)
                {
                    if (abs(bubbles[i].ddeltadisplacement[k]) < TOL_MIN)
                    {
                        bubbles[i].f_flag = 0;
                        break;
                    }
                }
            }
        }
    }

    return status/ bubbles[0].N;
}

double MD(const int &MD_option, const int &i, const double &dt, const int &option, double &f_norm, double &v_norm, double &P)
{
    switch (MD_option)
    {
    case 1:
        return simple_verlet_MD(i, dt, option, f_norm, v_norm, P);
        break;

    // case 2:
    //     return explicit_Euler_LD(i, dt, option, f_norm);
    //     break;

    // case 3:
    //     return RK_2_LD(i, dt, option);
    //     break;

    case 4:
        return explicit_Euler_thermal_LD(i, dt, option);
        break;

    // case 5:
    //     return RK_2_thermal_LD(i, dt, option);
    //     break;

    // case 6:
    //     return simple_verlet_thermal_NH_MD(i, dt, option, f_norm, v_norm, P);
    //     break;

    case 7:
        return simple_verlet_thermal_LD(i, dt, option, f_norm, v_norm, P);
        break;

    default:
        return simple_verlet_MD(i, dt, option, f_norm, v_norm, P);
        break;
    }
}

double simple_verlet_MD(const int &i, const double &dt, const int &option, double &f_norm, double &v_norm, double &P)
{
    if (option == 1)
    {
        bubbles[i].deltadisplacement = (bubbles[i].force / abs(bubbles[i].mass) * (dt / 2.0));

        bubbles[i].velocity += bubbles[i].deltadisplacement;

        bubbles[i].deltadisplacement *= dt;

        bubbles[i].ddeltadisplacement = bubbles[i].deltadisplacement;

        bubbles[i].deltadisplacement = (bubbles[i].velocity * dt);

        bubbles[i].position += bubbles[i].deltadisplacement;

        bubbles[i].position = r_BC(type_BC, bubbles[i].position, L_box);

        //bubbles[i].force.setZero(d);

        bubbles[i].bias_skin_distance += bubbles[i].deltadisplacement;

        bubbles[i].skin_distance += bubbles[i].deltadisplacement;

        bubbles[i].displacement += bubbles[i].deltadisplacement;
        bubbles[i].ddisplacement += bubbles[i].deltadisplacement;

        bubbles[i].ds = bubbles[i].deltadisplacement.norm();

        bubbles[i].contour += bubbles[i].ds;
        bubbles[i].dcontour += bubbles[i].ds;

        if (bubbles[0].nl_flag < 2)
        {
            bubbles[0].nl_flag = 2;
        }
        if ((bubbles[0].bl_flag > -1) && (bubbles[0].bl_flag < 2))
        {
            bubbles[0].bl_flag = 2;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].com_flag = 1;
        bubbles[0].z_flag = 1;

        return 0.5;
    }
    else if (option == 2)
    {
        bubbles[i].velocity += (bubbles[i].force / abs(bubbles[i].mass) * (dt / 2.0));

        if (f_norm != -1)
        {
            f_norm += bubbles[i].force.squaredNorm();

            v_norm += bubbles[i].velocity.squaredNorm();

            P += bubbles[i].force.dot(bubbles[i].velocity);
        }

        bubbles[i].ddeltadisplacement.setConstant(d, 1);

        return 0.5;
    }
    else
    {
        return 0;
    }
}

double explicit_Euler_thermal_LD(const int &i, const double &dt, const int &option)
{
    if (bubbles[i].mass < 0)
    {
        if (option == 1)
        {
            normal_distribution<double> distribution(0.0, 1.0);
            for (int k = 0; k < d; ++k)
            {
                bubbles[i].eta[k] = distribution(generator_3);
            }
            bubbles[i].eta = sqrt(2.0 * b * kbT * dt) * bubbles[i].eta / dt;

            bubbles[i].deltadisplacement = (bubbles[i].force + bubbles[i].eta) / b;

            bubbles[i].velocity = bubbles[i].deltadisplacement;

            bubbles[i].deltadisplacement *= dt;

            bubbles[i].ddeltadisplacement = bubbles[i].deltadisplacement;

            bubbles[i].position += bubbles[i].deltadisplacement;
            bubbles[i].position = r_BC(type_BC, bubbles[i].position, L_box);

            bubbles[i].bias_skin_distance += bubbles[i].deltadisplacement;

            bubbles[i].skin_distance += bubbles[i].deltadisplacement;

            bubbles[i].displacement += bubbles[i].deltadisplacement;
            bubbles[i].ddisplacement += bubbles[i].deltadisplacement;

            bubbles[i].ds = bubbles[i].deltadisplacement.norm();

            bubbles[i].contour += bubbles[i].ds;
            bubbles[i].dcontour += bubbles[i].ds;

            if (bubbles[0].nl_flag < 2)
            {
                bubbles[0].nl_flag = 2;
            }
            if ((bubbles[0].bl_flag > -1) && (bubbles[0].bl_flag < 2))
            {
                bubbles[0].bl_flag = 2;
            }
            bubbles[0].f_flag = 2;
            bubbles[0].t_flag = 2;
            h_flag = 3;
            r_flag = 2;
            bubbles[0].com_flag = 1;
            bubbles[0].z_flag = 1;

            return 1;
        }
        else if (option == 2)
        {
            bubbles[i].ddeltadisplacement.setConstant(d, 1);

            return 0;
        }
    }
    else
    {
        return 0;
    }
}

double simple_verlet_thermal_LD(const int &i, const double &dt, const int &option, double &f_norm, double &v_norm, double &P)
{
    if (option == 1)
    {
        normal_distribution<double> distribution(0.0, 1.0);
        for (int k = 0; k < d; ++k)
        {
            bubbles[i].eta[k] = distribution(generator_3);
        }
        bubbles[i].eta = 2.0 * sqrt(b * kbT * dt / abs(bubbles[i].mass)) * bubbles[i].eta / dt;

        bubbles[i].deltadisplacement = (bubbles[i].force / abs(bubbles[i].mass) - b * bubbles[i].velocity) + bubbles[i].eta;

        bubbles[i].deltadisplacement *= (dt / 2.0);

        bubbles[i].velocity += bubbles[i].deltadisplacement;

        bubbles[i].deltadisplacement *= dt;

        bubbles[i].ddeltadisplacement = bubbles[i].deltadisplacement;

        bubbles[i].deltadisplacement = (bubbles[i].velocity * dt);

        bubbles[i].position += bubbles[i].deltadisplacement;

        bubbles[i].position = r_BC(type_BC, bubbles[i].position, L_box);

        bubbles[i].bias_skin_distance += bubbles[i].deltadisplacement;

        bubbles[i].skin_distance += bubbles[i].deltadisplacement;

        bubbles[i].displacement += bubbles[i].deltadisplacement;
        bubbles[i].ddisplacement += bubbles[i].deltadisplacement;

        if (output_mode == 5)
        {
            bubbles[i].ds = bubbles[i].deltadisplacement.norm();
        }

        bubbles[i].contour += bubbles[i].ds;
        bubbles[i].dcontour += bubbles[i].ds;

        if (bubbles[0].nl_flag < 2)
        {
            bubbles[0].nl_flag = 2;
        }
        if ((bubbles[0].bl_flag > -1) && (bubbles[0].bl_flag < 2))
        {
            bubbles[0].bl_flag = 2;
        }
        bubbles[0].f_flag = 2;
        bubbles[0].t_flag = 2;
        h_flag = 3;
        r_flag = 2;
        bubbles[0].com_flag = 1;
        bubbles[0].z_flag = 1;

        return 0.5;
    }
    else if (option == 2)
    {
        normal_distribution<double> distribution(0.0, 1.0);
        for (int k = 0; k < d; ++k)
        {
            bubbles[i].eta[k] = distribution(generator_3);
        }
        bubbles[i].eta = 2.0 * sqrt(b * kbT * dt / abs(bubbles[i].mass)) * bubbles[i].eta / dt;

        bubbles[i].velocity += ((bubbles[i].force / abs(bubbles[i].mass) - b * bubbles[i].velocity) + bubbles[i].eta) * (dt / 2.0);

        if (f_norm != -1)
        {
            f_norm += bubbles[i].force.squaredNorm();

            v_norm += bubbles[i].velocity.squaredNorm();

            P += bubbles[i].force.dot(bubbles[i].velocity);
        }

        bubbles[i].ddeltadisplacement.setConstant(d, 1);

        return 0.5;
    }
    else
    {
        return 0;
    }
}

void tolerances_set_MD(const int &MD_option, const double &dt)
{
    switch (MD_option)
    {
    case 1:
        F_tol = max(2 * sqrt_d * TOL_MIN * abs(bubbles[0].mass) / ((gsl_pow_2(dt) / 2.0) * eff_F_0), F_tol);
        break;

    case 2:
        F_tol = max(2 * sqrt_d * TOL_MIN * b / (dt * eff_F_0), F_tol);
        break;

    case 3:
        F_tol = max(2 * sqrt_d * TOL_MIN * b / (dt * eff_F_0), F_tol);
        break;

    case 4:
        F_tol = max(2 * sqrt_d * TOL_MIN * b / (dt * eff_F_0), F_tol);
        break;

    case 5:
        F_tol = max(2 * sqrt_d * TOL_MIN * b / (dt * eff_F_0), F_tol);
        break;

    case 6:
        F_tol = max(2 * sqrt_d * TOL_MIN * abs(bubbles[0].mass) / ((gsl_pow_2(dt) / 2.0) * eff_F_0), F_tol);
        break;

    case 7:
        F_tol = max(2 * sqrt_d * TOL_MIN * abs(bubbles[0].mass) / ((gsl_pow_2(dt) / 2.0) * eff_F_0), F_tol);
        break;

    default:
        break;
    }


    // double U_tol_q = F_tol * eff_F_0 * bubbles[0].N * 100 * 10 * TOL_MIN / eff_U_0; // param_epsilon_ij(type_epsilon) * TOL_MIN / (gsl_pow_2(bubbles[0].radius) * F_tol * eff_F_0); // 4.0 * eff_F_0 / param_epsilon_ij(type_epsilon) * sqrt(d * bubbles[0].N) * TOL_MIN / (100 * F_tol);

    // double U_tol_r = ((bias_counter[1] >= -2) && (bias_U_0 > 0))? gsl_pow_2(100 * F_tol * eff_F_0) : 0;

    double F_tol_q = sqrt(U_tol);

    double F_tol_r = sqrt(1.0 - sqrt(U_tol * bubbles[0].N));

    F_tol_r = 4.0 * (1.0 - gsl_pow_2(F_tol_r)) * F_tol_r / sqrt(bubbles[0].N);

    if ((bias_counter[1] >= -1) || (bias_counter[0] > 0))
    {
        F_tol = max(F_tol_0, max(F_tol_q, F_tol_r));
    }
    else
    {
        F_tol = max(F_tol_0, F_tol_q);
    }

    U_tol = U_tol_0; //max(U_tol_0, max(U_tol_q, U_tol_r));
}