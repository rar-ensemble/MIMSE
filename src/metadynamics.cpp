/*
 * File: metadynamics.cpp
 * Project: Quenching
 * File Created: Wednesday, 26th June 2020 4:55:47 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:28:52 pm
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

void energy_landscape_exploration(const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt)
{
    metadyanmics_relaxation(relaxation_option, dynamics_option, MD_MC_option, dt, 1);

    for (int i = 1; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1))
        {
            bubbles_thermal[i].position = bubbles[i].position;
            bubbles_thermal[i].U = bubbles[i].U;
        }
    }

    bias_counter[1] = 1;
}

void metadyanmics_relaxation(const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt, const int &option)
{
    int count = 0;

    int min_counter_metadynamics = 0;

    counter_phyproc = count;
    time_phyproc = 0;

    do
    {
        output_analysis(2, counter_system, time_system, counter_phyproc);

        ++count;

        relaxation(relaxation_option, dynamics_option, MD_MC_option, dt);

        min_counter_metadynamics += min_counter;

        counter_phyproc = count;

        output_analysis(2, counter_system, time_system, counter_phyproc);

        if ((minima_check() == 1) || (bias_counter[0] == 0))
        {
            bubbles[0].state_flag -= 1;
            break;
        }
        else
        {
            bias_counter[1] = option;
        }

        if (bias_counter[1] != -4)
        {
            metadynamics_reset_system(bias_Z_cut_off);
        }

    } while ((count < METAD_STEP_MAX) && (bias_counter[1] >= -3));

    min_counter = min_counter_metadynamics;
}

void metadynamics_reset_system(int &Z_cut_off)
{
    update_Z();

    update_force_U(type_force, type_bias_force);

    bubbles[0].old_position.setZero(d);

    double random_push = 0.0;

    double random_displacement_push_max = 0;
    double random_displacement_push_max_comp = 0;

    normal_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            if (bias_counter[1] > 0)
            {
                if (bubbles[i].Z[0] >= Z_cut_off)
                {
                    bubbles[i].bias_U_position_list.push_back(bubbles[i].position);
                }
                else
                {
                    bubbles[i].bias_U_position_list.push_back(VectorXd::Constant(1, -1));
                }
            }

            if ((bubbles[i].nl_flag == -1) && (bubbles[i].Z[0] >= Z_cut_off))
            {
                // bubbles[i].old_position = VectorXd::Random(d);
                for (int k = 0; k < d; ++k)
                {
                    bubbles[i].old_position[k] = distribution(generator_4);
                }
            }
            else
            {
                bubbles[i].old_position.setZero(d);
                // mass_j += bubbles[i].mass;
            }
            bubbles[0].old_position += abs(bubbles[i].mass) * bubbles[i].old_position;
        }
        else if (bubbles[i].flag == 0)
        {
            if (bias_counter[1] > 0)
            {
                bubbles[i].bias_U_position_list.push_back(bubbles[i].position);
            }
        }
        else if (i == 0)
        {
            if (bias_counter[1] > 0)
            {
                VectorXd bias_U_position_list = VectorXd::Zero(d + 4);
                bias_U_position_list[0] = 1;
                if (bubbles[0].N_non_rattlers[Z_cut_off] > 0)
                {
                    bias_U_position_list[2] = bubbles[0].N_non_rattlers[Z_cut_off];
                    bias_U_position_list[3] = bubbles[0].mass_non_rattlers[Z_cut_off];
                }
                else
                {
                    cout << "System only has rattlers" << '\n'
                         << "Z_cut_off changed \n";

                    bias_Z_cut_off = 0;
                    Z_cut_off = 0;
                    bias_U_position_list[2] = bubbles[0].N_non_rattlers[Z_cut_off];
                    bias_U_position_list[3] = bubbles[0].mass_non_rattlers[Z_cut_off];
                }

                bubbles[0].bias_U_position_list.push_back(bias_U_position_list);

                ++bias_counter[0];
            }
        }
    }

    if (bias_counter[1] > 0)
    {
        int j = bubbles[0].bias_U_position_list.size() - 1;

        bubbles[0].bias_U_position_com_neighbor_list.push_back(j);
        bubbles[0].bias_U_position_neighbor_list.push_back(j);

        bias_counter[1] = 0;
    }

    bubbles[0].old_position = bubbles[0].old_position / abs(bubbles[0].mass_non_rattlers[Z_cut_off]);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            if (type_force_flag == 1)
            {
                if ((bubbles[i].nl_flag == -1) && (bubbles[i].Z[0] >= Z_cut_off))
                {
                    bubbles[i].old_position -= bubbles[0].old_position;

                    /* if (bubbles[i].z_flag != 0)
                { */
                    random_displacement_push_max = max(random_displacement_push_max, bubbles[i].old_position.norm());

                    random_displacement_push_max_comp = max(random_displacement_push_max_comp, bubbles[i].old_position.cwiseAbs().maxCoeff());
                    /* } */
                }
                else
                {
                    bubbles[i].old_position.setZero(d);
                }
            }
            random_push += bubbles[i].old_position.squaredNorm();
        }
    }

    random_push = sqrt(random_push);

    bubbles[0].old_position[0] = -1;

    double push_alpha = F_tol * eff_F_0 / (4.0 * bias_U_0 / bias_U_sigma / sqrt(bubbles[0].N));

    double push_factor = bias_U_sigma * cardano_real_cubic_equation_solver(-1, push_alpha)[1];

    // double push_factor = bias_U_sigma * sqrt(1.0 - sqrt(1.0 - 10.0 * bubbles[0].N * U_tol));

    // double push_factor = (-2.0 * TOL_MIN * sqrt(bubbles[0].N) + sqrt(4 * gsl_pow_2(TOL_MIN) * bubbles[0].N + gsl_pow_2(10 * U_tol)))/ (2.0 * 10 * U_tol);

    // double push_factor = max(1.0 - sqrt(1.0 - sqrt(bubbles[0].N * U_tol * eff_U_0 / bias_U_0)), sqrt(1.0 / (1 + 4.0 / (U_tol))));

    // push_factor = bias_U_sigma * max(sqrt(10 * TOL_MIN), push_factor);
    
    if (push_factor * random_displacement_push_max_comp / random_push < TOL_MIN)
    {
        push_factor = 2 * TOL_MIN * random_push / random_displacement_push_max_comp;
    }

    if (push_factor > bias_U_sigma)
    {
        cout << "bias_U_sigma too small";
        exit(0);
    }

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].old_position = push_factor * bubbles[i].old_position / random_push;
            bubbles[i].position += bubbles[i].old_position;

            bubbles[i].position = r_BC(type_BC, bubbles[i].position, L_box);

            bubbles[i].deltadisplacement = bubbles[i].old_position;

            bubbles[i].displacement += bubbles[i].old_position;
            bubbles[i].ddisplacement += bubbles[i].old_position;

            bubbles[i].bias_skin_distance += bubbles[i].old_position;

            bubbles[i].skin_distance += bubbles[i].old_position;
        }
    }

    if (bubbles[0].nl_flag < 2)
    {
        bubbles[0].nl_flag = 2;
    }
    if ((bias_counter[0] == 1) && (bubbles[0].bl_flag == -1))
    {
        bubbles[0].bl_flag = 3;
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

    bubbles[0].state_flag = 3;
    // update_center_of_mass();

    update_force_U(type_force, type_bias_force);
}

int minima_check()
{
    tolerances_set(type_dynamics_r, 0, 0);

    update_force_U(type_force, type_bias_force);

    if (bubbles[0].bias_U / eff_U_0 / bubbles[0].N < U_tol)
    {
        int force_flag = 1;
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                if ((bubbles[i].force.squaredNorm() > gsl_pow_2(F_tol * eff_F_0)) || ((bubbles[i].force - bubbles[i].bias_force).squaredNorm() > gsl_pow_2(F_tol * eff_F_0)) || (bubbles[i].bias_force.squaredNorm() > gsl_pow_2(F_tol * eff_F_0)))
                {
                    force_flag = 0;
                    break;
                }
            }
        }

        bubbles[0].dU_U = (1.0 - (bubbles[0].U - bubbles[0].bias_U) / bubbles_athermal[0].U);

        if ((force_flag == 1)) // && (abs(bubbles[0].dU_U) > U_tol))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

void bias_list(const int &bl_dest_flag, const double &bl_skin, const double &bl_end)
{
    if (bubbles[0].bl_flag > -1)
    {
        if (bubbles[0].bl_flag != 3)
        {
            bubbles[0].bias_skin_length = 0.0;

            for (int i = 1; i <= N_init; ++i)
            {
                if (bubbles[i].flag == 1)
                {
                    bubbles[0].bias_skin_length += bubbles[i].bias_skin_distance.squaredNorm();
                }
            }

            // bubbles[0].bias_skin_length = sqrt(bubbles[0].bias_skin_length);

            if ((bubbles[0].bias_skin_length >= gsl_pow_2(bl_skin)) || ((bubbles[0].bl_flag > 0) && (bl_skin == 0)))
            {
                bubbles[0].bl_flag = 3;
            }
        }

        if (bubbles[0].bl_flag > bl_dest_flag)
        {
            if (bubbles[0].bl_flag == 3)
            {
                bubbles[0].bias_skin_length = 0.0;

                bubbles[0].bias_U_position_com_neighbor_list.clear();
                bubbles[0].bias_U_position_neighbor_list.clear();

                for (unsigned int k = 0; k < bubbles[0].bias_U_position_list.size(); ++k)
                {
                    if (bubbles[0].bias_U_position_list[k][0] == 1)
                    {
                        bubbles[0].bias_U_position_list[k][1] = 0.0;
                        for (int i = 1; i <= N_init; ++i)
                        {
                            if (bubbles[i].flag == 1)
                            {
                                bubbles[i].bias_skin_distance.setZero(d);

                                if (bubbles[i].bias_U_position_list[k][0] != -1)
                                {
                                    bubbles[0].bias_U_position_list[k][1] += (periodic_BC(bubbles[i].position, bubbles[i].bias_U_position_list[k], L_box)).squaredNorm();
                                }
                            }
                        }
                        // bubbles[0].bias_U_position_list[k][1] = bubbles[0].bias_U_position_list[k][1] / gsl_pow_2(bias_U_sigma) - 1.0;

                        if (bubbles[0].bias_U_position_list[k][1] <= gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off))) + bl_end))
                        {
                            if (bubbles[0].bias_U_position_list[k][1] <= gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off))) + bl_skin)) // (bubbles[0].bias_U_position_list[k][1] <= bl_skin / bias_U_sigma * (2.0 + bl_skin / bias_U_sigma))
                            {
                                if (bubbles[0].bias_U_position_list[k][1] <= gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off)))))
                                {
                                    bubbles[0].bias_U_position_neighbor_list.push_back(k);
                                }
                                bubbles[0].bias_U_position_com_neighbor_list.push_back(k);
                            }
                        }
                        else
                        {
                            bubbles[0].bias_U_position_list[k][0] = 0;
                        }
                    }
                }
                bubbles[0].bl_flag = 0;
            }
            else if (bubbles[0].bl_flag == 2)
            {
                bubbles[0].bias_U_position_neighbor_list.clear();

                for (unsigned int k = 0; k < bubbles[0].bias_U_position_com_neighbor_list.size(); ++k)
                {
                    int j = bubbles[0].bias_U_position_com_neighbor_list[k];

                    bubbles[0].bias_U_position_list[j][1] = 0;
                    for (int i = 1; i <= N_init; ++i)
                    {
                        if (bubbles[i].flag == 1)
                        {
                            if (bubbles[i].bias_U_position_list[j][0] != -1)
                            {
                                bubbles[0].bias_U_position_list[j][1] += (periodic_BC(bubbles[i].position, bubbles[i].bias_U_position_list[j], L_box)).squaredNorm();
                            }
                        }
                    }
                    // bubbles[0].bias_U_position_list[j][1] = bubbles[0].bias_U_position_list[j][1] / gsl_pow_2(bias_U_sigma) - 1.0;

                    if (bl_dest_flag == 0)
                    {
                        if (bubbles[0].bias_U_position_list[j][1] <= gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off)))))
                        {
                            bubbles[0].bias_U_position_neighbor_list.push_back(j);
                        }
                    }
                }
                if (bl_dest_flag == 0)
                {
                    bubbles[0].bl_flag = 0;
                }
                else
                {
                    bubbles[0].bl_flag = 1;
                }
            }

            if ((bubbles[0].bl_flag == 1) && (bubbles[0].bl_flag > bl_dest_flag))
            {
                bubbles[0].bias_U_position_neighbor_list.clear();

                for (unsigned int k = 0; k < bubbles[0].bias_U_position_com_neighbor_list.size(); ++k)
                {
                    int j = bubbles[0].bias_U_position_com_neighbor_list[k];

                    if (bubbles[0].bias_U_position_list[j][1] <= gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off)))))
                    {
                        bubbles[0].bias_U_position_neighbor_list.push_back(j);
                    }
                }
                bubbles[0].bl_flag = 0;
            }
        }
    }
}