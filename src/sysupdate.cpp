/*
 * File: sysupdate.cpp
 * Project: Quenching
 * File Created: Friday, 29th June 2020 11:32:58 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 1:10:53 am
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

double update_pos_vel(const int &dynamics_option, const int &MD_MC_option, const double &dt, double &f_norm, double &v_norm, double &P)
{
    double status = 0;

    bubbles[0].s_displacement = 0.0;
    double dcontour = 0.0;
    double ds = 0.0;
    double ds_unit_vector = 0.0;
    double dcontour_nr = 0.0;
    double ds_nr = 0.0;

    double force_norm = 0.0;

    update_force_U(type_force, type_bias_force);

    for (int i = 0; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].dU = bubbles[i].U;

            force_norm = bubbles[i].force.norm();
            bubbles[i].F_max = max(bubbles[i].F_max, force_norm);
            bubbles[0].F_max = max(bubbles[0].F_max, force_norm);
        }
        else if (i == 0)
        {
            bubbles[i].dU = bubbles[i].U;
            // bubbles[0].F_max = 0.0;
        }
    }

    if (dt != 0)
    {
        status += dynamics(dynamics_option, MD_MC_option, dt, 1, f_norm, v_norm, P);

        // update_center_of_mass();
    }

    update_force_U(type_force, type_bias_force);

    if (dt != 0.0)
    {
        status += dynamics(dynamics_option, MD_MC_option, dt, 2, f_norm, v_norm, P);
    }

    update_force_U(type_force, type_bias_force);

    for (int i = 0; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1) || (i == 0))
        {
            if (i != 0)
            {
                bubbles[i].dU_U = (1.0 - bubbles[i].U / bubbles[i].dU);
                bubbles[i].dU -= bubbles[i].U;

                bubbles[0].s_displacement += bubbles[i].displacement.squaredNorm();

                double ds_i = bubbles[i].ddisplacement.squaredNorm();
                ds += ds_i;
                dcontour += gsl_pow_2(bubbles[i].ds);

                if (bubbles[i].z_flag == 1)
                {
                    ds_nr += ds_i;
                    dcontour_nr += gsl_pow_2(bubbles[i].ds);
                }

                if ((unit_vector_flag >= 1) && (bias_counter[0] > 0))
                {
                    int j = bias_counter[0] - 1;

                    if (bubbles[i].bias_U_position_list[j][0] != -1)
                    {
                        ds_unit_vector += ds_i;
                    }
                }

                force_norm = bubbles[i].force.norm();
                bubbles[i].F_max = max(bubbles[i].F_max, force_norm);
                bubbles[0].F_max = max(bubbles[0].F_max, force_norm);
            }
            if (i == 0)
            {
                // bubbles[0].F_max = 0.0;
            }

            bubbles[i].U_non_bias_min = min(bubbles[i].U_non_bias_min, bubbles[0].U - (bias_counter[0] > 0) * bubbles[0].bias_U);
            bubbles[i].U_non_bias_max = max(bubbles[i].U_non_bias_max, bubbles[0].U - (bias_counter[0] > 0) * bubbles[0].bias_U);
        }
    }
    bubbles[0].dU_U = (1.0 - bubbles[0].U / bubbles[0].dU);
    bubbles[0].dU -= bubbles[0].U;

    bubbles[0].s_displacement = sqrt(bubbles[0].s_displacement);
    dcontour = sqrt(dcontour);
    ds = sqrt(ds);
    dcontour_nr = sqrt(dcontour_nr);
    ds_nr = sqrt(ds_nr);

    bubbles[0].dcontour += dcontour;
    bubbles[0].contour += dcontour;

    bubbles[0].s = ds;
    bubbles[0].ds = dcontour;

    bubbles[0].s_non_rattlers = ds_nr;
    bubbles[0].ds_non_rattlers = dcontour_nr;

    if ((unit_vector_flag >= 1) && (bias_counter[0] > 0))
    {
        ds_unit_vector = sqrt(ds_unit_vector);

        while (ds_unit_vector >= unit_vector_r[unit_vector_flag - 1])
        {
            output_displacement_vector(displacement_vector_filename, counter_system, time_system, counter_phyproc, time_phyproc, bias_counter[0] - 1, unit_vector_r[unit_vector_flag - 1]);

            if (unit_vector_flag < unit_vector_r.size())
            {
                ++unit_vector_flag;
            }
            else
            {
                unit_vector_flag = 0;
                break;
            }
        }
    }

    return status;
}

void update_force_U(const int &force_option, const int &bias_force_option)
{
    if (bubbles[0].f_flag >= 1)
    {
        switch (force_option)
        {
        case 1:
            U_0 = param_epsilon_ij(type_epsilon) / force_para_1; //
            F_0 = U_0 * force_para_1 / param_sigma_ij(type_sigma);
            break;

        case 2:
            U_0 = param_epsilon_ij(type_epsilon) / force_para_1; //
            F_0 = U_0 * force_para_1 / param_sigma_ij(type_sigma);
            break;

        case 3:
            static const double k_2 = force_para_2 / (force_para_1 - force_para_2) * pow(force_para_2 / force_para_1, force_para_1 / (force_para_2 - force_para_1));
            U_0 = k_2 * param_epsilon_ij(type_epsilon); //
            F_0 = U_0 / param_sigma_ij(type_sigma);
            break;

        case 4:
            static const double k_3 = force_para_2 / (force_para_1 - force_para_2) * pow(force_para_2 / force_para_1, force_para_1 / (force_para_2 - force_para_1));
            U_0 = k_3 * param_epsilon_ij(type_epsilon); //
            F_0 = U_0 / param_sigma_ij(type_sigma);
            break;

        default:
            U_0 = param_epsilon_ij(type_epsilon) / force_para_1; //
            F_0 = U_0 * force_para_1 / param_sigma_ij(type_sigma);
            break;
        }

        bias_F_0 = bias_U_0 / bias_U_sigma;

        eff_U_0 = ((bias_counter[0] > 0)) ? min(U_0, bias_U_0 / bubbles[0].N) : U_0;
        eff_F_0 = ((bias_counter[0] > 0)) ? min(F_0, bias_F_0) : F_0;

        if (bubbles[0].f_flag == 2)
        {
            for (int i = 0; i <= N_init; ++i)
            {
                if ((bubbles[i].flag == 1) || (i == 0))
                {
                    bubbles[i].force.setZero(d);
                    bubbles[i].U = 0.0;

                    bubbles[i].bias_force.setZero(d);
                    bubbles[i].bias_U = 0.0;
                }
            }
            bubbles[0].f_flag = 1;

            if (bias_counter[0] > 0)
            {
                for (unsigned int k = 0; k < bubbles[0].bias_U_position_list.size(); ++k)
                {
                    if (bubbles[0].bias_U_position_list[k][0] == 1)
                    {
                        bubbles[0].bias_U_position_list[k].tail(d) = VectorXd::Zero(d);
                    }
                }
            }
        }

        if ((particle_fix_flag != 0) || (type_force_flag == 1))
        {
            update_Z();
        }

        neighbor_list(type_neighbor_list, 2, nl_para);

        if (bias_counter[0] > 0)
        {
            bias_list(1, bl_para[0], bl_para[1]);
        }

        for (int i = 0; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1) || (i == 0))
            {
                update_force_U_i(force_option, i, bubbles[0].nl_flag);
            }

            if (bias_counter[0] > 0)
            {
                if ((bubbles[i].flag == 1) || (i == 0))
                {
                    update_bias_force_U_i(bias_force_option, i, bubbles[0].bl_flag);
                }
            }
        }

        update_force_correction();

        bubbles[0].f_flag = 0;
    }
}

void update_force_U_i(const int &force_option, const int &i, const int &nl_flag)
{
    VectorXd force_U_r_ij(2 * d + 1);

    if (i != 0)
    {
        if ((nl_flag == 0) || (nl_flag == 1))
        {
            for (unsigned int k = 0; k < bubbles[i].non_red_neighbor_list.size(); ++k)
            {
                int j = bubbles[i].non_red_neighbor_list[k];

                if (bubbles[j].flag == 1)
                {
                    force_U_r_ij = force(force_option, i, j, force_para_1, force_para_2, r_force_cut_off);

                    bubbles[0].U += force_U_r_ij[d];
                    bubbles[i].U += force_U_r_ij[d];
                    bubbles[j].U += force_U_r_ij[d];

                    bubbles[i].force += force_U_r_ij.head(d);
                    bubbles[j].force += (-1.0 * force_U_r_ij.head(d));

                    bubbles[0].force += force_U_r_ij.head(d);
                    bubbles[0].force += (-1.0 * force_U_r_ij.head(d));
                }
            }
        }
        else if (nl_flag == 2)
        {
            for (unsigned int k = 0; k < bubbles[i].non_red_com_neighbor_list.size(); ++k)
            {
                int j = bubbles[i].non_red_com_neighbor_list[k];

                if (bubbles[j].flag == 1)
                {
                    force_U_r_ij = force(force_option, i, j, force_para_1, force_para_2, r_force_cut_off);

                    bubbles[0].U += force_U_r_ij[d];
                    bubbles[i].U += force_U_r_ij[d];
                    bubbles[j].U += force_U_r_ij[d];

                    bubbles[i].force += force_U_r_ij.head(d);
                    bubbles[j].force += (-1.0 * force_U_r_ij.head(d));

                    bubbles[0].force += force_U_r_ij.head(d);
                    bubbles[0].force += (-1.0 * force_U_r_ij.head(d));
                }
            }
        }
        else if (nl_flag == -1)
        {
            for (unsigned int k = 0; k < bubbles[i].neighbor_list.size(); ++k)
            {
                int j = bubbles[i].neighbor_list[k];

                if (bubbles[j].flag == 1)
                {
                    force_U_r_ij = force(force_option, i, j, force_para_1, force_para_2, r_force_cut_off);

                    bubbles[i].U += force_U_r_ij[d];

                    bubbles[i].force += force_U_r_ij.head(d);
                }
            }
        }
        else if (nl_flag == -2)
        {
            for (unsigned int k = 0; k < bubbles[i].com_neighbor_list.size(); ++k)
            {
                int j = bubbles[i].com_neighbor_list[k];

                if (bubbles[j].flag == 1)
                {
                    force_U_r_ij = force(force_option, i, j, force_para_1, force_para_2, r_force_cut_off);

                    bubbles[i].U += force_U_r_ij[d];

                    bubbles[i].force += force_U_r_ij.head(d);
                }
            }
        }
    }
    else
    {
        force(force_option, i, 0, force_para_1, force_para_2, r_force_cut_off);
    }
}

void update_bias_force_U_i(const int &bias_force_option, const int &i, const int &bl_flag)
{
    VectorXd force_U_i(d + 1);

    if (bl_flag == 0)
    {
        if (i == 0)
        {
            for (unsigned int k = 0; k < bubbles[0].bias_U_position_neighbor_list.size(); ++k)
            {
                int j = bubbles[0].bias_U_position_neighbor_list[k];

                force_U_i = force(bias_force_option, i, j);

                bubbles[0].U += force_U_i[d];
                bubbles[0].bias_U += force_U_i[d];
            }
        }
        else
        {
            for (unsigned int k = 0; k < bubbles[0].bias_U_position_neighbor_list.size(); ++k)
            {
                int j = bubbles[0].bias_U_position_neighbor_list[k];

                force_U_i = force(bias_force_option, i, j);

                bubbles[i].force += force_U_i.head(d);
                bubbles[i].bias_force += force_U_i.head(d);

                bubbles[0].force += force_U_i.head(d);
                bubbles[0].bias_force += force_U_i.head(d);

                bubbles[0].bias_U_position_list[j].tail(d) += force_U_i.head(d);
            }
        }
    }
    else if (bl_flag == 1)
    {
        if (i == 0)
        {
            for (unsigned int k = 0; k < bubbles[0].bias_U_position_com_neighbor_list.size(); ++k)
            {
                int j = bubbles[0].bias_U_position_com_neighbor_list[k];

                force_U_i = force(bias_force_option, i, j);

                bubbles[0].U += force_U_i[d];
                bubbles[0].bias_U += force_U_i[d];
            }
        }
        else
        {
            for (unsigned int k = 0; k < bubbles[0].bias_U_position_com_neighbor_list.size(); ++k)
            {
                int j = bubbles[0].bias_U_position_com_neighbor_list[k];

                force_U_i = force(bias_force_option, i, j);

                bubbles[i].force += force_U_i.head(d);
                bubbles[i].bias_force += force_U_i.head(d);

                bubbles[0].force += force_U_i.head(d);
                bubbles[0].bias_force += force_U_i.head(d);

                bubbles[0].bias_U_position_list[j].tail(d) += force_U_i.head(d);
            }
        }
    }
}

void update_force_correction()
{
    update_N();

    if ((particle_fix_flag != 0) || (type_force_flag == 1))
    {
        update_Z();
    }

    if (particle_fix_flag != 0)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1) && (bubbles[i].nl_flag == 0))
            {
                if (bubbles[i].Z[0] >= 1)
                {
                    bubbles[0].force -= bubbles[i].force / (1.0 - abs(bubbles[i].mass) / abs(bubbles[0].mass_non_rattlers[1]));
                    bubbles[i].force.setZero(d);
                }

                if (bias_counter[0] > 0)
                {
                    bubbles[0].bias_force -= bubbles[i].bias_force / (1.0 - abs(bubbles[i].mass) / abs(bubbles[0].mass_non_rattlers[d + 1]));
                    bubbles[i].bias_force.setZero(d);
                }
                break;
            }
        }
    }

    if (type_force_flag == 1)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1) && (bubbles[i].nl_flag != 0))
            {
                if (bubbles[i].Z[0] >= 1)
                {
                    bubbles[i].force -= abs(bubbles[i].mass) * (bubbles[0].force - bubbles[0].bias_force) / abs(bubbles[0].mass_non_rattlers[1]);
                }

                if (bias_counter[0] > 0)
                {
                    bias_list(1, bl_para[0], bl_para[1]);

                    if (bubbles[0].bl_flag == 0)
                    {
                        if (bubbles[i].flag == 1)
                        {
                            for (unsigned int k = 0; k < bubbles[0].bias_U_position_neighbor_list.size(); ++k)
                            {
                                int j = bubbles[0].bias_U_position_neighbor_list[k];

                                if (bubbles[i].bias_U_position_list[j][0] != -1)
                                {
                                    bubbles[i].force -= abs(bubbles[i].mass) * bubbles[0].bias_U_position_list[j].tail(d) / abs(bubbles[0].bias_U_position_list[j][3]);
                                    bubbles[i].bias_force -= abs(bubbles[i].mass) * bubbles[0].bias_U_position_list[j].tail(d) / abs(bubbles[0].bias_U_position_list[j][3]);
                                }
                            }
                        }
                    }
                    else if (bubbles[0].bl_flag == 1)
                    {
                        if (bubbles[i].flag == 1)
                        {
                            for (unsigned int k = 0; k < bubbles[0].bias_U_position_com_neighbor_list.size(); ++k)
                            {
                                int j = bubbles[0].bias_U_position_com_neighbor_list[k];

                                if (bubbles[i].bias_U_position_list[j][0] != -1)
                                {
                                    bubbles[i].force -= abs(bubbles[i].mass) * bubbles[0].bias_U_position_list[j].tail(d) / abs(bubbles[0].bias_U_position_list[j][3]);
                                    bubbles[i].bias_force -= abs(bubbles[i].mass) * bubbles[0].bias_U_position_list[j].tail(d) / abs(bubbles[0].bias_U_position_list[j][3]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void update_N()
{
    if (bubbles[0].n_flag == 1)
    {
        bubbles[0].N = 0;

        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                ++bubbles[0].N;
                bubbles[i].real_index[0] = bubbles[0].N;
            }
            else if ((bubbles[i].flag == 0) && (bubbles[i].n_flag == 0))
            {
                bubbles[i].n_flag = -1;

                bubbles[i].radius = 0;

                bubbles[i].mass = 0;

                bubbles[i].velocity.setZero(d);

                bubbles[i].force.setZero(d);

                bubbles[i].U = 0;

                bubbles[i].bias_force.setZero(d);

                bubbles[i].bias_U = 0;

                bubbles[i].ddisplacement.setZero(d);
                bubbles[i].ds = 0.0;
                bubbles[i].dcontour = 0.0;

                bubbles[i].KE = 0.0;
                bubbles[i].E = 0.0;

                bubbles[i].Tau.setZero(d, d);

                bubbles[i].Q_flux = 0.0;

                bubbles[i].skin_distance.setZero(d);
                bubbles[i].skin_length = 0.0;
                bubbles[i].skin_radius_length = 0.0;

                bubbles[i].neighbor_list.clear();
                bubbles[i].non_red_neighbor_list.clear();

                bubbles[i].com_neighbor_list.clear();
                bubbles[i].non_red_com_neighbor_list.clear();

                bubbles[i].Z.assign(1, 0);

                bubbles[i].real_index.assign(1, -1);
                bubbles[i].grid_index.setConstant(d, -1);

                bubbles[i].flag = 0;
            }
        }

        bubbles[0].n_flag = 0;

        bubbles[0].z_flag = 1;
    }
}

void update_kbT(const int &system_counter, const double &system_time)
{
    switch (type_kbT)
    {
    case 1:
        dkbT = 0.0;
        break;

    case 2:
        if (system_counter <= NUMBER_OF_EQB_STEPS)
        {
            kbT = kbTi;
        }
        else
        {
            kbT = kbTf;
        }
        dkbT = 0.0;
        break;

    case 3:
        if (system_counter <= NUMBER_OF_EQB_STEPS)
        {
            dkbT = (kbTf - kbTi) / NUMBER_OF_EQB_STEPS;
            kbT += dkbT;
        }
        else
        {
            dkbT = 0.0;
            kbT = kbTf;
        }
        break;

    case 4:
        if (system_counter <= NUMBER_OF_EQB_STEPS)
        {
            dkbT = 0.0;
            kbT = kbTi;
        }
        else
        {
            dkbT = (kbTf - kbTi) / NUMBER_OF_QUENCHES;
            kbT += dkbT;
        }
        break;

    default:
        dkbT = 0.0;
        break;
    }
}

void update_Z()
{
    if (bubbles[0].z_flag == 1)
    {
        neighbor_list(type_neighbor_list, 0, nl_para);

        bubbles[0].N_non_rattlers.assign(2 * d + 1, 0);
        bubbles[0].mass_non_rattlers.assign(2 * d + 1, 0);
        bubbles[0].Z.assign(2 * d + 1, 0.0);
        vector<int> rattlers(2 * d + 1, 0);

        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                bubbles[i].real_index.assign(2 * d + 1, -1);
                bubbles[i].Z.assign(1, bubbles[i].neighbor_list.size());

                for (int j = 0; j <= 2 * d; ++j)
                {
                    if (bubbles[i].Z[0] >= j)
                    {
                        bubbles[0].Z[j] += bubbles[i].Z[0];

                        ++(bubbles[0].N_non_rattlers[j]);

                        bubbles[0].mass_non_rattlers[j] += bubbles[i].mass;

                        bubbles[i].real_index[j] = bubbles[0].N_non_rattlers[j];
                    }
                }
            }
        }

        for (int j = 0; j < (2 * d + 1); ++j)
        {
            bubbles[0].Z[j] = (bubbles[0].Z[j] / (bubbles[0].N_non_rattlers[j])) - 2.0 * d;
        }

        bubbles[0].z_flag = 0;
    }
}

void mass_conservation()
{
    mass_N = 0.0;
    double mass = 0.0;
    double rho = 1.0;

    double delta_R;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            mass_N += (rho * bubbles[i].volume);
        }
    }

    if (mass_N != mass_init)
    {
        for (int i = 1; i <= N_init; ++i)
        {
            if (bubbles[i].flag == 1)
            {
                delta_R = -1.0 * bubbles[i].radius;
                mass = rho * bubbles[i].volume + ((mass_init - mass_N) / bubbles[0].N);
                if (d == 2)
                {
                    bubbles[i].radius = sqrt((d / (2.0 * (d - 1) * pi)) * mass / rho);
                }
                else if (d == 3)
                {
                    bubbles[i].radius = cbrt((d / (2.0 * (d - 1) * pi)) * mass / rho);
                }
                else
                {
                    bubbles[i].radius = pow((d / (2.0 * (d - 1) * pi)) * mass / rho, 1.0 / d);
                }
                bubbles[i].volume = mass / rho;
                // bubbles[i].mass = 1.0;
                if (mass <= 0)
                {
                    cout << "error in mass cons";
                    exit(0);
                }
                delta_R += bubbles[i].radius;
                bubbles[i].skin_radius_length += delta_R;

                if (delta_R != 0.0)
                {
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
                    bubbles[0].com_flag = 1;
                    bubbles[0].z_flag = 1;
                }
            }
        }
    }

    bubbles[0].radius = mean_radius();
}

void update_center_of_mass()
{
    if (bubbles[0].com_flag == 1)
    {
        update_N();

        for (int i = 0; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1))
            {
                bubbles[0].position += abs(bubbles[i].mass) * bubbles[i].position / (abs(bubbles[0].mass) * bubbles[0].N);
                bubbles[0].displacement += abs(bubbles[i].mass) * bubbles[i].deltadisplacement / (abs(bubbles[0].mass) * bubbles[0].N);
                bubbles[0].deltadisplacement += abs(bubbles[i].mass) * bubbles[i].deltadisplacement / (abs(bubbles[0].mass) * bubbles[0].N);
            }
            else if (i == 0)
            {
                bubbles[0].position.setZero(d);
                bubbles[0].deltadisplacement.setZero(d);
            }
        }

        for (int i = 1; i <= N_init; ++i)
        {
            if ((bubbles[i].flag == 1))
            {
                bubbles[i].position += -1.0 * bubbles[0].deltadisplacement;
                bubbles[i].position = r_BC(type_BC, bubbles[i].position, L_box);

                bubbles[i].bias_skin_distance += -1.0 * bubbles[0].deltadisplacement;

                bubbles[i].skin_distance += -1.0 * bubbles[0].deltadisplacement;

                bubbles[i].displacement += -1.0 * bubbles[0].deltadisplacement;
                bubbles[i].ddisplacement += -1.0 * bubbles[0].deltadisplacement;
            }
        }

        bubbles[0].displacement += -1.0 * bubbles[0].deltadisplacement;
        bubbles[0].deltadisplacement.setZero(d);

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
        bubbles[0].z_flag = 1;

        bubbles[0].com_flag = 0;
    }
}

void update_vol_frac()
{
    double vol_box = 1.0;

    for (int k = 0; k < d; ++k)
    {
        vol_box *= L_box[k];
    }
    vol_frac = bubbles[0].volume / vol_box;
}

void update_p()
{
    update_N();
    
    neighbor_list(type_neighbor_list, 1, nl_para);
    int nl_flag = bubbles[0].nl_flag;

    for (int i = 0; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1) || (i == 0))
        {
            bubbles[i].p = 0.0;
        }
    }

    for (int i = 1; i <= N_init; ++i)
    {
        if ((bubbles[i].flag == 1))
        {
            if ((nl_flag == 0) || (nl_flag == 1))
            {
                for (unsigned int k = 0; k < bubbles[i].non_red_neighbor_list.size(); ++k)
                {
                    int j = bubbles[i].non_red_neighbor_list[k];

                    if (bubbles[j].flag == 1)
                    {
                        bubbles[i].p += 2.0 * gsl_pow_3(bubbles[i].radius + bubbles[j].radius);
                        bubbles[j].p += 2.0 * gsl_pow_3(bubbles[i].radius + bubbles[j].radius);
                        bubbles[0].p += 2.0 * gsl_pow_3(bubbles[i].radius + bubbles[j].radius);
                    }
                }
            }
            else if (nl_flag == -1)
            {
                for (unsigned int k = 0; k < bubbles[i].neighbor_list.size(); ++k)
                {
                    int j = bubbles[i].neighbor_list[k];

                    if (bubbles[j].flag == 1)
                    {
                        bubbles[i].p += gsl_pow_3(bubbles[i].radius + bubbles[j].radius);
                        bubbles[0].p += gsl_pow_3(bubbles[i].radius + bubbles[j].radius);
                    }
                }
            }
        }
    }

    double vol_box = 1.0;

    for (int k = 0; k < d; ++k)
    {
        vol_box *= L_box[k];
    }

    bubbles[0].p = 1.0 + 2.0 * pi / (3.0 * bubbles[0].N * vol_box) * bubbles[0].p;
}