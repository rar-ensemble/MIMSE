/*
 * File: sysint.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:16:56 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 12:56:17 am
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

VectorXd force(const int &force_option, const int &i, const int &j, double &para_1_force, double &para_2_force, const double &r_cut_off)
{
    switch (force_option)
    {
    case -1:
        if (int(para_1_force) == para_1_force)
        {
            para_1_force = int(para_1_force);
        }
        if (int(para_2_force) == para_2_force)
        {
            para_2_force = int(para_2_force);
        }
        return VectorXd::Zero(2 * d + 1);
        break;

    case 1:
        return soft_sphere_force(i, j, para_1_force, r_cut_off);
        break;

    case 2:
        return hard_sphere_force(i, j, para_1_force, r_cut_off);
        break;

    case 3:
        return LJ_force(i, j, para_1_force, para_2_force, r_cut_off);
        break;

    case 4:
        return LJ_quadratic_cut_off_force(i, j, para_1_force, para_2_force, r_cut_off);
        break;

    case 5:
        return bias_force(i, j);
        break;

    default:
        return soft_sphere_force(i, j, para_1_force, r_cut_off);
        break;
    }
}

VectorXd soft_sphere_force(const int &i, const int &j, double &alpha, const double &r_cut_off)
{
    if (i != 0)
    {
        VectorXd force_U_r_ij = VectorXd::Zero(2 * d + 1);
        VectorXd force_r = VectorXd::Zero(d);
        VectorXd U = VectorXd::Zero(1);

        double sigma_i = param_sigma_ij(type_sigma, i);
        double sigma_j = param_sigma_ij(type_sigma, j);

        VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);
        double r_ij_norm = r_ij.norm();

        double epsilon_ij = param_epsilon_ij(type_epsilon, i, j);
        double sigma_ij = param_sigma_ij(type_sigma, i, j);

        if (r_ij_norm < r_cut_off * sigma_ij)
        {
            if (alpha == 2)
            {
                force_r = 1.0 * epsilon_ij / sigma_ij * (1.0 - (r_ij_norm / sigma_ij)) * r_ij / r_ij_norm;
                U[0] = epsilon_ij / alpha * gsl_pow_2((1.0 - (r_ij_norm / sigma_ij)));
                force_U_r_ij << force_r, U, r_ij;
            }
            else if (alpha == int(alpha))
            {
                force_r = 1.0 * epsilon_ij / sigma_ij * gsl_pow_int((1.0 - (r_ij_norm / sigma_ij)), (alpha - 1)) * r_ij / r_ij_norm;
                U[0] = epsilon_ij / alpha * gsl_pow_int((1.0 - (r_ij_norm / sigma_ij)), alpha);
                force_U_r_ij << force_r, U, r_ij;
            }
        }

        if (r_ij_norm <= r_cut_off * sigma_ij)
        {
            choke_check(i, j, r_ij, sigma_i, sigma_j, sigma_ij, r_cut_off, L_box);
            // cons_force_physicality(counter_system, time_system, counter_phyproc, time_phyproc, i, j, r_ij, U[0]);
        }

        return force_U_r_ij;
    }
    else
    {
        double epsilon_ij = param_epsilon_ij(type_epsilon, i, j);
        if (alpha == 2)
        {
            if (sigma_tol_0 == -1)
            {
                sigma_tol = sqrt(eff_U_0 * U_tol * alpha / epsilon_ij);
            }
            else if (sigma_tol_0 == -2)
            {
                sigma_tol = sqrt(eff_U_0 * U_tol * alpha / epsilon_ij);
                sigma_tol = max(sigma_tol_0, sigma_tol);
            }
            else
            {
                sigma_tol = sigma_tol_0;
            }
        }
        else if (alpha == int(alpha))
        {
            if (sigma_tol_0 == -1)
            {
                sigma_tol = pow(eff_U_0 * U_tol * alpha / epsilon_ij, 1.0 / alpha);
            }
            else if (sigma_tol_0 == -2)
            {
                sigma_tol = pow(U_0 * U_tol * alpha / epsilon_ij, 1.0 / alpha);
                sigma_tol = max(sigma_tol_0, sigma_tol);
            }
            else
            {
                sigma_tol = sigma_tol_0;
            }
        }

        return VectorXd::Zero(d);
    }
}

VectorXd hard_sphere_force(const int &i, const int &j, double &alpha, const double &r_cut_off)
{
    VectorXd force_U_r_ij = VectorXd::Zero(2 * d + 1);
    VectorXd force_r = VectorXd::Zero(d);
    VectorXd U = VectorXd::Zero(1);

    VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);
    double r_ij_norm = r_ij.norm();

    double epsilon_ij = param_epsilon_ij(type_epsilon, i, j);
    double sigma_ij = param_sigma_ij(type_sigma, i, j);

    if (r_ij_norm < r_cut_off * sigma_ij)
    {
        if (alpha == 2)
        {
            U[0] = epsilon_ij / alpha * gsl_pow_2((1.0 - (r_ij_norm / sigma_ij)));
        }
        else if (alpha == int(alpha))
        {
            U[0] = epsilon_ij / alpha * gsl_pow_int((1.0 - (r_ij_norm / sigma_ij)), alpha);
        }
    }

    if (U[0] / eff_U_0 < U_tol)
    {
        U[0] = 0;
    }
    else
    {
        U[0] = INFINITY;
    }

    force_U_r_ij << force_r, U, r_ij;

    return force_U_r_ij;
}

void cons_force_physicality(const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, int i, int j, const VectorXd &r_ij, const double &U)
{
    double r_ij_norm = r_ij.norm();

    double sigma_ij = param_sigma_ij(type_sigma, i, j);

    if ((1.0 - (r_ij_norm / sigma_ij)) > 1 / 4.0)
    {
        analysis_files.open((output_PATH + unphysical_filename).c_str(), ios::out | ios::app);
        if (analysis_files.tellp() == 0)
        {
            analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag i j U";
        }

        analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                       << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << i << " " << j << " " << U;

        analysis_files.close();
    }
}

VectorXd LJ_force(const int &i, const int &j, double &m, double &n, const double &r_LJ_cut_off)
{
    VectorXd force_U_r_ij = VectorXd::Zero(2 * d + 1);
    VectorXd force_r = VectorXd::Zero(d);
    VectorXd U = VectorXd::Zero(1);

    static const double k = n / (m - n) * pow(n / m, m / (n - m));

    double epsilon_ij = param_epsilon_ij(type_epsilon, i, j);
    double sigma_ij = param_sigma_ij(type_sigma, i, j);

    VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);
    double r_ij_norm = r_ij.norm();

    double sigma_r_ij = sigma_ij / r_ij_norm;

    double F_max = 1e3 * k * epsilon_ij / sigma_ij * (m - n);

    int unphysical_flag = 0;

    if (r_ij_norm < r_LJ_cut_off * sigma_ij)
    {
        if ((m == 12) && (n == 6))
        {
            double sigma_r_ij_n = gsl_pow_6(sigma_r_ij);
            double sigma_r_ij_m = gsl_pow_2(sigma_r_ij_n);

            double r_LJ_cut_off_n = gsl_pow_6(1.0 / r_LJ_cut_off);
            double r_LJ_cut_off_m = gsl_pow_2(r_LJ_cut_off_n);

            double force_coeff = k * epsilon_ij / sigma_ij * (m * sigma_r_ij_m * sigma_r_ij - n * sigma_r_ij_n * sigma_r_ij);
            if (force_coeff >= F_max)
            { 
                force_coeff = F_max;
                unphysical_flag = 1;
            }
            force_r = force_coeff * r_ij / r_ij_norm;
            U[0] = k * epsilon_ij * ((sigma_r_ij_m - sigma_r_ij_n) - (r_LJ_cut_off_m - r_LJ_cut_off_n));
            force_U_r_ij << force_r, U, r_ij;
        }
        else
        {
            double force_coeff = k * epsilon_ij / sigma_ij * (m * pow(sigma_ij / r_ij_norm, m + 1) - n * pow(sigma_ij / r_ij_norm, n + 1));
            if (force_coeff >= F_max)
            { 
                force_coeff = F_max;
                unphysical_flag = 1;
            }
            force_r = force_coeff * r_ij / r_ij_norm;
            U[0] = k * epsilon_ij * (pow(sigma_ij / r_ij_norm, m) - pow(sigma_ij / r_ij_norm, n) - (pow(r_LJ_cut_off, -m) - pow(r_LJ_cut_off, -n)));
            force_U_r_ij << force_r, U, r_ij;
        }
    }

    if (unphysical_flag == 1)
    {
        LJ_force_physicality(counter_system, time_system, counter_phyproc, time_phyproc, i, j, r_ij, U[0]);
    }

    return force_U_r_ij;
}

VectorXd LJ_quadratic_cut_off_force(const int &i, const int &j, double &m, double &n, const double &r_LJ_cut_off)
{
    VectorXd force_U_r_ij = VectorXd::Zero(2 * d + 1);
    VectorXd force_r = VectorXd::Zero(d);
    VectorXd U = VectorXd::Zero(1);

    static const double k = n / (m - n) * pow(n / m, m / (n - m));

    double epsilon_ij = param_epsilon_ij(type_epsilon, i, j);
    double sigma_ij = param_sigma_ij(type_sigma, i, j);

    VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);
    double r_ij_norm = r_ij.norm();

    double sigma_r_ij = sigma_ij / r_ij_norm;

    double F_max = 1e2 * k * epsilon_ij / sigma_ij * (m - n);

    int unphysical_flag = 0;

    if (r_ij_norm < r_LJ_cut_off * sigma_ij)
    {
        if ((m == 12) && (n == 6))
        {
            double sigma_r_ij_n = gsl_pow_6(sigma_r_ij);
            double sigma_r_ij_m = gsl_pow_2(sigma_r_ij_n);

            double r_LJ_cut_off_n = gsl_pow_6(1.0 / r_LJ_cut_off);
            double r_LJ_cut_off_m = gsl_pow_2(r_LJ_cut_off_n);

            double force_coeff = k * epsilon_ij / sigma_ij * ((m * sigma_r_ij_m * sigma_r_ij - n * sigma_r_ij_n * sigma_r_ij) - (m * r_LJ_cut_off_m / r_LJ_cut_off - n * r_LJ_cut_off_n / r_LJ_cut_off) * (1.0 / (r_LJ_cut_off * sigma_r_ij)));
            if (force_coeff >= F_max)
            { 
                force_coeff = F_max;
                unphysical_flag = 1;
            }
            force_r = force_coeff * r_ij / r_ij_norm;
            U[0] = k * epsilon_ij * ((sigma_r_ij_m - sigma_r_ij_n) + (m / 2 * r_LJ_cut_off_m - n / 2 * r_LJ_cut_off_n) * gsl_pow_2(1.0 / (r_LJ_cut_off * sigma_r_ij)) - ((m / 2 + 1) * r_LJ_cut_off_m - (n / 2 + 1) * r_LJ_cut_off_n));
            force_U_r_ij << force_r, U, r_ij;
        }
        else
        {
            double force_coeff = k * epsilon_ij / sigma_ij * ((m * pow(sigma_ij / r_ij_norm, m + 1) - n * pow(sigma_ij / r_ij_norm, n + 1)) - (m * pow(r_LJ_cut_off, -(m + 1)) - n * pow(r_LJ_cut_off, -(n + 1))) * r_ij_norm / (r_LJ_cut_off * sigma_ij));
            if (force_coeff >= F_max)
            { 
                force_coeff = F_max;
                unphysical_flag = 1;
            }
            force_r = force_coeff * r_ij / r_ij_norm;
            U[0] = k * epsilon_ij * ((pow(sigma_ij / r_ij_norm, m) - pow(sigma_ij / r_ij_norm, n)) + (m / 2 * pow(r_LJ_cut_off, -m) - n / 2 * pow(r_LJ_cut_off, -n)) * pow(r_ij_norm / (r_LJ_cut_off * sigma_ij), 2) - ((m / 2 + 1) * pow(r_LJ_cut_off, -m) - (n / 2 + 1) * pow(r_LJ_cut_off, -n)));
            force_U_r_ij << force_r, U, r_ij;
        }
    }

    if (unphysical_flag == 1)
    {
        LJ_force_physicality(counter_system, time_system, counter_phyproc, time_phyproc, i, j, r_ij, U[0]);
    }

    return force_U_r_ij;
}

void LJ_force_physicality(const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, const int &i, const int &j, const VectorXd &r_ij, const double &U)
{
    analysis_files.open((output_PATH + unphysical_filename).c_str(), ios::out | ios::app);
    if (analysis_files.tellp() == 0)
    {
        analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag i j U";
    }

    analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                    << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << i << " " << j << " " << U;

    analysis_files.close();
}

double param_epsilon_ij(const int &epsilon_option, const int &i, const int &j)
{
    switch (epsilon_option)
    {
    case 1:
        return 2.0;
        break;

    case 2:
        return 1.0 / (bubbles[i].sigma + bubbles[j].sigma);
        break;

    case 3:
        return bubbles[0].radius / 2.0;
        break;

    case 4:
        return 0;
        break;

    case 5:
        return sqrt(bubbles[i].epsilon * bubbles[j].epsilon);
        break;

    case 6:
        if (j == 0)
        {
            return 1.0;
        }
        else if (bubbles[i].type == bubbles[j].type)
        {
            if (bubbles[i].type == 1)
            {
                return 1.0;
            }
            else if (bubbles[i].type == 2)
            {
                return 0.5;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 1.5;
        }

    default:
        return 1.0;
        break;
    }
}

double param_sigma_ij(const int &sigma_option, const int &i, const int &j)
{
    switch (sigma_option)
    {
    case 1:
        if (j == 0)
        {
            return 2.0 * bubbles[i].radius;
        }
        else
        {
            return bubbles[i].radius + bubbles[j].radius;
        }
        break;

    case 2:
        if (j == 0)
        {
            if (i == 0)
            {
                return 1.0;
            }
            else if (bubbles[i].type == 1)
            {
                return 1.0;
            }
            else if (bubbles[i].type == 2)
            {
                return 0.88;
            }
            else
            {
                return 0;
            }
        }
        else if (bubbles[i].type == bubbles[j].type)
        {
            if (bubbles[i].type == 1)
            {
                return 1.0;
            }
            else if (bubbles[i].type == 2)
            {
                return 0.88;
            }
            else
            {
                return 0;
            }
        }
        else if (bubbles[i].type != bubbles[j].type)
        {
            return 0.8;
        }
        break;

    default:
        return 1.0;
        break;
    }
}

VectorXd bias_force(const int &i, const int &j)
{
    VectorXd force_U_i = VectorXd::Zero(d + 1);
    VectorXd force_r = VectorXd::Zero(d);
    VectorXd U = VectorXd::Zero(1);

    double biaspointfactor = bubbles[0].bias_U_position_list[j][1];
    double wellpointfactor = 0.0;

    if (biaspointfactor < gsl_pow_2(max(bias_U_sigma, ((bias_counter[2] == 1) * (well_U_sigma + well_cut_off)))))
    {
        if (i == 0)
        {
            if (biaspointfactor < gsl_pow_2(bias_U_sigma))
            {
                biaspointfactor = bubbles[0].bias_U_position_list[j][1] / gsl_pow_2(bias_U_sigma) - 1.0;

                U[0] = bias_U_0 * gsl_pow_2(biaspointfactor);
            }

            if (bias_counter[2] == 1)
            {
                wellpointfactor = sqrt(bubbles[0].bias_U_position_list[j][1]);

                if ((wellpointfactor > well_cut_off) && (wellpointfactor < (well_U_sigma + well_cut_off)))
                {
                    U[0] += well_U_0 * gsl_pow_2((wellpointfactor - well_cut_off) / (well_U_sigma));
                }
            }
        }
        else
        {
            if (bubbles[i].bias_U_position_list[j][0] != -1)
            {
                if (biaspointfactor < gsl_pow_2(bias_U_sigma))
                {
                    biaspointfactor = bubbles[0].bias_U_position_list[j][1] / gsl_pow_2(bias_U_sigma) - 1.0;
                    U[0] = bias_U_0 * gsl_pow_2(biaspointfactor);
                    force_r = -4.0 * bias_U_0 / gsl_pow_2(bias_U_sigma) * biaspointfactor * periodic_BC(bubbles[i].position, bubbles[i].bias_U_position_list[j], L_box);
                }

                if (bias_counter[2] == 1)
                {
                    wellpointfactor = sqrt(bubbles[0].bias_U_position_list[j][1]);

                    if ((wellpointfactor > well_cut_off) && (wellpointfactor < (well_U_sigma + well_cut_off)))
                    {
                        U[0] += well_U_0 * gsl_pow_2((wellpointfactor - well_cut_off) / (well_U_sigma));
                        force_r += (-2.0 * well_U_0 / well_U_sigma * ((1.0 - (well_cut_off / wellpointfactor)) / well_U_sigma) * periodic_BC(bubbles[i].position, bubbles[i].bias_U_position_list[j], L_box));
                    }
                }
            }
        }
        force_U_i << force_r, U;
    }

    return force_U_i;
}

void neighbor_criteria(const int &i, const int &j, const double &nl_skin, const double &r_cut_off)
{
    if ((i != j) && (bubbles[i].flag == 1) && (bubbles[j].flag == 1))
    {
        double sigma_i = param_sigma_ij(type_sigma, i);
        double sigma_j = param_sigma_ij(type_sigma, j);

        VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);
        double r_ij_norm_2 = r_ij.squaredNorm();

        double sigma_ij = param_sigma_ij(type_sigma, i, j);

        if (nl_skin != -1.0)
        {
            if (((sigma_i > sigma_j) || ((sigma_i == sigma_j) && (i > j))) && (r_ij_norm_2 <= gsl_pow_2(r_cut_off * sigma_ij + sigma_tol + nl_skin)))
            {
                if (nl_skin == 0.0)
                {
                    bubbles[i].non_red_neighbor_list.push_back(j);

                    bubbles[i].neighbor_list.push_back(j);
                    bubbles[j].neighbor_list.push_back(i);
                }

                bubbles[i].non_red_com_neighbor_list.push_back(j);

                bubbles[i].com_neighbor_list.push_back(j);
                bubbles[j].com_neighbor_list.push_back(i);
            }
        }
        else if (nl_skin == -1.0)
        {
            if (((sigma_i > sigma_j) || ((sigma_i == sigma_j) && (i > j))) && (r_ij_norm_2 <= gsl_pow_2(r_cut_off * sigma_ij + sigma_tol)))
            {
                bubbles[i].non_red_neighbor_list.push_back(j);

                bubbles[i].neighbor_list.push_back(j);
                bubbles[j].neighbor_list.push_back(i);
            }
            else if (((sigma_i < sigma_j) || ((sigma_i == sigma_j) && (i < j))) && (r_ij_norm_2 <= gsl_pow_2(r_cut_off * sigma_ij + sigma_tol)))
            {
                bubbles[j].non_red_neighbor_list.push_back(i);

                bubbles[j].neighbor_list.push_back(i);
                bubbles[i].neighbor_list.push_back(j);
            }
        }
    }
}