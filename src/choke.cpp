/*
 * File: choke.cpp
 * Project: Quenching
 * File Created: Sunday, 1st July 2020 7:22:29 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:03:11 pm
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

void choke_check(const int &i, const int &j, const VectorXd &r_ij, const double &sigma_i, const double &sigma_j, const double &sigma_ij, const double &r_cut_off, const VectorXd &L)
{
    if ((2.0 * r_cut_off * sigma_ij >= L.minCoeff()) || (r_cut_off * sigma_i >= L.minCoeff()) || (r_cut_off * sigma_j >= L.minCoeff()))
    {
        VectorXd delta_r_pbc_1 = L.array() * (L - (2.0 * r_ij)).array();
        VectorXd delta_r_pbc_2 = L.array() * (L + (2.0 * r_ij)).array();
        double delta_r_pbc = min(delta_r_pbc_1.minCoeff(), delta_r_pbc_2.minCoeff());
        double r_ij_norm_2 = r_ij.squaredNorm();
        double r_ij_pbc_norm_2 = r_ij.squaredNorm() + delta_r_pbc;

        if (((r_ij_norm_2 <= gsl_pow_2(r_cut_off * sigma_ij)) && (r_ij_pbc_norm_2 <= gsl_pow_2(r_cut_off * sigma_ij))) || (r_cut_off * sigma_i >= L.minCoeff()) || (r_cut_off * sigma_j >= L.minCoeff())) // || (2.0 * sigma_ij >= L.minCoeff())) 
        {
            cout << "System choked for particles: " << i << " and " << j << '\n'
                 << "program run stopped \n";
            string filename = "choke.lammpstrj";

            create_lammpstrj(filename, counter_system, time_system);

            time_final = clock();
            run_time = double(time_final - time_initial) / CLOCKS_PER_SEC;

            cout << "time taken :" << run_time << " seconds";

            exit(0);
        }
    }
}