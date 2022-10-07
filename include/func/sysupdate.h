/*
 * File: sysupdate.h
 * Project: Quenching
 * File Created: Friday, 29th June 2020 11:34:51 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 1:10:59 am
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef SYSUPDATE_H
#define SYSUPDATE_H

#include <headers.h>
#include <globals.h>

double update_pos_vel(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = 0.0, double &f_norm = REF, double &v_norm = REF, double &P = REF);

void update_force_U(const int &force_option = type_force, const int &bias_force_option = type_bias_force);

void update_force_U_i(const int &, const int &, const int &nl_flag = bubbles[0].nl_flag);

void update_bias_force_U_i(const int &, const int &, const int &bl_flag = bubbles[0].bl_flag);

void update_force_correction();

void update_R();

void update_N();

void update_Tau(const int &force_option = type_force);

void update_Hessian(const int &force_option = type_force, double &para_1_force = force_para_1, double &para_2_force = force_para_2, const double &r_cut_off = r_force_cut_off);

void update_min_position();

void update_r_ij();

void update_Hessian_eigenvalues(int num_of_eigenvalues = d * bubbles[0].N, const int &eig_option = 1, const bool &print_eigvec = print_eigenvector);

void update_kbT(const int &system_counter = counter_system, const double &system_time = time_system);

void update_apollonian_order_parameter();

void update_Z();

void mass_conservation();

void update_center_of_mass();

void update_vol_frac();

void update_p();

#endif