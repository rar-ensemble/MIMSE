/*
 * File: resetsys.h
 * Project: Quenching
 * File Created: Friday, 28th June 2020 4:48:01 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 11:11:24 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef RESETSYS_H
#define RESETSYS_H

#include <headers.h>
#include <globals.h>

void start_system();

void reset_system(const int &reset_system_option = type_reset_system, const int &relaxation_option = type_relaxation_r, const int &dynamics_option = type_dynamics_r, const int &MD_MC_option = type_MD_MC_r, const double &dt = dt_reset_system, const int &distribution_option = type_distribution, const int &initialization_option = type_initialization);

void position_reset_system(const int &initialization_option = type_initialization);

void dispalcement_reset_system(double &push_factor = bias_U_sigma, const double &bias_U = bias_U_0, const int &Z_cut_off = bias_Z_cut_off);

#endif