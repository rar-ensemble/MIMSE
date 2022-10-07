/*
 * File: metadynamics.h
 * Project: Quenching
 * File Created: Wednesday, 26th June 2020 8:43:36 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:27:28 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include <headers.h>
#include <globals.h>

void energy_landscape_exploration(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

void minimum_exploration(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

void basin_exploration(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

void metadyanmics_relaxation(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer, const int &option = 1);

void metadynamics_reset_system(int &Z_cut_off = bias_Z_cut_off);

int minima_check();

void bias_list(const int &bl_dest_flag = 1, const double &bl_skin = 0.0, const double &bl_end = INFINITY);

#endif