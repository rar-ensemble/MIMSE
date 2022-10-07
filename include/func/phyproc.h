/*
 * File: phyproc.h
 * Project: Quenching
 * File Created: Saturday, 30th June 2020 1:59:25 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:37:40 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef PHYPROC_H
#define PHYPROC_H

#include <headers.h>
#include <globals.h>

void quenching(const int &quenching_option = type_quenching, const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer, const int &option = 1, const double &control_variable_final = control_variable_final);

void relaxation(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = 0.002);

double dynamics(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer, const int &option = 1, double &f_norm = REF, double &v_norm = REF, double &P = REF);

void tolerances_set(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

#endif