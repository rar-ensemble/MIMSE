/*
 * File: relaxquenchreset.h
 * Project: Quenching
 * File Created: Thursday, 20th June 2020 3:23:08 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:43:10 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef RELAXRESETQUENCH_H
#define RELAXRESETQUENCH_H

#include <headers.h>
#include <globals.h>

void relaxation_quench(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

void relaxation_reset_system(const int &relaxation_option = type_relaxation_r, const int &dynamics_option = type_dynamics_r, const int &MD_MC_option = type_MD_MC_r, const double &dt = dt_reset_system);

#endif