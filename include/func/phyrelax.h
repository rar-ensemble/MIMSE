/*
 * File: phyrelax.h
 * Project: Quenching
 * File Created: Tuesday, 24th July 2020 6:35:04 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:41:07 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef PHYRELAX_H
#define PHYRELAX_H

#include <headers.h>
#include <globals.h>

void physical_relaxation(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

#endif