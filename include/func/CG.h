/*
 * File: CG.h
 * Project: func
 * File Created: Monday, 9th August 2021 6:53:10 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 7:57:15 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef CG_H
#define CG_H

#include <headers.h>
#include <globals.h>

void CG(const int &dynamics_option = type_dynamics_r, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer);

#endif