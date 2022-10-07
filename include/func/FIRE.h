/*
 * File: FIRE.h
 * Project: Quenching
 * File Created: Friday, 29th June 2020 10:23:36 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:06:15 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef FIRE_H
#define FIRE_H

#include <headers.h>
#include <globals.h>

void FIRE(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = 0.002);

void hybrid_relaxation(const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = 0.002);

#endif