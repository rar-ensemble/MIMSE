/*
 * File: phyrelax.cpp
 * Project: Quenching
 * File Created: Tuesday, 24th July 2020 6:20:10 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:42:12 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#include <headers.h>

#include <globals.h>
#include <functions.h>

void physical_relaxation(const int &dynamics_option, const int &MD_MC_option, const double &dt)
{
    update_force_U(type_force, type_bias_force);

    min_counter_SUCCESSFUL = update_pos_vel(dynamics_option, MD_MC_option, dt);
    min_counter = min_counter_SUCCESSFUL;

    bubbles[0].state_flag = 0;
}