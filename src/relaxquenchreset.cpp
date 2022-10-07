/*
 * File: relaxquenchreset.cpp
 * Project: Quenching
 * File Created: Thursday, 20th June 2020 2:22:07 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:43:49 pm
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

void relaxation_quench(const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt)
{
    relaxation(relaxation_option, dynamics_option, MD_MC_option, dt);

    bubbles[0].state_flag -= 1;
}

void relaxation_reset_system(const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt)
{    
    relaxation(relaxation_option, dynamics_option, MD_MC_option, dt);
}