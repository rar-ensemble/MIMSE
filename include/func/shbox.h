/*
 * File: shbox.h
 * Project: Quenching
 * File Created: Friday, 19th April 2020 7:15:58 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 11:50:10 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef SHBOX_H
#define SHBOX_H

#include <headers.h>
#include <globals.h>

void shrinking_box_relaxation(const int &relaxation_option = type_relaxation_q, const int &dynamics_option = type_dynamics_q, const int &MD_MC_option = type_MD_MC_q, const double &dt = dt_minimizer, const int &option = 1, const double &control_variable_target = control_variable_final);

void update_scale_position(const VectorXd &L_old = L_box + dL_box, const VectorXd &L = L_box);

#endif