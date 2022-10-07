/*
 * File: BC.h
 * Project: Quenching
 * File Created: Thursday, 28th June 2020 5:28:43 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:33:02 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef BC_H
#define BC_H

#include <headers.h>
#include <globals.h>

VectorXd delta_r(int , VectorXd , VectorXd , VectorXd );

VectorXd r_BC(int , VectorXd , VectorXd );

void fix_wall_particles();

VectorXd fix_wall_BC(const VectorXd &pos_i, const VectorXd &pos_j, const VectorXd &L);

VectorXd fix_wall_BC(const VectorXd &pos_i, const VectorXd &L);

VectorXd periodic_BC(const VectorXd &pos_i, const VectorXd &pos_j, const VectorXd &L);

VectorXd periodic_BC(const VectorXd &pos_i, const VectorXd &L);

void index_correction(VectorXd &index_1, VectorXd &index_2, const VectorXd &N_grid);

#endif
