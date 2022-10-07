/*
 * File: choke.h
 * Project: Quenching
 * File Created: Sunday, 1st July 2020 11:12:17 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:03:46 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef CHOKE_H
#define CHOKE_H

#include <headers.h>

void choke_check(const int &i, const int &j, const VectorXd &r_ij, const double &sigma_i, const double &sigma_j, const double &sigma_ij, const double &r_cut_off, const VectorXd &L);

#endif