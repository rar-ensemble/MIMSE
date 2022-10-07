/*
 * File: N_ngbhr.h
 * Project: Quenching
 * File Created: Tuesday, 3rd July 2020 12:42:16 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:31:44 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef NNGBHR_H
#define NNGBHR_H

#include <headers.h>
#include <globals.h>

void N_neighbor(const double &nl_skin = 0.0, const int &nl_dest_flag = 2, const double &r_cut_off = r_force_cut_off);

#endif