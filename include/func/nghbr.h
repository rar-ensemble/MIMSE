/*
 * File: nghbr.h
 * Project: Quenching
 * File Created: Thursday, 28th June 2020 11:24:47 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:30:14 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef NGBHR_H
#define NGBHR_H

#include <headers.h>
#include <globals.h>

void neighbor_list(const int &neighbor_list_option = type_neighbor_list, const int &nl_dest_flag = 2, const double &nl_skin = 0.0);

void exact_neighbor_list();

#endif