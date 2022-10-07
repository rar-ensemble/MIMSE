/*
 * File: inout.h
 * Project: Quenching
 * File Created: Sunday, 1st July 2020 5:24:44 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:09:08 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef INOUT_H
#define INOUT_H

#include <headers.h>
#include <globals.h>

void input();

void output();

void file_init();

void create_lammpstrj(string , const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, int option = 0.0);

void output_utils(const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void input_utils();

void output_analysis(int output_flag = output_mode, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, int &avalanchecounter = iREF, int force_option = type_force);

#endif