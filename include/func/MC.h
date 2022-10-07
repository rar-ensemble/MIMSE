/*
 * File: MC.h
 * Project: Quenching
 * File Created: Thursday, 5th December 2020 2:08:58 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:45:32 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef MC_H
#define MC_H

#include <headers.h>
#include <globals.h>

double MC(const int &MC_option, const double &dt, const int &option);

double Metropolis_MC(const int &option);

double particle_swap_MC(const int &option);

double radius_swap_MC(const int &option);

double MC_random_particle(const int &a = 1, const int &b = bubbles[0].N);

void MC_particle_displacement_move(vector<int> &indices);

void MC_particle_swap_move(const int &i, const int &j, vector<particle> &MC_list);

void MC_radius_swap_move(const int &i, const int &j, vector<particle> &MC_list);

double Metropolis_check(vector<int> &indices, vector<particle> &MC_list);

void tolerances_set_MC(const int &MC_option);

#endif