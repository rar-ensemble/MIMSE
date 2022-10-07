/*
 * File: init.h
 * Project: Quenching
 * File Created: Saturday, 30th June 2020 1:56:38 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:07:28 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef INIT_H
#define INIT_H

#include <headers.h>
#include <globals.h>

void initialization();

void size_distribution(const int &distribution_option = type_distribution);

void position_initialization(const int &initialization_option = type_initialization);

#endif