/*
 * File: mathf.h
 * Project: Quenching
 * File Created: Wednesday, 27th June 2020 10:22:06 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:31:18 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef MATHF_H
#define MATHF_H

#include <headers.h>
#include <globals.h>

void uniform_real_dist(const double &, const double &);

void real_gaussian_dist(const double &, const double &);

void gaussian_dist(const double &, const double &);

void mod_gaussian_dist(const double &, const double &);

void real_weibull_dist(const double &, const double &);

void weibull_dist(const double &, const double &);

double weibull_func(const double &, const double &, const double &, const double &, int);

void bidisperse_dist(const double & , const double & );

void monodisperse_dist(const double & );

void binary_dist(const double &a, const double &b, const double &x1);

void truncated_power_law(const double & , const double & );

void swap_power_law(double & , double & );

void random_position_init();

void indexing();

void mass_initial(const int &);

double mean_radius();

double mean_radius2();

double mean_radius3();

double max_radius();

double stddev_radius();

vector<int> largest_radius(const int &n = 3);

void init_run();

vector<double> cardano_real_cubic_equation_solver(const double &p, const double &q);

template <typename T> int sgn(T val);

#endif