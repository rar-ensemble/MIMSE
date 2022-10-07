/*
 * File: sysint.h
 * Project: Quenching
 * File Created: Friday, 29th June 2020 12:00:19 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 12:53:08 am
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef SYSINT_H
#define SYSINT_H

#include <headers.h>
#include <globals.h>
#include <functions.h>

VectorXd force(const int &force_option = type_force, const int &i = 0, const int &j = 0, double &para_1_force = force_para_1, double &para_2_force = force_para_2, const double &r_cut_off = r_force_cut_off);

VectorXd soft_sphere_force(const int & , const int & , double &alpha = force_para_1, const double &r_cut_off = r_force_cut_off);

VectorXd hard_sphere_force(const int & , const int & , double &alpha = force_para_1, const double &r_cut_off = r_force_cut_off);

void cons_force_physicality(const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, const int &i = 0, const int &j = 0, const VectorXd &r_ij = VectorXd::Zero(d), const double &U = bubbles[0].U);

VectorXd LJ_force(const int &i, const int &j, double &m = force_para_1, double &n = force_para_2, const double &r_LJ_cut_off = r_force_cut_off);

VectorXd LJ_quadratic_cut_off_force(const int &i, const int &j, double &m = force_para_1, double &n = force_para_2, const double &r_LJ_cut_off = r_force_cut_off);

void LJ_force_physicality(const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, const int &i = 0, const int &j = 0, const VectorXd &r_ij = VectorXd::Zero(d), const double &U = bubbles[0].U);

double param_epsilon_ij(const int &epsilon_option = type_epsilon, const int &i = 0, const int &j = 0);

double param_sigma_ij(const int &sigma_option = type_sigma, const int &i = 0, const int &j = 0);

VectorXd bias_force(const int & , const int & );

void neighbor_criteria(const int & , const int & , const double &nl_skin = 0.0, const double &r_cut_off = r_force_cut_off);

MatrixXd stress_tensor(const int &, const int &, const int &);

MatrixXd Hessian_ij(const int &force_option = type_force, const int &i = 0, const int &j = 0, double &para_1_force = force_para_1, double &para_2_force = force_para_2, const double &r_cut_off = r_force_cut_off);

#endif