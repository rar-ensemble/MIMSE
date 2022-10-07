/*
 * File: MD.h
 * Project: Quenching
 * File Created: Friday, 29th June 2020 1:06:34 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:56:45 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef MD_H
#define MD_H

#include <headers.h>
#include <globals.h>

double MD(const int &MD_option, const double &dt, const int &option, double &f_norm = REF, double &v_norm = REF, double &P = REF);

double MD(const int &MD_option, const int &i, const double &dt, const int &option, double &f_norm = REF, double &v_norm = REF, double &P = REF);

double simple_verlet_MD(const int &i, const double &dt, const int &option, double &f_norm = REF, double &v_norm = REF, double &P = REF);

double explicit_Euler_LD(const int &i, const double &dt, const int &option, double &f_norm = REF);

double RK_2_LD(const int &i, const double &dt, const int &option);

double explicit_Euler_thermal_LD(const int &i, const double &dt, const int &option);

double RK_2_thermal_LD(const int &i, const double &dt, const int &option);

double simple_verlet_thermal_NH_MD(const int &i, const double &dt, const int &option, double &f_norm = REF, double &v_norm = REF, double &P = REF);

double simple_verlet_thermal_LD(const int &i, const double &dt, const int &option, double &f_norm = REF, double &v_norm = REF, double &P = REF);

void tolerances_set_MD(const int &MD_option, const double &dt);

#endif