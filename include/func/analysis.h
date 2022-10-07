/*
 * File: analysis.h
 * Project: Quenching
 * File Created: Saturday, 30th June 2020 2:01:16 am
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 10:41:36 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <headers.h>
#include <globals.h>

using namespace std;

void output_U(string filename = U_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, const double &min_count_SUCCESSFUL = min_counter_SUCCESSFUL, const int &min_count = min_counter);

void output_Z(string filename = Z_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_R(string filename = R_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_s_contour(string filename = s_contour_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_displacement(string filename = displacement_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_order_parameter(string filename = order_parameter_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_r_ij_parameter(string filename = r_ij_parameter_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_r_ij(string filename = r_ij_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, int n = 3);

void output_U_i(string filename = U_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_Z_i(string filename = Z_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_R_i(string filename = R_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_s_i(string filename = s_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_contour_i(string filename = contour_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_ds_i(string filename = ds_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_dconotur_i(string filename = dcontour_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_displacement_i(string filename = displacement_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_ddisplacement_i(string filename = ddisplacement_i_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_val(double , string );

void output_avalanche(int avalanche_counter, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_structure(string filename = structure_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, int force_option = type_force);

void output_eigenvalue(string filename = eigenvalues_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_eigenvector(string filename = eigenvectors_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_structure_factor(string filename = structure_factor_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_tolerances(string filename = tolerances_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1);

void output_displacement_vector(string filename = displacement_vector_filename, const int &system_counter = counter_system, const double &system_time = time_system, const int &phyproc_counter = -1, const double &phyproc_time = -1, const unsigned int &j = bias_counter[0] - 1, const double &displacement_r = unit_vector_r[unit_vector_flag - 1]);

#endif