/*
 * File: globals.h
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 5:02:55 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Wednesday, 26th January 2022 9:47:07 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include <headers.h>
#include <classes.h>

#define pi 3.14159265358979
#define d 3

#define N_MAX 10000

#define FIRE_STEP_MAX 1e20
#define SD_STEP_MAX 1e20
#define CG_STEP_MAX 1e20
#define SHBOX_STEP_MAX 1e20
#define METAD_STEP_MAX 1e20
#define TOL_MIN 1e-15
/*#define COARSEN_STEP_MAX 1e20*/

/* extern vector<double> random_seed;
extern int random_seed_counter; */

extern double sqrt_d;

extern int type_force_flag;
extern int particle_fix_flag;

extern double random_seed_size_dist;
extern double random_seed_pos_init;
extern double random_seed_pos_init_discard;
extern double random_seed_random_number;

extern default_random_engine generator_1;
extern default_random_engine generator_2;
extern default_random_engine generator_3;
extern default_random_engine generator_4;

extern int N_init;

extern double mass_init;
extern double mass_N;

extern vector<particle> bubbles;
extern vector<particle> bubbles_athermal;
extern vector<particle> bubbles_thermal;

extern VectorXd L_box;

extern VectorXd dL_box;

extern double vol_frac;

extern double dvol_frac_0;

extern double b;

extern int type_kbT;
extern double kbT;
extern double kbTi;
extern double kbTf;
extern double dkbT;

extern int counter_system;
extern double time_system;

extern int counter_phyproc;
extern double time_phyproc;

extern int NUMBER_OF_QUENCHES;
extern int NUMBER_OF_EQB_STEPS;

extern SparseMatrix<double> Hessian_system;

extern SparseMatrix<double> r_ij_system;

extern VectorXd Hessianeigenvalues_system;
extern MatrixXd Hessianeigenvectors_system;

extern vector<Triplet<double> > Hessian_tripletList;

extern vector<Triplet<double> > r_ij_tripletList;

extern VectorXd N_grid_sf;
extern VectorXd N_grid_sf_max;
extern vector<vector<vector<VectorXd> > > structure_factor_position;
extern vector<vector<vector<int> > > structure_factor_value;
extern vector<vector<vector<double> > > structure_factor_value_weighted;
extern vector<vector<vector<int> > > I_value;

extern double r_ij_parameter_system;

extern double apollonian_order_parameter_system;

extern double x_ij_system;

extern int type_BC;

extern int type_distribution;
extern double para_1;
extern double para_2;
extern double para_3;

extern int type_mass;

extern int type_initialization;

extern int type_pre_start;

extern int type_quenching;
extern int type_reset_system;

extern int type_relaxation_q;
extern int type_relaxation_r;

extern int type_dynamics_q;
extern int type_dynamics_r;

extern int type_MD_MC_q;
extern int type_MD_MC_r;

extern double p_swap;
extern double delta_r_MC;

extern int type_force;
extern int type_bias_force;
extern int type_epsilon;
extern int type_sigma;

extern int type_neighbor_list;
extern double nl_para;

extern vector<double> bl_para;

extern double dt_minimizer;
extern double dt_reset_system;

extern double U_tol_0;
extern double F_tol_0;

extern double U_tol;
extern double F_tol;
extern int min_counter;
extern double min_counter_SUCCESSFUL;

extern double sigma_tol;
extern double sigma_tol_0;

extern vector<int> bias_counter;

extern double force_para_1;
extern double force_para_2;
extern double r_force_cut_off;

extern double F_0;
extern double U_0;

extern double bias_U_0;
extern double bias_U_sigma;
extern int bias_Z_cut_off;
extern double well_U_0;
extern double well_U_sigma;
extern double well_cut_off;

extern double bias_F_0;

extern double eff_U_0;
extern double eff_F_0;

extern int h_flag;
extern int r_flag;

/* extern int job_id; */

extern fstream analysis_files;
extern fstream io_files;

extern string input_filename;
extern string init_filename;

extern string utils_filename;
extern string frame_filename;

extern string analysis_filename;
/* extern string output_foldername; */

extern int frame_print_frequency;
extern int utils_print_frequency;
extern int tau_print_frequency;
extern int hessian_print_frequency;
extern int structure_print_frequency;
extern int structure_factor_print_frequency;
extern int data_print_frequency;
extern int individual_data_print_frequency;

extern double ds_print_frequency;
extern double dcontour_print_frequency;

extern int ds_print_counter;
extern int dcontour_print_counter;

extern int avalanche_count;

extern string U_filename;
extern string Z_filename;
extern string R_filename;
extern string s_contour_filename;
extern string displacement_filename;
extern string order_parameter_filename;
extern string r_ij_parameter_filename;
extern string r_ij_filename;
extern string structure_filename;
extern string U_i_filename;
extern string Z_i_filename;
extern string R_i_filename;
extern string s_i_filename;
extern string contour_i_filename;
extern string ds_i_filename;
extern string dcontour_i_filename;
extern string displacement_i_filename;
extern string ddisplacement_i_filename;
extern string eigenvalues_filename;
extern string eigenvectors_filename;
extern string structure_factor_filename;
extern string tolerances_filename;
extern string unphysical_filename;
extern string displacement_vector_filename;

extern string run_parameters_filename;

extern string input_PATH;
extern string output_PATH;

extern int output_mode; // 1 = outside, 2 = inside

extern string analysis_PATH;

extern double REF;
extern int iREF;
extern vector<double> vdREF;
extern vector<VectorXd> vvREF;

extern clock_t time_initial;
extern clock_t time_final;
extern double run_time;

extern int avalanche_counter;

extern int max_iters_eig;
extern double eigs_tol;
extern int number_of_eigenvalues;
extern int type_eigenvalue_solver;
extern bool print_eigenvector;

extern int type_controller;
extern double control_variable_final;
extern double cv_tol;

extern int unit_vector_flag;
extern vector<double> unit_vector_r;

#endif
