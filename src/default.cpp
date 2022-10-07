/*
 * File: deafult.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 5:05:12 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Wednesday, 26th January 2022 9:46:26 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */

#include <headers.h>
#include <classes.h>
#include <globals.h>

using namespace std;

double sqrt_d = sqrt(d);

clock_t time_initial = 0;
clock_t time_final = 0;
double run_time = 0;

int counter_system = 0;
double time_system = 0.0;

int counter_phyproc = 0;
double time_phyproc = 0.0;

//n_flag = 1 = need to change, = 0 = all good
int N_init = 2000;

vector<particle> bubbles(N_init + 1);
vector<particle> bubbles_athermal;
vector<particle> bubbles_thermal;

VectorXd L_box = VectorXd::Zero(d);
VectorXd dL_box = VectorXd::Zero(d);

double vol_frac = 0.75;
double dvol_frac_0 = 0;

double b = 0.00001;

int type_kbT = 1.0;
double kbT = 1.0;
double kbTi = 1.0;
double kbTf = 0.0;
double dkbT = 0.0;

// bubbles[0].f_flag = 2 = forces needed and wrong values assigned, = 1 = forces zero and correct values needed, = 0 = all good
int type_force = 1;
int type_epsilon = 1;
int type_sigma = 1;
double force_para_1 = 2.0;
double force_para_2 = 2.0;
double r_force_cut_off = 2.0;

double U_0 = 1.0;
double F_0 = 1.0;

double mass_init = 0.0;
double mass_N = 0.0;
int type_distribution = 2;
double para_1 = 1.0; //1.0; //1.0; //2.1;//1;//0.72;//1;//5;//2;//1.0;
double para_2 = 0.5; //1.13;//0.4;//0.4;//0.5;//1;//2/sqrt(pi);//0.5;
double para_3 = -1;

int type_mass = -1;

int r_flag = -1; // r_flag = 2 = change needed, = 1 = zero values, = 0 = all good
int type_initialization = 1;
SparseMatrix<double> r_ij_system; // (d * N_init, N_init);
vector<Triplet<double> > r_ij_tripletList;
double r_ij_parameter_system = 0.0;

double x_ij_system = INFINITY;

int type_BC = 2;

VectorXd N_grid_sf = VectorXd::Constant(d, 1);
VectorXd N_grid_sf_max = VectorXd::Constant(d, 100);
vector<vector<vector<VectorXd> > > structure_factor_position;
vector<vector<vector<int> > > structure_factor_value;
vector<vector<vector<double> > > structure_factor_value_weighted;
vector<vector<vector<int> > > I_value;

int type_pre_start = 1;

int type_quenching = 1;
int min_counter = 0;
double min_counter_SUCCESSFUL = 0;
int type_reset_system = 1;

int type_relaxation_q = 1;
int type_relaxation_r = 1;

int type_controller = 1;
double control_variable_final = vol_frac;

// bubbles[0].bl_flag = 3 approx list needed, 2 = values needed, 1 = exact list needed, 0 = all good
vector<double> bl_para(2, 0);

int type_bias_force = 5;
double bias_U_0 = 0.0;
double bias_U_sigma = 1.0;
int bias_Z_cut_off = d + 1;
double well_U_0 = 0.0;
double well_U_sigma = 1.0;
double well_cut_off = 1.0;

double bias_F_0 = bias_U_0 / bias_U_sigma;

vector<int> bias_counter(3, -4);

double eff_U_0 = (bias_counter[1] >= -1) ? min(U_0, bias_U_0) : U_0;
double eff_F_0 = (bias_counter[1] >= -1) ? min(F_0, bias_F_0) : F_0;

int type_dynamics_q = 1;
int type_dynamics_r = 1;

int type_MD_MC_q = 1;
int type_MD_MC_r = 1;

double p_swap = 0.1;
double delta_r_MC = bubbles[0].radius;

// bubbles[0].nl_flag = 3 = both lists needed, = 2 = exact list needed, 1 = exact list has del particles, 0 = all good
int type_neighbor_list = 1;
double nl_para = 1.0;

// bubbles[0].z_flag = 1 = change needed, = 0 = all good

double apollonian_order_parameter_system = 0;

// bubbles[0].t_flag = 2 = change needed, 1 = z val, = 0 = all good

double dt_minimizer = 0.01;
double dt_reset_system = 0.01;

int NUMBER_OF_QUENCHES = 1;
int NUMBER_OF_EQB_STEPS = 0;

double U_tol_0 = 1e-10;
double F_tol_0 = 1e-10;

double U_tol = 1e-10;
double F_tol = 1e-10;
double cv_tol = 0.01;

double sigma_tol_0 = 0;
double sigma_tol = 0;

int h_flag = -1; // h_flag = 3 = change needed, = 2 = zero values, 1 = good with Hessian,  = 0 = all good
SparseMatrix<double> Hessian_system; // (d * N_init, d * N_init);
vector<Triplet<double> > Hessian_tripletList;
VectorXd Hessianeigenvalues_system; // = VectorXd::Zero(d * N_init);
MatrixXd Hessianeigenvectors_system; // = MatrixXd::Zero(d * N_init, d * N_init);
int type_eigenvalue_solver = 1;
bool print_eigenvector = false;
int number_of_eigenvalues = d * N_init;
int max_iters_eig = 1000;
double eigs_tol = 1e-12;

int type_force_flag = 1;
int particle_fix_flag = 0;

double random_seed_size_dist = 1.0;
double random_seed_pos_init = 1.0;
double random_seed_pos_init_discard = 0;
double random_seed_random_number = 1.0;
default_random_engine generator_1;
default_random_engine generator_2;
default_random_engine generator_3;
default_random_engine generator_4;

string input_PATH = "input/";
string input_filename = "init.input";
string init_filename = "init.lammpstrj";

string run_parameters_filename = "run_config.txt";

fstream analysis_files;
fstream io_files;
int output_mode = 1; // 1 = outside, 2 = inside

int frame_print_frequency = 0;
int utils_print_frequency = 0;
int tau_print_frequency = 0;
int hessian_print_frequency = 0;
int structure_print_frequency = 0;
int structure_factor_print_frequency = 0;
int data_print_frequency = 1;
int individual_data_print_frequency = 0;

double ds_print_frequency = 0;
double dcontour_print_frequency = 0;

int ds_print_counter = 1;
int dcontour_print_counter = 1;

int avalanche_count = 10;

string output_PATH = "output/";
string U_filename = "U.txt";
string Z_filename = "Z.txt";
string R_filename = "R.txt";
string s_contour_filename = "s.txt";
string displacement_filename = "displacement.txt";
string order_parameter_filename = "op.txt";
string r_ij_parameter_filename = "r_ij_p.txt";
string r_ij_filename = "r_ij.txt";
string structure_filename = "structure_";
string U_i_filename = "U_i.txt";
string Z_i_filename = "Z_i.txt";
string R_i_filename = "R_i.txt";
string s_i_filename = "s_i.txt";
string contour_i_filename = "contour_i.txt";
string ds_i_filename = "ds_i.txt";
string dcontour_i_filename = "dconotur_i.txt";
string displacement_i_filename = "displacement_i.txt";
string ddisplacement_i_filename = "ddisplacement_i.txt";
string eigenvalues_filename = "eigenvalue.txt";
string eigenvectors_filename = "eigenvector.txt";
string frame_filename = "config_";
string utils_filename = "config_";
string structure_factor_filename = "struc_";
string tolerances_filename = "tol.txt";
string unphysical_filename = "unphysical.txt";
string displacement_vector_filename = "d_v.txt";

string analysis_PATH = "analysis/input/";
string analysis_filename = "analysis.input";

int unit_vector_flag = -1;
vector<double> unit_vector_r;

double REF = -1;
int iREF = -1;
vector<double> vdREF;
vector<VectorXd> vvREF;