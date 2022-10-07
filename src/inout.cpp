/*
 * File: inout.cpp
 * Project: Quenching
 * File Created: Sunday, 1st July 2020 4:35:51 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 1st February 2022 6:04:27 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#include <headers.h>

#include <globals.h>
#include <classes.h>
#include <functions.h>

void input()
{
  string line;

  double temp;
  double n_size;

  io_files.open((input_PATH + input_filename).c_str(), ios::in);
  io_files.clear();
  io_files.seekg(0, ios::beg);

  io_files >> N_init;
  getline(io_files, line);
  getline(io_files, line);

  for (int k = 0; k < d; ++k)
  {
    io_files >> L_box[k];
  }
  getline(io_files, line);
  for (int k = 0; k < d; ++k)
  {
    io_files >> dL_box[k];
  }
  getline(io_files, line);
  getline(io_files, line);

  io_files >> vol_frac;
  getline(io_files, line);
  io_files >> dvol_frac_0;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> b;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_kbT;
  getline(io_files, line);
  io_files >> kbTi;
  getline(io_files, line);
  io_files >> kbTf;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_force;
  getline(io_files, line);
  io_files >> type_epsilon;
  getline(io_files, line);
  io_files >> type_sigma;
  getline(io_files, line);
  io_files >> force_para_1;
  getline(io_files, line);
  io_files >> force_para_2;
  getline(io_files, line);
  io_files >> r_force_cut_off;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_distribution;
  getline(io_files, line);
  io_files >> para_1;
  getline(io_files, line);
  io_files >> para_2;
  getline(io_files, line);
  io_files >> para_3;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_mass;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_BC;
  getline(io_files, line);
  getline(io_files, line);

  for (int k = 0; k < d; ++k)
  {
    io_files >> N_grid_sf[k];
  }
  getline(io_files, line);
  for (int k = 0; k < d; ++k)
  {
    io_files >> N_grid_sf_max[k];
  }
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_initialization;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_pre_start;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_quenching;
  getline(io_files, line);
  io_files >> type_reset_system;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_relaxation_q;
  getline(io_files, line);
  io_files >> type_relaxation_r;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_controller;
  getline(io_files, line);
  io_files >> control_variable_final;
  getline(io_files, line);
  getline(io_files, line);

  for (int k = 0; k < 2; ++k)
  {
    io_files >> bl_para[k];
  }
  getline(io_files, line);
  getline(io_files, line);

  io_files >> bias_U_0;
  getline(io_files, line);
  io_files >> bias_U_sigma;
  getline(io_files, line);
  io_files >> bias_Z_cut_off;
  getline(io_files, line);
  io_files >> well_U_0;
  getline(io_files, line);
  io_files >> well_U_sigma;
  getline(io_files, line);
  io_files >> well_cut_off;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_dynamics_q;
  getline(io_files, line);
  io_files >> type_dynamics_r;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_MD_MC_q;
  getline(io_files, line);
  io_files >> type_MD_MC_r;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> p_swap;
  getline(io_files, line);
  io_files >> delta_r_MC;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_neighbor_list;
  getline(io_files, line);
  io_files >> nl_para;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> dt_minimizer;
  getline(io_files, line);
  io_files >> dt_reset_system;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> NUMBER_OF_QUENCHES;
  getline(io_files, line);
  io_files >> NUMBER_OF_EQB_STEPS;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> U_tol_0;
  getline(io_files, line);
  io_files >> F_tol_0;
  getline(io_files, line);
  io_files >> cv_tol;
  getline(io_files, line);
  io_files >> sigma_tol_0;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_eigenvalue_solver;
  getline(io_files, line);
  io_files >> print_eigenvector;
  getline(io_files, line);
  io_files >> number_of_eigenvalues;
  getline(io_files, line);
  io_files >> max_iters_eig;
  getline(io_files, line);
  io_files >> eigs_tol;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> type_force_flag;
  getline(io_files, line);
  io_files >> particle_fix_flag;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> random_seed_size_dist;
  getline(io_files, line);
  io_files >> random_seed_pos_init;
  getline(io_files, line);
  io_files >> random_seed_pos_init_discard;
  getline(io_files, line);
  io_files >> random_seed_random_number;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> init_filename;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> run_parameters_filename;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> output_mode;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> frame_print_frequency;
  getline(io_files, line);
  io_files >> utils_print_frequency;
  getline(io_files, line);
  io_files >> tau_print_frequency;
  getline(io_files, line);
  io_files >> hessian_print_frequency;
  getline(io_files, line);
  io_files >> structure_print_frequency;
  getline(io_files, line);
  io_files >> structure_factor_print_frequency;
  getline(io_files, line);
  io_files >> data_print_frequency;
  getline(io_files, line);
  io_files >> individual_data_print_frequency;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> ds_print_frequency;
  getline(io_files, line);
  io_files >> dcontour_print_frequency;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> avalanche_count;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> output_PATH;
  getline(io_files, line);
  io_files >> U_filename;
  getline(io_files, line);
  io_files >> Z_filename;
  getline(io_files, line);
  io_files >> R_filename;
  getline(io_files, line);
  io_files >> s_contour_filename;
  getline(io_files, line);
  io_files >> displacement_filename;
  getline(io_files, line);
  io_files >> order_parameter_filename;
  getline(io_files, line);
  io_files >> r_ij_parameter_filename;
  getline(io_files, line);
  io_files >> r_ij_filename;
  getline(io_files, line);
  io_files >> structure_filename;
  getline(io_files, line);
  io_files >> U_i_filename;
  getline(io_files, line);
  io_files >> Z_i_filename;
  getline(io_files, line);
  io_files >> R_i_filename;
  getline(io_files, line);
  io_files >> s_i_filename;
  getline(io_files, line);
  io_files >> contour_i_filename;
  getline(io_files, line);
  io_files >> ds_i_filename;
  getline(io_files, line);
  io_files >> dcontour_i_filename;
  getline(io_files, line);
  io_files >> displacement_i_filename;
  getline(io_files, line);
  io_files >> ddisplacement_i_filename;
  getline(io_files, line);
  io_files >> eigenvalues_filename;
  getline(io_files, line);
  io_files >> eigenvectors_filename;
  getline(io_files, line);
  io_files >> frame_filename;
  getline(io_files, line);
  io_files >> utils_filename;
  getline(io_files, line);
  io_files >> structure_factor_filename;
  getline(io_files, line);
  io_files >> tolerances_filename;
  getline(io_files, line);
  io_files >> unphysical_filename;
  getline(io_files, line);
  io_files >> displacement_vector_filename;
  getline(io_files, line);
  getline(io_files, line);

  io_files >> analysis_PATH;
  getline(io_files, line);
  io_files >> analysis_filename;
  getline(io_files, line);

  io_files >> unit_vector_flag;
  getline(io_files, line);
  io_files >> n_size;
  unit_vector_r.assign(n_size, -1);

  for (int k = 0; k < n_size; ++k)
  {
    io_files >> temp;
    unit_vector_r[k] = temp;
  }
  getline(io_files, line);

  io_files.close();
}

void output()
{
  io_files.open((output_PATH + run_parameters_filename).c_str(), ios::out | ios::trunc);

  analysis_files << setprecision(numeric_limits<double>::digits10) << scientific;

  io_files << N_init << " ";
  io_files << "\t\t\t\tN_init";
  io_files << '\n';
  io_files << '\n';

  for (int k = 0; k < d; ++k)
  {
    io_files << L_box[k] << " ";
  }
  io_files << "\t\t\t\tL_box";
  io_files << '\n';
  for (int k = 0; k < d; ++k)
  {
    io_files << dL_box[k] << " ";
  }
  io_files << "\t\t\t\tdL_box";
  io_files << '\n';
  io_files << '\n';

  io_files << vol_frac << " ";
  io_files << "\t\t\t\tvol_frac";
  io_files << '\n';
  io_files << dvol_frac_0 << " ";
  io_files << "\t\t\t\tdvol_frac";
  io_files << '\n';
  io_files << '\n';

  io_files << b << " ";
  io_files << "\t\t\t\tb";
  io_files << '\n';
  io_files << '\n';

  io_files << type_kbT << " ";
  io_files << "\t\t\t\ttype_kbT";
  io_files << '\n';
  io_files << kbTi << " ";
  io_files << "\t\t\t\tkbTi";
  io_files << '\n';
  io_files << kbTf << " ";
  io_files << "\t\t\t\tkbTf";
  io_files << '\n';
  io_files << '\n';

  io_files << type_force << " ";
  io_files << "\t\t\t\ttype_force";
  io_files << '\n';
  io_files << type_epsilon << " ";
  io_files << "\t\t\t\ttype_epsilon";
  io_files << '\n';
  io_files << type_sigma << " ";
  io_files << "\t\t\t\ttype_sigma";
  io_files << '\n';
  io_files << force_para_1 << " ";
  io_files << "\t\t\t\tforce_para_1";
  io_files << '\n';
  io_files << force_para_2 << " ";
  io_files << "\t\t\t\tforce_para_2";
  io_files << '\n';
  io_files << r_force_cut_off << " ";
  io_files << "\t\t\t\tr_force_cut_off";
  io_files << '\n';
  io_files << '\n';

  io_files << type_distribution << " ";
  io_files << "\t\t\t\ttype_distribution";
  io_files << '\n';
  io_files << para_1 << " ";
  io_files << "\t\t\t\tpara_1";
  io_files << '\n';
  io_files << para_2 << " ";
  io_files << "\t\t\t\tpara_2";
  io_files << '\n';
  io_files << para_3 << " ";
  io_files << "\t\t\t\tpara_3";
  io_files << '\n';
  io_files << '\n';

  io_files << type_mass << " ";
  io_files << "\t\t\t\ttype_mass";
  io_files << '\n';
  io_files << '\n';

  io_files << type_BC << " ";
  io_files << "\t\t\t\ttype_BC";
  io_files << '\n';
  io_files << '\n';

  for (int k = 0; k < d; ++k)
  {
    io_files << N_grid_sf[k] << " ";
  }
  io_files << "\t\t\t\tN_grid_sf";
  io_files << '\n';
  for (int k = 0; k < d; ++k)
  {
    io_files << N_grid_sf_max[k] << " ";
  }
  io_files << "\t\t\t\tN_grid_sf_max";
  io_files << '\n';
  io_files << '\n';

  io_files << type_initialization << " ";
  io_files << "\t\t\t\ttype_initialization";
  io_files << '\n';
  io_files << '\n';

  io_files << type_pre_start << " ";
  io_files << "\t\t\t\ttype_pre_start";
  io_files << '\n';
  io_files << '\n';

  io_files << type_quenching << " ";
  io_files << "\t\t\t\ttype_quenching";
  io_files << '\n';
  io_files << type_reset_system << " ";
  io_files << "\t\t\t\ttype_reset_system";
  io_files << '\n';
  io_files << '\n';

  io_files << type_relaxation_q << " ";
  io_files << "\t\t\t\ttype_relaxation_q";
  io_files << '\n';
  io_files << type_relaxation_r << " ";
  io_files << "\t\t\t\ttype_relaxation_r";
  io_files << '\n';
  io_files << '\n';

  io_files << type_controller << " ";
  io_files << "\t\t\t\ttype_controller";
  io_files << '\n';
  io_files << control_variable_final << " ";
  io_files << "\t\t\t\tcontrol_variable_final";
  io_files << '\n';
  io_files << '\n';

  for (int k = 0; k < 2; ++k)
  {
    io_files << bl_para[k] << " ";
  }
  io_files << "\t\t\t\tbl_para";
  io_files << '\n';
  io_files << '\n';

  io_files << bias_U_0 << " ";
  io_files << "\t\t\t\tbias_U_0";
  io_files << '\n';
  io_files << bias_U_sigma << " ";
  io_files << "\t\t\t\tbias_U_sigma";
  io_files << '\n';
  io_files << bias_Z_cut_off << " ";
  io_files << "\t\t\t\tbias_Z_cut_off";
  io_files << '\n';
  io_files << well_U_0 << " ";
  io_files << "\t\t\t\twell_U_0";
  io_files << '\n';
  io_files << well_U_sigma << " ";
  io_files << "\t\t\t\twell_U_sigma";
  io_files << '\n';
  io_files << well_cut_off << " ";
  io_files << "\t\t\t\twell_cut_off";
  io_files << '\n';
  io_files << '\n';

  io_files << type_dynamics_q << " ";
  io_files << "\t\t\t\ttype_dynamics_q";
  io_files << '\n';
  io_files << type_dynamics_r << " ";
  io_files << "\t\t\t\ttype_dynamics_r";
  io_files << '\n';
  io_files << '\n';

  io_files << type_MD_MC_q << " ";
  io_files << "\t\t\t\ttype_MD_MC_q";
  io_files << '\n';
  io_files << type_MD_MC_r << " ";
  io_files << "\t\t\t\ttype_MD_MC_r";
  io_files << '\n';
  io_files << '\n';

  io_files << p_swap << " ";
  io_files << "\t\t\t\tp_swap";
  io_files << '\n';
  io_files << delta_r_MC << " ";
  io_files << "\t\t\t\tdelta_r_MC";
  io_files << '\n';
  io_files << '\n';

  io_files << type_neighbor_list << " ";
  io_files << "\t\t\t\ttype_neighbor_list";
  io_files << '\n';
  io_files << nl_para << " ";
  io_files << "\t\t\t\tnl_para";
  io_files << '\n';
  io_files << '\n';

  io_files << dt_minimizer << " ";
  io_files << "\t\t\t\tdt_minimizer";
  io_files << '\n';
  io_files << dt_reset_system << " ";
  io_files << "\t\t\t\tdt_reset_system";
  io_files << '\n';
  io_files << '\n';

  io_files << NUMBER_OF_QUENCHES << " ";
  io_files << "\t\t\t\tNUMBER_OF_QUENCHES";
  io_files << '\n';
  io_files << NUMBER_OF_EQB_STEPS << " ";
  io_files << "\t\t\t\tNUMBER_OF_EQB_STEPS";
  io_files << '\n';
  io_files << '\n';

  io_files << U_tol_0 << " ";
  io_files << "\t\t\t\tU_tol_0";
  io_files << '\n';
  io_files << F_tol_0 << " ";
  io_files << "\t\t\t\tF_tol_0";
  io_files << '\n';
  io_files << cv_tol << " ";
  io_files << "\t\t\t\tcv_tol";
  io_files << '\n';
  io_files << sigma_tol_0 << " ";
  io_files << "\t\t\t\tsigma_tol_0";
  io_files << '\n';
  io_files << '\n';

  io_files << type_eigenvalue_solver << " ";
  io_files << "\t\t\t\ttype_eigenvalue_solver";
  io_files << '\n';
  io_files << print_eigenvector << " ";
  io_files << "\t\t\t\tprint_eigenvector";
  io_files << '\n';
  io_files << number_of_eigenvalues << " ";
  io_files << "\t\t\t\tnumber_of_eigenvalues";
  io_files << '\n';
  io_files << max_iters_eig << " ";
  io_files << "\t\t\t\tmax_iters_eig";
  io_files << '\n';
  io_files << eigs_tol << " ";
  io_files << "\t\t\t\teigs_tol";
  io_files << '\n';
  io_files << '\n';

  io_files << type_force_flag << " ";
  io_files << "\t\t\t\ttype_force_flag";
  io_files << '\n';
  io_files << particle_fix_flag << " ";
  io_files << "\t\t\t\tparticle_fix_flag";
  io_files << '\n';
  io_files << '\n';

  io_files << random_seed_size_dist << " ";
  io_files << "\t\t\t\trandom_seed_size_dist";
  io_files << '\n';
  io_files << random_seed_pos_init << " ";
  io_files << "\t\t\t\trandom_seed_pos_init";
  io_files << '\n';
  io_files << random_seed_pos_init_discard << " ";
  io_files << "\t\t\t\trandom_seed_pos_init_discard";
  io_files << '\n';
  io_files << random_seed_random_number << " ";
  io_files << "\t\t\t\trandom_seed_random_number";
  io_files << '\n';
  io_files << '\n';

  io_files << init_filename << " ";
  io_files << "\t\t\t\tinit_filename";
  io_files << '\n';
  io_files << '\n';

  io_files << run_parameters_filename << " ";
  io_files << "\t\t\t\trun_parameters_filename";
  io_files << '\n';
  io_files << '\n';

  io_files << output_mode << " ";
  io_files << "\t\t\t\toutput_mode";
  io_files << '\n';
  io_files << '\n';

  io_files << frame_print_frequency << " ";
  io_files << "\t\t\t\tframe_print_frequency";
  io_files << '\n';
  io_files << utils_print_frequency << " ";
  io_files << "\t\t\t\tutils_print_frequency";
  io_files << '\n';
  io_files << tau_print_frequency << " ";
  io_files << "\t\t\t\ttau_print_frequency";
  io_files << '\n';
  io_files << hessian_print_frequency << " ";
  io_files << "\t\t\t\thessian_print_frequency";
  io_files << '\n';
  io_files << structure_print_frequency << " ";
  io_files << "\t\t\t\tstructure_print_frequency";
  io_files << '\n';
  io_files << structure_factor_print_frequency << " ";
  io_files << "\t\t\t\tstructure_factor_print_frequency";
  io_files << '\n';
  io_files << data_print_frequency << " ";
  io_files << "\t\t\t\tdata_print_frequency";
  io_files << '\n';
  io_files << individual_data_print_frequency << " ";
  io_files << "\t\t\t\tindividual_data_print_frequency";
  io_files << '\n';
  io_files << '\n';

  io_files << ds_print_frequency << " ";
  io_files << "\t\t\t\tds_print_frequency";
  io_files << '\n';
  io_files << dcontour_print_frequency << " ";
  io_files << "\t\t\t\tdcontour_print_frequency";
  io_files << '\n';
  io_files << '\n';

  io_files << avalanche_count << " ";
  io_files << "\t\t\t\tavalanche_count";
  io_files << '\n';
  io_files << '\n';

  io_files << output_PATH << " ";
  io_files << "\t\t\t\toutput_PATH";
  io_files << '\n';
  io_files << U_filename << " ";
  io_files << "\t\t\t\tU_filename";
  io_files << '\n';
  io_files << Z_filename << " ";
  io_files << "\t\t\t\tZ_filename";
  io_files << '\n';
  io_files << R_filename << " ";
  io_files << "\t\t\t\tR_filename";
  io_files << '\n';
  io_files << s_contour_filename << " ";
  io_files << "\t\t\t\ts_contour_filename";
  io_files << '\n';
  io_files << displacement_filename << " ";
  io_files << "\t\t\t\tdisplacement_filename";
  io_files << '\n';
  io_files << order_parameter_filename << " ";
  io_files << "\t\t\t\torder_parameter_filename";
  io_files << '\n';
  io_files << r_ij_parameter_filename << " ";
  io_files << "\t\t\t\tr_ij_parameter_filename";
  io_files << '\n';
  io_files << r_ij_filename << " ";
  io_files << "\t\t\t\tr_ij_filename";
  io_files << '\n';
  io_files << structure_filename << " ";
  io_files << "\t\t\t\tstructure_filename";
  io_files << '\n';
  io_files << U_i_filename << " ";
  io_files << "\t\t\t\tU_i_filename";
  io_files << '\n';
  io_files << Z_i_filename << " ";
  io_files << "\t\t\t\tZ_i_filename";
  io_files << '\n';
  io_files << R_i_filename << " ";
  io_files << "\t\t\t\tR_i_filename";
  io_files << '\n';
  io_files << s_i_filename << " ";
  io_files << "\t\t\t\ts_i_filename";
  io_files << '\n';
  io_files << contour_i_filename << " ";
  io_files << "\t\t\t\tcontour_i_filename";
  io_files << '\n';
  io_files << ds_i_filename << " ";
  io_files << "\t\t\t\tds_i_filename";
  io_files << '\n';
  io_files << dcontour_i_filename << " ";
  io_files << "\t\t\t\tdcontour_i_filename";
  io_files << '\n';
  io_files << displacement_i_filename << " ";
  io_files << "\t\t\t\tdisplacement_i_filename";
  io_files << '\n';
  io_files << ddisplacement_i_filename << " ";
  io_files << "\t\t\t\tddisplacement_i_filename";
  io_files << '\n';
  io_files << eigenvalues_filename << " ";
  io_files << "\t\t\t\teigenvalues_filename";
  io_files << '\n';
  io_files << eigenvectors_filename << " ";
  io_files << "\t\t\t\teigenvectors_filename";
  io_files << '\n';
  io_files << frame_filename << " ";
  io_files << "\t\t\t\tframe_filename";
  io_files << '\n';
  io_files << utils_filename << " ";
  io_files << "\t\t\t\tutils_filename";
  io_files << '\n';
  io_files << structure_factor_filename << " ";
  io_files << "\t\t\t\tstructure_factor_filename";
  io_files << '\n';
  io_files << tolerances_filename << " ";
  io_files << "\t\t\t\ttolerances_filename";
  io_files << '\n';
  io_files << unphysical_filename << " ";
  io_files << "\t\t\t\tunphysical_filename";
  io_files << '\n';
  io_files << displacement_vector_filename << " ";
  io_files << "\t\t\t\tdisplacement_vector_filename";
  io_files << '\n';
  io_files << '\n';

  io_files << analysis_PATH << " ";
  io_files << "\t\t\t\tanalysis_PATH";
  io_files << '\n';
  io_files << analysis_filename << " ";
  io_files << "\t\t\t\tanalysis_filename";
  io_files << '\n';
  io_files << '\n';

  io_files << unit_vector_flag << " ";
  io_files << "\t\t\t\tunit_vector_flag";
  io_files << '\n';
  for (unsigned int k = 0; k < unit_vector_r.size(); ++k)
  {
    io_files << unit_vector_r[k] << " ";
  }
  io_files << "\t\t\t\tunit_vector_r";

  io_files.close();
}

void file_init()
{
  string line;
  double temp;

  VectorXd L_box_old = VectorXd::Zero(d);

  bubbles[0].volume = 0.0;

  io_files.open((input_PATH + init_filename).c_str(), ios::in);
  io_files.clear();
  io_files.seekg(0, ios::beg);

  getline(io_files, line);
  io_files >> time_system;

  getline(io_files, line);

  getline(io_files, line);
  io_files >> bubbles[0].N;

  getline(io_files, line);

  getline(io_files, line);
  for (int k = 0; k < d; ++k)
  {
    io_files >> temp;
    io_files >> L_box_old[k];
  }

  getline(io_files, line);
  getline(io_files, line);

  particle p_temp;

  int i = 1;
  while ((!io_files.eof()) && (i <= N_init))
  {
    io_files >> p_temp.real_index[0];
    io_files >> p_temp.type;
    io_files >> p_temp.particle_index;

    while (i < p_temp.particle_index)
    {
      bubbles[i].flag = 0;
      ++i;
    }
    bubbles[i].flag = 1;
    bubbles[i].mass = 1.0;
    bubbles[i].real_index[0] = p_temp.real_index[0];
    bubbles[i].type = p_temp.type;
    bubbles[i].particle_index = p_temp.particle_index;

    for (int k = 0; k < d; ++k)
    {
      io_files >> bubbles[i].position[k];
    }

    io_files >> bubbles[i].radius;

    bubbles[i].volume = (2.0 * (d - 1) / d) * pi;
    if (d == 2)
    {
      bubbles[i].volume *= gsl_pow_2(bubbles[i].radius);
    }
    else if (d == 3)
    {
      bubbles[i].volume *= gsl_pow_3(bubbles[i].radius);
    }
    else
    {
      bubbles[i].volume *= pow(bubbles[i].radius, d);
    }
    bubbles[0].volume += bubbles[i].volume;

    ++i;
  }
  while (i <= N_init)
  {
    bubbles[i].flag = 0;
    ++i;
  }

  //type_initialization = -1;

  if ((L_box[0] == -1) && (vol_frac != -1))
  {
    if (d == 2)
    {
      L_box.setConstant(d, sqrt(bubbles[0].volume / vol_frac));
    }
    else if (d == 3)
    {
      L_box.setConstant(d, cbrt(bubbles[0].volume / vol_frac));
    }
    else
    {
      L_box.setConstant(d, pow(bubbles[0].volume / vol_frac, 1.0 / d));
    }
  }
  else if ((vol_frac == -1) && (L_box[0] != -1))
  {
    if (d == 2)
    {
      vol_frac = bubbles[0].volume / (L_box[0] * L_box[1]);
    }
    else if (d == 3)
    {
      vol_frac = bubbles[0].volume / (L_box[0] * L_box[1] * L_box[2]);
    }
    else
    {
      L_box.setConstant(d, pow(bubbles[0].volume / vol_frac, 1.0 / d));
    }
  }
  else if ((vol_frac == -1) && (L_box[0] == -1))
  {
    L_box = L_box_old;

    if (d == 2)
    {
      vol_frac = bubbles[0].volume / (L_box[0] * L_box[1]);
    }
    else if (d == 3)
    {
      vol_frac = bubbles[0].volume / (L_box[0] * L_box[1] * L_box[2]);
    }
    else
    {
      L_box.setConstant(d, pow(bubbles[0].volume / vol_frac, 1.0 / d));
    }
  }
  else
  {
    cout << "No box length and volume fraction available" << '\n' << "program run stopped \n";

    time_final = clock();
    run_time = double(time_final - time_initial) / CLOCKS_PER_SEC;

    cout << "time taken :" << run_time << " seconds";

    exit(0);
  }

  for (int i = 0; i <= N_init; ++i)
  {
    bubbles[i].position = bubbles[i].position.array() * L_box.array() / L_box_old.array();
  }

  io_files.close();
}

void create_lammpstrj(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, int option)
{
  filename += convert_int(system_counter) + "_" + convert_double(system_time) + "_" + convert_int(phyproc_counter) + "_" + convert_double(phyproc_time) + "_" + convert_double(bubbles[0].state_flag) + ".lammpstrj";

  io_files.open((output_PATH + filename).c_str(), ios::out | ios::trunc);
  io_files << setprecision(numeric_limits<double>::digits10) << fixed << "ITEM: TIMESTEP" << '\n';
  io_files << system_time << '\n';
  io_files << "ITEM: NUMBER OF ATOMS" << '\n';
  io_files << bubbles[0].N << '\n';
  io_files << "ITEM: BOX BOUNDS pp pp pp" << '\n';
  io_files << 0.0 << " " << L_box[0] << '\n';
  io_files << 0.0 << " " << L_box[1] << '\n';

  if (d == 3)
  {
    io_files << 0.0 << " " << L_box[2] << '\n';
  }

  io_files << "ITEM: ATOMS id type mol x y ";
  if (d == 3)
  {
    io_files << "z ";
  }
  io_files << "radius" << '\n';

  switch (option)
  {
  case 0:
    for (int i = 1; i <= N_init; ++i)
    {
      if (bubbles[i].flag == 1)
      {
        if (bubbles[i].real_index[0] > 1)
        {
          io_files << '\n';
        }

        io_files << bubbles[i].real_index[0] << " " << bubbles[i].type << " " << bubbles[i].particle_index;
        for (int k = 0; k < d; ++k)
        {
          io_files << " " << bubbles[i].position[k];
        }

        io_files << " " << bubbles[i].radius;
      }
    }
    break;

  case 1:
    for (int i = 1; i <= N_init; ++i)
    {
      if (bubbles_thermal[i].flag == 1)
      {
        if (bubbles_thermal[i].real_index[0] > 1)
        {
          io_files << '\n';
        }

        io_files << bubbles_thermal[i].real_index[0] << " " << 1 << " " << bubbles_thermal[i].particle_index;
        for (int k = 0; k < d; ++k)
        {
          io_files << " " << bubbles_thermal[i].position[k];
        }

        io_files << " " << bubbles_thermal[i].radius;
      }
    }

    break;
  }

  io_files.close();
}

void output_utils(const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  string filename = utils_filename + convert_int(system_counter) + "_" + convert_double(system_time) + "_" + convert_int(phyproc_counter) + "_" + convert_double(phyproc_time) + "_" + convert_double(bubbles[0].state_flag) + ".dat";

  io_files.open((output_PATH + filename).c_str(), ios::out | ios::binary | ios::trunc);
  io_files.clear();

  for (int i = 0; i <= N_init; ++i)
  {
    io_files.write((char *)&system_counter, sizeof(int));

    io_files.write((char *)&system_time, sizeof(double));

    io_files.write((char *)&phyproc_counter, sizeof(int));

    io_files.write((char *)&phyproc_time, sizeof(double));

    io_files.write((char *)&bubbles[0].state_flag, sizeof(double));

    io_files.write((char *)&bubbles[0].N, sizeof(int));

    io_files.write((char *)&vol_frac, sizeof(double));

    for (int k = 0; k < d; ++k)
    {
      io_files.write((char *)&L_box[k], sizeof(double));
    }
    io_files.write((char *)&bubbles[i].particle_index, sizeof(int));

    io_files.write((char *)&bubbles[i].real_index[0], sizeof(int));

    io_files.write((char *)&bubbles[i].flag, sizeof(int));

    io_files.write((char *)&bubbles[i].radius, sizeof(double));

    for (int k = 0; k < d; ++k)
    {
      io_files.write((char *)&bubbles[i].position[k], sizeof(double));
    }

    for (int k = 0; k < d; ++k)
    {
      io_files.write((char *)&bubbles[i].displacement[k], sizeof(double));
    }

    for (int k = 0; k < d; ++k)
    {
      io_files.write((char *)&bubbles[i].ddisplacement[k], sizeof(double));
    }

    io_files.write((char *)&bubbles[i].s, sizeof(double));

    io_files.write((char *)&bubbles[i].ds, sizeof(double));

    io_files.write((char *)&bubbles[i].contour, sizeof(double));

    io_files.write((char *)&bubbles[i].dcontour, sizeof(double));

    io_files.write((char *)&bubbles[i].U, sizeof(double));

    io_files.write((char *)&bubbles[i].dU, sizeof(double));

    io_files.write((char *)&bubbles[i].dU_U, sizeof(double));

    io_files.write((char *)&bubbles[i].Z[0], sizeof(double));

    for (int l = 0; l < d; ++l)
    {
      for (int m = 0; m < d; ++m)
      {
        io_files.write((char *)&bubbles[i].Tau(l, m), sizeof(double));
      }
    }
  }

  //io_files.write((char *)&bubbles, sizeof(particle)*(N_init + 1));

  io_files.close();
}

void input_utils()
{
  bubbles[0].volume = 0.0;

  io_files.open((input_PATH + utils_filename).c_str(), ios::in | ios::binary);
  io_files.clear();
  io_files.seekg(0, ios::beg);

  bubbles[0].N = 0;

  double temp;

  io_files.write((char *)&temp, sizeof(int));

  io_files.write((char *)&time_system, sizeof(double));

  io_files.write((char *)&temp, sizeof(int));

  io_files.write((char *)&temp, sizeof(double));

  io_files.read((char *)&bubbles[0].state_flag, sizeof(double));

  io_files.read((char *)&N_init, sizeof(int));

  io_files.read((char *)&vol_frac, sizeof(double));

  for (int k = 0; k < d; ++k)
  {
    io_files.read((char *)&L_box[k], sizeof(double));
  }

  for (int i = 0; i <= N_init; ++i)
  {
    io_files.read((char *)&bubbles[i].particle_index, sizeof(int));

    io_files.read((char *)&bubbles[i].real_index[0], sizeof(int));

    io_files.read((char *)&bubbles[i].flag, sizeof(int));
    if (bubbles[i].flag == 1)
    {
      ++bubbles[0].N;
    }

    io_files.read((char *)&bubbles[i].radius, sizeof(double));

    bubbles[i].volume = (2.0 * (d - 1) / d) * pi;
    if (d == 2)
    {
      bubbles[i].volume *= gsl_pow_2(bubbles[i].radius);
    }
    else if (d == 3)
    {
      bubbles[i].volume *= gsl_pow_3(bubbles[i].radius);
    }
    else
    {
      bubbles[i].volume *= pow(bubbles[i].radius, d);
    }
    bubbles[0].volume += bubbles[i].volume;

    for (int k = 0; k < d; ++k)
    {
      io_files.read((char *)&bubbles[i].position[k], sizeof(double));
    }

    for (int k = 0; k < d; ++k)
    {
      io_files.read((char *)&bubbles[i].displacement[k], sizeof(double));
    }

    for (int k = 0; k < d; ++k)
    {
      io_files.read((char *)&bubbles[i].ddisplacement[k], sizeof(double));
    }

    io_files.read((char *)&bubbles[i].s, sizeof(double));

    io_files.read((char *)&bubbles[i].ds, sizeof(double));

    io_files.read((char *)&bubbles[i].contour, sizeof(double));

    io_files.read((char *)&bubbles[i].dcontour, sizeof(double));

    io_files.read((char *)&bubbles[i].U, sizeof(double));

    io_files.read((char *)&bubbles[i].dU, sizeof(double));

    io_files.read((char *)&bubbles[i].dU_U, sizeof(double));

    io_files.read((char *)&bubbles[i].Z[0], sizeof(int));

    for (int l = 0; l < d - 1; ++l)
    {
      for (int m = l + 1; m < d; ++m)
      {
        io_files.read((char *)&bubbles[i].Tau(l, m), sizeof(double));
      }
    }
  }
  io_files.close();
}

void output_analysis(int output_flag, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, int &avalanchecounter, int force_option)
{
  string filename;
  if ((output_flag <= 0) || ((output_flag <= output_mode) && (system_time > 0)))
  {
    if ((output_flag <= 0) && (system_time <= 0))
    {
      bubbles[0].state_flag = output_flag;
    }
    int step_counter;
    if ((bubbles[0].state_flag == 0) || (bubbles[0].state_flag == 4))
    {
      step_counter = system_counter;
    }
    else
    {
      step_counter = phyproc_counter;
    }
    if (data_print_frequency != 0)
    {
      if (((step_counter % data_print_frequency == 0) /*  && (step_counter != 0) */) || (output_flag <= 0))
      {
        output_tolerances(tolerances_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_U(U_filename, system_counter, system_time, phyproc_counter, phyproc_time, min_counter_SUCCESSFUL, min_counter);
        output_Z(Z_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_R(R_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_s_contour(s_contour_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_displacement(displacement_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        // output_order_parameter(order_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        output_r_ij_parameter(r_ij_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_r_ij(r_ij_filename, system_counter, system_time, phyproc_counter, phyproc_time, 3);
      }
    }

    if (data_print_frequency != 0)
    {
      while ((((output_flag == 1) && (bubbles[0].s >= ds_print_counter * ds_print_frequency)) || ((output_flag >= 2) && (bubbles[0].dcontour >= ds_print_counter * ds_print_frequency))) && (ds_print_frequency != 0))
      {
        output_tolerances("ds_" + tolerances_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_U("ds_" + U_filename, system_counter, system_time, phyproc_counter, phyproc_time, min_counter_SUCCESSFUL, min_counter);
        output_Z("ds_" + Z_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_R("ds_" + R_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_s_contour("ds_" + s_contour_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_s_contour("ds_" + displacement_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        // output_order_parameter("ds_" + order_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        output_r_ij_parameter("ds_" + r_ij_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_r_ij("ds_" + r_ij_filename, system_counter, system_time, phyproc_counter, phyproc_time, 3);

        if (individual_data_print_frequency != 0)
        {
          output_U_i("ds_" + U_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_Z_i("ds_" + Z_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_R_i("ds_" + R_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_s_i("ds_" + s_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);

          output_displacement_i("ds_" + displacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_ddisplacement_i("ds_" + ddisplacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        }

        ds_print_counter++;
      }

      while ((((output_flag == 1) && (bubbles[0].contour >= dcontour_print_counter * dcontour_print_frequency)) || ((output_flag >= 2) && (bubbles[0].dcontour >= dcontour_print_counter * dcontour_print_frequency))) && (dcontour_print_frequency != 0))
      {
        output_tolerances("dcontour_" + tolerances_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_U("dcontour_" + U_filename, system_counter, system_time, phyproc_counter, phyproc_time, min_counter_SUCCESSFUL, min_counter);
        output_Z("dcontour_" + Z_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_R("dcontour_" + R_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_s_contour("dcontour_" + s_contour_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_s_contour("dcontour_" + displacement_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        // output_order_parameter("dcontour_" + order_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);

        output_r_ij_parameter("dcontour_" + r_ij_parameter_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        output_r_ij("dcontour_" + r_ij_filename, system_counter, system_time, phyproc_counter, phyproc_time, 3);

        if (individual_data_print_frequency != 0)
        {
          output_U_i("dcontour_" + U_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_Z_i("dcontour_" + Z_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_R_i("dcontour_" + R_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_s_i("dcontour_" + s_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);

          output_displacement_i("dcontour_" + displacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_ddisplacement_i("dcontour_" + ddisplacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        }

        dcontour_print_counter++;
      }
    }

    if ((frame_print_frequency != 0) && (bubbles[0].state_flag == 0))
    {
      if (((step_counter % frame_print_frequency == 0) && (step_counter != 0)) || (output_flag <= 0))
      {
        int option = 0;
        // if (bubbles[0].state_flag == 1)
        // {
        //   option = 1;
        // }
        create_lammpstrj(frame_filename, system_counter, system_time, phyproc_counter, phyproc_time, option);
      }
    }
    else if (output_flag <= 0)
    {
      create_lammpstrj(frame_filename, system_counter, system_time, phyproc_counter, phyproc_time, 0);
    }

    if ((output_flag == 1) || ((output_mode == 5) && (output_flag >= 0)))
    {
      if (individual_data_print_frequency != 0)
      {
        if ((step_counter % individual_data_print_frequency == 0) /*  && (step_counter != 0) */ && (bubbles[0].state_flag == 0))
        {
          output_U_i(U_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_Z_i(Z_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_R_i(R_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_s_i(s_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);

          output_displacement_i(displacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
          output_ddisplacement_i(ddisplacement_i_filename, system_counter, system_time, phyproc_counter, phyproc_time);
        }
      }

      if (utils_print_frequency != 0)
      {
        if ((step_counter % utils_print_frequency == 0) /*  && (step_counter != 0) */  && (bubbles[0].state_flag == 0))
        {
          // if (tau_print_frequency != 0)
          // {
          //   if ((step_counter % tau_print_frequency == 0) /*  && (step_counter != 0) */)
          //   {
          //     update_Tau(force_option);
          //   }
          // }
          output_utils(system_counter, system_time, phyproc_counter, phyproc_time);
        }
      }

    //   if (hessian_print_frequency != 0)
    //   {
    //     if ((step_counter % hessian_print_frequency == 0) /*  && (step_counter != 0) */)
    //     {
    //       update_Hessian(force_option);
    //       update_Hessian_eigenvalues(number_of_eigenvalues, type_eigenvalue_solver);
    //       update_min_position();
    //       output_eigenvalue(eigenvalues_filename, system_counter, system_time, phyproc_counter, phyproc_time);
    //       if (print_eigenvector == 1)
    //       {
    //         output_eigenvector(eigenvectors_filename, system_counter, system_time, phyproc_counter, phyproc_time);
    //       }
    //     }
    //   }

    //   if (structure_print_frequency != 0)
    //   {
    //     if ((step_counter % structure_print_frequency == 0) /*  && (step_counter != 0) */)
    //     {
    //       output_structure(structure_filename, system_counter, system_time, phyproc_counter, phyproc_time);
    //     }
    //   }

    //   if (structure_factor_print_frequency != 0)
    //   {
    //     if ((step_counter % structure_factor_print_frequency == 0) /*  && (step_counter != 0) */)
    //     {
    //       output_structure_factor(structure_factor_filename, system_counter, system_time, phyproc_counter, phyproc_time);
    //     }
    //   }
    }

    if (avalanchecounter != -1)
    {
      if ((avalanchecounter < avalanche_count) && (bubbles[0].dU_U >= 10 * bubbles_athermal[0].dU_U) && (bubbles[0].dU_U > 0) && (bubbles_athermal[0].dU_U > 0))
      {
        output_avalanche(++avalanchecounter, system_counter, system_time, phyproc_counter, phyproc_time);
      }
    }
  }
}