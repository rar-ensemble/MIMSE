/*
 * File: analysis.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:16:20 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 10:41:48 pm
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

void output_U(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, const double &min_count_SUCCESSFUL, const int &min_count)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);
  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag U dU dU_U N min_counter_SUCCESSFUL min_counter U_min U_max F_max";

    if (type_kbT > 1)
    {
      analysis_files << " kbT";
    }

    if (type_quenching >= 3)
    {
      analysis_files << " bias_U number_of_biases";
    }
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << bubbles[0].U - (bias_counter[0] > 0) * bubbles[0].bias_U << " " << bubbles[0].dU << " " << bubbles[0].dU_U << " " << bubbles[0].N << " " << min_count_SUCCESSFUL << " " << min_count << " " << bubbles[0].U_non_bias_min << " " << bubbles[0].U_non_bias_max << " " << bubbles[0].F_max;

  if (type_kbT > 1)
  {
    analysis_files << " " << kbT;
  }

  if (type_quenching >= 3)
  {
    analysis_files << " " << bubbles[0].bias_U / bias_U_0 << " " << bias_counter[0];
  }

  analysis_files.close();
}

void output_Z(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  update_Z();

  update_p();

  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag Nc >= 0 Nc >= 1 Nc >= 2 Nc >= 3 Nc >= 4 Nc >= 5 Nc >= 6 N_non_rattlers[0] [1] [2] [3] [4] [5] [6] vol_frac p N";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int j = 0; j < (2 * d + 1); ++j)
  {
    analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[0].Z[j] << " ";
  }

  for (int j = 0; j < (2 * d + 1); ++j)
  {
    analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[0].N_non_rattlers[j] << " ";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << vol_frac << " " << bubbles[0].p << " " << bubbles[0].N;

  analysis_files.close();
}

void output_R(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  vector<int> indices = largest_radius(3);

  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag N (mass_init-mass_N) mean_radius() pow(mean_radius(),2) stddev_radius() R_1 R_2 R_3 R12 < L";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << bubbles[0].N << " " << (mass_init - mass_N) << " " << bubbles[0].radius << " " << gsl_pow_2(bubbles[0].radius) << " " << stddev_radius();

  for (unsigned int k = 0; k < indices.size(); ++k)
  {
    analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << " " << bubbles[indices[k]].radius;
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << " " << (2.0 * (bubbles[indices[0]].radius + bubbles[indices[1]].radius) < L_box.minCoeff());

  analysis_files.close();
}

void output_s_contour(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag bubbles[0].s_displacement bubbles[0].s bubbles[0].ds bubbles[0].contour bubbles[0].dcontour N bubbles[0].s_non_rattlers bubbles[0].ds_non_rattlers N_non_rattlers";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << bubbles[0].s_displacement << " " << bubbles[0].s << " " << bubbles[0].ds << " " << bubbles[0].contour << " " << bubbles[0].dcontour << " " << bubbles[0].N << " " << bubbles[0].s_non_rattlers << " " << bubbles[0].ds_non_rattlers << " " << bubbles[0].N_non_rattlers[2 * d - 2];

  analysis_files.close();
}

void output_displacement(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag bubbles[0].displacement_norm bubbles[0].ddisplacement_norm N";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << bubbles[0].displacement.norm() << " " << bubbles[0].ddisplacement.norm() << " " << bubbles[0].N;

  analysis_files.close();
}

void output_order_parameter(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  // update_apollonian_order_parameter();

  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag apollonian_order_parameter_system N";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << apollonian_order_parameter_system << " " << bubbles[0].N;

  analysis_files.close();
}

void output_r_ij_parameter(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  // update_r_ij();

  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag r_ij_parameter_system x_ij_system N";
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << r_ij_parameter_system << " " << x_ij_system << " " << bubbles[0].N;

  analysis_files.close();
}

void output_r_ij(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, int n)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int l = 1; l < n; ++l)
    {
      for (int m = l + 1; m < n + 1; ++m)
      {
        analysis_files << " r_" << l << m;
      }
    }
  }

  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag;

  vector<int> n_indices(n);

  n_indices = largest_radius(n);

  for (int l = 0; l < n - 1; ++l)
  {
    for (int m = l + 1; m < n; ++m)
    {
      int i = n_indices[l];
      int j = n_indices[m];

      VectorXd r_ij = periodic_BC(bubbles[i].position, bubbles[j].position, L_box);

      analysis_files << " " << r_ij.norm() / (0.5 * L_box.norm());
      // analysis_files << " " << r_ij_system.block((n_indices[i] - 1) * d, n_indices[j] - 1, d, 1).norm();
    }
  }

  analysis_files.close();
}

void output_val(double value, string filename)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  analysis_files << setprecision(numeric_limits<double>::digits10) << " " << value;

  analysis_files.close();
}

void output_avalanche(int avalanche_counter, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  string filename = "avalanche_";
  create_lammpstrj(filename, system_counter, system_time, avalanche_counter, -1, 1);
}

void output_U_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);
  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }
  analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].U - (bias_counter[0] > 0) * bubbles[i].bias_U << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_Z_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  update_Z();

  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << bubbles[i].Z[0] << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_R_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].radius << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_s_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].s << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_ds_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].ds << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_contour_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].contour << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_dcontour_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].dcontour << " ";
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_displacement_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      for (int k = 0; k < d; ++k)
      {
        analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].displacement[k] << " ";
      }
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_ddisplacement_i(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);

  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag";
    for (int i = 1; i <= N_init; ++i)
    {
      analysis_files << " " << i;
    }
  }

  analysis_files << setprecision(4) << fixed << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " ";

  for (int i = 1; i <= N_init; ++i)
  {
    if (bubbles[i].flag == 1)
    {
      for (int k = 0; k < d; ++k)
      {
        analysis_files << setprecision(numeric_limits<double>::digits10) << fixed << bubbles[i].ddisplacement[k] << " ";
      }
    }
    else
    {
      analysis_files << "nan ";
    }
  }

  analysis_files.close();
}

void output_tolerances(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time)
{
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);
  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag U_tol eff_U_0 F_tol eff_F_0";
  }
  analysis_files << setprecision(numeric_limits<double>::digits10) << scientific << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << U_tol << " " << eff_U_0 << " " << F_tol << " " << eff_F_0;

  analysis_files.close();
}

void output_displacement_vector(string filename, const int &system_counter, const double &system_time, const int &phyproc_counter, const double &phyproc_time, const unsigned int &j, const double &displacement_r)
{
  double ds_unit_vector = 0;
  analysis_files.open((output_PATH + filename).c_str(), ios::out | ios::app);
  if (analysis_files.tellp() == 0)
  {
    analysis_files << "system_counter system_time phyproc_counter phyproc_time status_flag j r";

    for (int i = 1; i <= N_init; ++i)
    {
      if ((bubbles[i].flag == 1) && (bubbles[i].bias_U_position_list[j][0] != -1))
      {
        for (int k = 0; k < d; ++k)
        {
          analysis_files << " " << i << "_" << k + 1;
        }
      }
    }

    analysis_files << " displacement_vector_r N";
  }
  analysis_files << setprecision(numeric_limits<double>::digits10) << scientific << '\n'
                 << system_counter << " " << system_time << " " << phyproc_counter << " " << phyproc_time << " " << bubbles[0].state_flag << " " << j << " " << displacement_r;

  for (int i = 1; i <= N_init; ++i)
  {
    if ((bubbles[i].flag == 1) && (bubbles[i].bias_U_position_list[j][0] != -1))
    {
      for (int k = 0; k < d; ++k)
      {
        analysis_files << " " << bubbles[i].ddisplacement[k];
      }

      ds_unit_vector += bubbles[i].ddisplacement.squaredNorm();
    }
  }

  ds_unit_vector = sqrt(ds_unit_vector);

  analysis_files << " " << ds_unit_vector << " " << bubbles[0].bias_U_position_list[j][1];

  analysis_files.close();
}