/*
 * File: classes.h
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 5:02:55 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 9th August 2021 6:23:19 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */

#ifndef CLASSES_H
#define CLASSES_H

#include <headers.h>

class particle
{
public:
    double state_flag;

    int N;

    int n_flag;

    int type;

    double radius;

    double volume;

    int com_flag;

    double mass;

    double sigma;

    double epsilon;

    VectorXd position;
    //VectorXd real_position;

    VectorXd old_position;
    //VectorXd real_old_position;

    VectorXd min_position;

    VectorXd deltadisplacement;
    VectorXd displacement;
    double s_displacement;
    double s;
    double contour;

    double s_non_rattlers;

    VectorXd velocity;

    VectorXd force;

    VectorXd conj_force;

    int f_flag;

    double U;

    double dU;
    double dU_U;

    VectorXd bias_force;

    double bias_U;

    double U_non_bias_min;
    double U_non_bias_max;

    double F_max;

    vector<VectorXd> bias_U_position_list;

    VectorXd bias_skin_distance;
    double bias_skin_length;

    vector<int> bias_U_position_com_neighbor_list;
    vector<int> bias_U_position_neighbor_list;

    int bl_flag;

    double KE;
    double dKE;
    double dKE_KE;

    double E;
    double dE;
    double dE_E;

    VectorXd eta;

    VectorXd ddeltadisplacement;
    VectorXd ddisplacement;
    double ds;
    double dcontour;

    double ds_non_rattlers;

    MatrixXd Tau;

    int t_flag;

    double Q_flux;

    VectorXd skin_distance;
    double skin_length;
    double skin_radius_length;

    vector<int> neighbor_list;
    vector<int> non_red_neighbor_list;

    vector<int> com_neighbor_list;
    vector<int> non_red_com_neighbor_list;

    int nl_flag;

    vector<int> N_non_rattlers;

    vector<double> mass_non_rattlers;

    vector<double> Z;

    double p;

    int z_flag;

    int particle_index;
    vector<int> real_index;
    VectorXd grid_index;

    int flag;

    particle();
    ~particle();
};

struct radii
{
    int particle_index;
    double radius;

    bool operator<(const radii &R) const
    {
        return (radius < R.radius);
    }

    bool operator>(const radii &R) const
    {
        return (radius > R.radius);
    }
};

#endif