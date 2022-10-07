/*
 * File: classes.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 5:05:12 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 9th August 2021 6:23:22 pm
 * Modified By: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Copyright (c) 2020-2022 Amru, University of Pennsylvania
 * 
 * Summary: Fill In
 */
#include <headers.h>

#include <globals.h>
#include <classes.h>

particle::particle()
{
    state_flag = -1;

    N = -1;

    n_flag = -1;

    type = -1;

    radius = 0.0;

    volume = 0.0;

    com_flag = -1;

    mass = 0.0;

    sigma = -1;

    epsilon = 0;

    position.setZero(d);
    //real_position.setZero(d);

    old_position.setZero(d);
    //real_old_position.setZero(d);

    min_position.setZero(d);

    deltadisplacement.setZero(d);
    displacement.setZero(d);
    s_displacement = 0.0;
    s = 0.0;
    contour = 0.0;

    s_non_rattlers = 0.0;

    velocity.setZero(d);

    force.setZero(d);

    conj_force.setZero(d);

    f_flag = -1;

    U = 0.0;

    dU = 0.0;
    dU_U = 0.0;

    bias_force.setZero(d);

    bias_U = 0.0;

    U_non_bias_min = 0.0;
    U_non_bias_max = 0.0;

    F_max = 0.0;

    bias_U_position_list.clear();

    bias_skin_distance.setZero(d);
    bias_skin_length = 0.0;

    bias_U_position_com_neighbor_list.clear();
    bias_U_position_neighbor_list.clear();

    bl_flag = -1;

    KE = 0.0;

    dKE = 0.0;
    dKE_KE = 0.0;

    E = 0.0;

    dE = 0.0;
    dE_E = 0.0;

    eta.setZero(d);

    ddeltadisplacement.setZero(d);
    ddisplacement.setZero(d);
    ds = 0.0;
    dcontour = 0.0;

    ds_non_rattlers = 0.0;

    Tau.setZero(d, d);

    t_flag = -1;

    Q_flux = 0.0;

    skin_distance.setZero(d);
    skin_length = 0.0;
    skin_radius_length = 0.0;

    neighbor_list.clear();
    non_red_neighbor_list.clear();

    com_neighbor_list.clear();
    non_red_com_neighbor_list.clear();

    nl_flag = -1;

    N_non_rattlers.assign(1, -1);

    mass_non_rattlers.assign(1, -1);

    Z.assign(1, 0);

    p = 0.0;

    z_flag = -1;

    particle_index = -1;
    real_index.assign(1, -1);
    grid_index.setConstant(d, -1);

    flag = -1;
}

particle::~particle()
{
}