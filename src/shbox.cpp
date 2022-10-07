/*
 * File: shbox.cpp
 * Project: Quenching
 * File Created: Friday, 19th April 2020 7:02:37 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Tuesday, 26th April 2022 11:30:01 pm
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

void shrinking_box_relaxation(const int &relaxation_option, const int &dynamics_option, const int &MD_MC_option, const double &dt, const int &option, const double &control_variable_target)
{
    VectorXd L_box_target;

    int count = 0;

    counter_phyproc = count;
    time_phyproc = 0;

    double control_variable = 0.0;

    VectorXd L_box_old = L_box;

    double dvol_frac = dvol_frac_0;

    int comp_decomp_flag = 0;

    switch (option)
    {
    case 1:
        cv_tol = max(cv_tol, (dL_box.array() / L_box.array()).sum() * control_variable_target);
        if (d == 2)
        {
            L_box_target.setConstant(d, sqrt(bubbles[0].volume / control_variable_target));
        }
        else if (d == 3)
        {
            L_box_target.setConstant(d, cbrt(bubbles[0].volume / control_variable_target));
        }
        else
        {
            L_box_target.setConstant(d, pow(bubbles[0].volume / control_variable_target, 1.0 / d));
        }

        if (vol_frac == -1)
        {
            if (d == 2)
            {
                L_box.setConstant(d, sqrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else if (d == 3)
            {
                L_box.setConstant(d, cbrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else
            {
                L_box.setConstant(d, pow(bubbles[0].volume / (control_variable_target / 2.0), 1.0 / d));
            }
            update_scale_position(L_box_old, L_box);

            update_vol_frac();
        }
        control_variable = vol_frac;
        break;

    case 2:
        dL_box = dL_box.array().min(dvol_frac / vol_frac * L_box.array() / 3.0);
        cv_tol = max(cv_tol, (dL_box.array() / L_box.array()).sum() * control_variable_target);
        if (d == 2)
        {
            L_box_target.setConstant(d, sqrt(bubbles[0].volume / control_variable_target));
        }
        else if (d == 3)
        {
            L_box_target.setConstant(d, cbrt(bubbles[0].volume / control_variable_target));
        }
        else
        {
            L_box_target.setConstant(d, pow(bubbles[0].volume / control_variable_target, 1.0 / d));
        }

        if (vol_frac == -1)
        {
            if (d == 2)
            {
                L_box.setConstant(d, sqrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else if (d == 3)
            {
                L_box.setConstant(d, cbrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else
            {
                L_box.setConstant(d, pow(bubbles[0].volume / (control_variable_target / 2.0), 1.0 / d));
            }
            update_scale_position(L_box_old, L_box);

            update_vol_frac();
        }
        control_variable = vol_frac;
        break;

    case 3:
        dL_box = dL_box.array().min(dvol_frac / vol_frac * L_box.array() / 3.0);
        cv_tol = max(cv_tol, (dL_box.array() / L_box.array()).sum() * control_variable_target);
        if (d == 2)
        {
            L_box_target.setConstant(d, sqrt(bubbles[0].volume / control_variable_target));
        }
        else if (d == 3)
        {
            L_box_target.setConstant(d, cbrt(bubbles[0].volume / control_variable_target));
        }
        else
        {
            L_box_target.setConstant(d, pow(bubbles[0].volume / control_variable_target, 1.0 / d));
        }

        if (vol_frac == -1)
        {
            if (d == 2)
            {
                L_box.setConstant(d, sqrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else if (d == 3)
            {
                L_box.setConstant(d, cbrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else
            {
                L_box.setConstant(d, pow(bubbles[0].volume / (control_variable_target / 2.0), 1.0 / d));
            }
            update_scale_position(L_box_old, L_box);

            update_vol_frac();
        }
        control_variable = vol_frac;
        break;

    case 4:
        dL_box = dL_box.array().min(dvol_frac / vol_frac * L_box.array() / 3.0);
        cv_tol = max(cv_tol, (dL_box.array() / L_box.array()).sum() * control_variable_target);
        if (d == 2)
        {
            L_box_target.setConstant(d, sqrt(bubbles[0].volume / control_variable_target));
        }
        else if (d == 3)
        {
            L_box_target.setConstant(d, cbrt(bubbles[0].volume / control_variable_target));
        }
        else
        {
            L_box_target.setConstant(d, pow(bubbles[0].volume / control_variable_target, 1.0 / d));
        }

        if (vol_frac == -1)
        {
            if (d == 2)
            {
                L_box.setConstant(d, sqrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else if (d == 3)
            {
                L_box.setConstant(d, cbrt(bubbles[0].volume / (control_variable_target / 2.0)));
            }
            else
            {
                L_box.setConstant(d, pow(bubbles[0].volume / (control_variable_target / 2.0), 1.0 / d));
            }
            update_scale_position(L_box_old, L_box);

            update_vol_frac();
        }
        control_variable = vol_frac;
        break;

    default:
        cv_tol = max(cv_tol, (L_box.array() / dL_box.array()).sum() * control_variable_target);
        break;
    }

    update_force_U();
    
    int min_counter_shbox = 0;

    do
    {
        output_analysis(2, counter_system, time_system, counter_phyproc);

        ++count;

        relaxation(relaxation_option, dynamics_option, MD_MC_option, dt);

        min_counter_shbox += min_counter;

        counter_phyproc = count;

        if (abs(control_variable - control_variable_target) < cv_tol)
        {
            break;
        }

        L_box_old = L_box;
        switch (option)
        {
        case 1:
            output_analysis(2, counter_system, time_system, counter_phyproc);
                
            cout << "Relaxed with: " << min_counter << " steps in a box of " << L_box[0] << " at volume fraction: " << vol_frac << "\n";
            
            if ((abs((L_box - L_box_target).minCoeff()) <= dL_box.maxCoeff())) // || ((control_variable_target - control_variable) <= cv_tol))
            {
                L_box = L_box_target;
            }
            else
            {
                L_box -= dL_box;
            }
            update_vol_frac();
            control_variable = vol_frac;
            break;

        case 2:
            output_analysis(2, counter_system, time_system, counter_phyproc);
                
            cout << "Relaxed with: " << min_counter << " steps in a box of " << L_box[0] << " at volume fraction: " << vol_frac << "\n";

            if ((abs((L_box - L_box_target).minCoeff()) <= dL_box.maxCoeff())) // || ((control_variable_target - control_variable) <= cv_tol))
            {
                L_box = L_box_target;
            }
            else
            {
                dL_box = dvol_frac / vol_frac * L_box / 3.0;
                L_box -= dL_box;
            }
            update_vol_frac();
            control_variable = vol_frac;
            break;

        case 3:
            if (abs(bubbles[0].U / eff_U_0) / bubbles[0].N < U_tol)
            {
                output_analysis(2, counter_system, time_system, counter_phyproc);
                
                cout << "Relaxed with: " << min_counter << " steps in a box of " << L_box[0] << " at volume fraction: " << vol_frac << "\n";
                
                if ((abs((L_box - L_box_target).minCoeff()) <= dL_box.maxCoeff())) // || ((control_variable_target - control_variable) <= cv_tol))
                {
                    L_box = L_box_target;
                }
                else
                {
                    dvol_frac = dvol_frac_0;
                    dL_box = dvol_frac / vol_frac * L_box / 3.0;
                    L_box -= dL_box;
                }
            }
            else
            {
                bubbles[0].state_flag = 3;
                output_analysis(2, counter_system, time_system, counter_phyproc);
                
                dvol_frac /= 2.0;
                dL_box = dvol_frac / vol_frac * L_box / 3.0;
                L_box += dL_box;
            }
            update_vol_frac();
            control_variable = vol_frac;
            break;

        case 4:
            if (abs(bubbles[0].U / eff_U_0) / bubbles[0].N < U_tol)
            {   
                cout << "Relaxed with: " << min_counter << " steps in a box of " << L_box[0] << " at volume fraction: " << vol_frac << "\n";

                if (comp_decomp_flag == 2)
                {
                    dvol_frac /= 2.0;
                }
                // dvol_frac = dvol_frac_0;
                comp_decomp_flag = 1;
            }
            else
            {
                // dvol_frac /= 2.0;
                bubbles[0].state_flag = 3;
                if (comp_decomp_flag == 1)
                {
                    dvol_frac /= 2.0;
                }
                comp_decomp_flag = 2; 
            }

            output_analysis(2, counter_system, time_system, counter_phyproc);

            if (comp_decomp_flag == 1)
            {   
                dL_box = dvol_frac / vol_frac * L_box / 3.0;
                if ((abs((L_box - L_box_target).minCoeff()) <= dL_box.maxCoeff())) // || ((control_variable_target - control_variable) <= cv_tol))
                {
                    L_box = L_box_target;
                }
                else
                {
                    L_box -= dL_box;
                }
            }
            else if (comp_decomp_flag == 2)
            {
                // dvol_frac /= 2.0;
                dL_box = dvol_frac / vol_frac * L_box / 3.0;
                L_box += dL_box;
            }
            update_vol_frac();
            control_variable = vol_frac;
            break;

        default:
            update_vol_frac();
            control_variable = vol_frac;
            break;
        }

        if (dvol_frac < 1e-6)
        {
            cout << "dvol_frac is too small\n";
            break;
        }
        
        update_scale_position(L_box_old, L_box);
        
        update_force_U(type_force, type_bias_force);

    } while ((count <= SHBOX_STEP_MAX));

    min_counter = min_counter_shbox;
}

void update_scale_position(const VectorXd &L_old, const VectorXd &L)
{
    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].old_position = bubbles[i].position;
            bubbles[i].position = bubbles[i].position.array() * L.array() / L_old.array();
            bubbles[i].ddeltadisplacement = bubbles[i].position - bubbles[i].old_position;
            bubbles[i].deltadisplacement = bubbles[i].ddeltadisplacement;

            bubbles[i].bias_skin_distance += bubbles[i].deltadisplacement;

            bubbles[i].skin_distance += bubbles[i].deltadisplacement;

            bubbles[i].displacement += bubbles[i].deltadisplacement;
            bubbles[i].ddisplacement += bubbles[i].deltadisplacement;

            bubbles[i].ds = bubbles[i].deltadisplacement.norm();

            bubbles[i].contour += bubbles[i].ds;
            bubbles[i].dcontour += bubbles[i].ds;
        }
    }

    if (bubbles[0].nl_flag < 2)
    {
        bubbles[0].nl_flag = 2;
    }
    if ((bubbles[0].bl_flag > -1) && (bubbles[0].bl_flag < 2))
    {
        bubbles[0].bl_flag = 2;
    }
    bubbles[0].f_flag = 2;
    bubbles[0].t_flag = 2;
    h_flag = 3;
    r_flag = 2;
    bubbles[0].com_flag = 1;
    bubbles[0].z_flag = 1;

    if (bubbles[0].state_flag <= 0.5)
    {
        bubbles[0].state_flag = 4;
    }
    else if (bubbles[0].state_flag <= 1.5)
    {
        bubbles[0].state_flag = 3;
    }
    else if (bubbles[0].state_flag != 3)
    {
        bubbles[0].state_flag = 4;
    }
}