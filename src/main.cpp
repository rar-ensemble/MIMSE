/*
 * File: main.cpp
 * Project: Quenching
 * File Created: Monday, 18th June 2020 12:11:53 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Wednesday, 20th April 2022 9:59:23 pm
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

int main()
{
    time_initial = clock();

    // ios_base::sync_with_stdio(0);

    int avalanche_counter = 0;

    initialization();

    output_analysis(-2, counter_system, counter_system);

    if (type_pre_start == 1)
    {
        relaxation(1, 1, 1);
    }
    else if (type_pre_start == 2)
    {
        reset_system(1, 0, 0, 0, -1, 1);
    }
    
    output_analysis(-1, counter_system, counter_system);

    start_system();

    while (counter_system < (NUMBER_OF_QUENCHES + NUMBER_OF_EQB_STEPS))
    {
        ++counter_system;
        time_system += dt_reset_system;
        
        update_kbT(counter_system, time_system);

        reset_system(type_reset_system, type_relaxation_r, type_dynamics_r, type_MD_MC_r, dt_reset_system, -1, 1);

        output_analysis(1, counter_system, time_system);

        if (((type_reset_system == 2) && (counter_system % data_print_frequency == 0)) || (type_reset_system != 2))
        {
            quenching(type_quenching, type_relaxation_q, type_dynamics_q, type_MD_MC_q, dt_minimizer, type_controller, control_variable_final);

            output_analysis(1, counter_system, time_system, -1, -1, avalanche_counter);
        }
    }
    
    output_analysis(-1, counter_system, time_system);

    cout << "program run complete \n";

    time_final = clock();
    run_time = double(time_final - time_initial) / CLOCKS_PER_SEC;

    cout << "time taken :" << run_time << " seconds";
}