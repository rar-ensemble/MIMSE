/*
 * File: mathf.cpp
 * Project: Quenching
 * File Created: Tuesday, 26th June 2020 8:18:14 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 8:25:05 pm
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

void uniform_real_dist(const double &a, const double &b)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    uniform_real_distribution<double> distribution(a, b);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].radius = 0.5 * abs(distribution(generator_1));

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void real_gaussian_dist(const double &mean, const double &sigma)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    normal_distribution<double> distribution(mean, sigma);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].radius = 0.5 * abs(distribution(generator_1));

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void real_weibull_dist(const double &k, const double &lambda)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    weibull_distribution<double> distribution(k, lambda);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].radius = 0.5 * distribution(generator_1);

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void bidisperse_dist(const double &a, const double &b)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    // discrete_distribution<int> distribution{1, 1};

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            /* int index = distribution(generator_1);

            switch (index)
            {
            case 0:
                bubbles[i].radius = 0.5 * a;
                break;

            case 1:
                bubbles[i].radius = 0.5 * b;
                break;
            } */

            if (bubbles[i].real_index[0] <= 0.5 * bubbles[0].N)
            {
                bubbles[i].radius = 0.5 * a;
            }
            else
            {
                bubbles[i].radius = 0.5 * b;
            }

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void monodisperse_dist(const double &a)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    //discrete_distribution<int> distribution (1, 1);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].radius = 0.5 * a;

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void binary_dist(const double &a, const double &b, const double &x1)
{
    bubbles[0].volume = 0.0;

    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    //discrete_distribution<int> distribution (1, 1);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            if (bubbles[i].real_index[0] <= round(x1 * bubbles[0].N))
            {
                bubbles[i].type = 1;

                bubbles[i].mass = 1;

                bubbles[i].radius = 0.5 * a;
            }
            else
            {
                bubbles[i].type = 2;

                bubbles[i].mass = 1;

                bubbles[i].radius = 0.5 * b;
            }

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }
}

void swap_power_law(double &xmin, double &xmax)
{
    if (xmin > xmax)
    {
        swap(xmin, xmax);
    }

    bubbles[0].volume = 0.0;

    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            double p = distribution(generator_1);

            bubbles[i].radius = 1.0 / (sqrt(1.0 / gsl_pow_2(xmin) - p * (1.0 / gsl_pow_2(xmin) - 1.0 / gsl_pow_2(xmax))));

            bubbles[i].type = 1;

            bubbles[i].mass = 1;

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
        }
    }

    if (L_box[0] == -1)
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
    else
    {
        update_vol_frac();
    }

    double delta = sqrt(log(xmin/ xmax)/ ((xmin/ xmax - 1)/(xmin/ xmax + 1)) - 1);
}

void random_position_init()
{
    //srand(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    //default_random_engine generator_1(random_seed[min(int(random_seed.size()) - 1, random_seed_counter++)]);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            for (int k = 0; k < d; ++k)
            {
                bubbles[i].position[k] = L_box[k] * distribution(generator_2); //double(rand()) / RAND_MAX;
            }
            bubbles[i].velocity.setZero(d);
        }
    }
}

void indexing()
{
    for (int i = 1; i <= N_init; ++i)
    {
        bubbles[i].particle_index = i;
        bubbles[i].real_index.assign(1, i);
        if (i != 0)
        {
            bubbles[i].flag = 1;
        }
    }

    bubbles[0].state_flag = 4;
    bubbles[0].n_flag = 0;
    bubbles[0].nl_flag = 3;
    if (bubbles[0].bl_flag > -1)
    {
        bubbles[0].bl_flag = 3;
    }
    bubbles[0].f_flag = 1;
    bubbles[0].t_flag = 1;
    bubbles[0].z_flag = 1;
    bubbles[0].com_flag = 1;
}

void mass_initial(const int &mass_option)
{
    double rho = 1.0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            bubbles[i].mass *= mass_option;
            bubbles[0].mass += bubbles[i].mass;
            mass_init += (bubbles[i].volume * rho);
        }
    }

    bubbles[0].mass /= bubbles[0].N;
}

double mean_radius()
{
    double mean_r = 0.0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            mean_r += bubbles[i].radius;
        }
    }
    mean_r /= bubbles[0].N;

    return mean_r;
}

double mean_radius2()
{
    double sum_r2 = 0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            sum_r2 += gsl_pow_2(bubbles[i].radius);
        }
    }
    sum_r2 /= bubbles[0].N;

    return sum_r2;
}

double mean_radius3()
{
    double sum_r3 = 0.0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            sum_r3 += gsl_pow_3(bubbles[i].radius);
        }
    }
    sum_r3 /= bubbles[0].N;

    return sum_r3;
}

double max_radius()
{
    double max_r = 0.0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            max_r = max(max_r, bubbles[i].radius);
        }
    }

    return max_r;
}

double stddev_radius()
{
    double sum_r = 0;
    double sum_r2 = 0;

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            sum_r += bubbles[i].radius;
            sum_r2 += gsl_pow_2(bubbles[i].radius);
        }
    }

    double stddev = sqrt(sum_r2 / (bubbles[0].N - 1) - (bubbles[0].N / (bubbles[0].N - 1)) * gsl_pow_2(sum_r / bubbles[0].N));

    return stddev;
}

vector<int> largest_radius(const int &n)
{
    vector<int> largest_n(n, 1);

    vector<radii> R(bubbles[0].N);

    for (int i = 1; i <= N_init; ++i)
    {
        if (bubbles[i].flag == 1)
        {
            R[bubbles[i].real_index[0] - 1].radius = bubbles[i].radius;
            R[bubbles[i].real_index[0] - 1].particle_index = bubbles[i].particle_index;
        }
    }

    partial_sort(R.begin(), R.begin() + n, R.end(), greater<radii>());

    for (int k = 0; k < n; ++k)
    {
        largest_n[k] = R[k].particle_index;
    }

    return largest_n;
}

void init_run()
{
    force(-1);

    update_kbT();

    kbT = kbTi;

    bubbles[0].com_flag = 0;

    nl_para = L_box.norm() / 10;

    for (unsigned int k = 0; k < bl_para.size(); ++k)
    {
        if (bl_para[k] < 0)
        {
            bl_para[k] = INFINITY;
        }
    }

    switch (particle_fix_flag)
    {
    case 1:
        bubbles[largest_radius(1)[0]].nl_flag = 0;
        break;
    }

    if (bias_U_0 > 0)
    {
        bias_counter[1] = 1;
    }
    /* else if (bias_U_sigma)
    {
        bias_counter[1] = -1;
    } */
    bias_counter[0] = 0;

    if (well_U_0 > 0)
    {
        bias_counter[2] = 1;
    }

    update_force_U(type_force, type_bias_force);

    number_of_eigenvalues = (number_of_eigenvalues > d * bubbles[0].N) ? d * bubbles[0].N : number_of_eigenvalues;

    U_tol = U_tol_0;
    F_tol = F_tol_0;
    tolerances_set(type_dynamics_q, type_MD_MC_q, dt_minimizer);

    bubbles_thermal = bubbles;
}

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

vector<double> cardano_real_cubic_equation_solver(const double &p, const double &q)
{
    double D = -27 * gsl_pow_2(q) - 4 * gsl_pow_3(p);

    vector<double> root(3, 0);

    if ((D > 0) && (p < 0))
    {
        for (int k = 0; k < 3; ++k)
        {
            root[2 - k] = 2.0 * sqrt(-p / 3.0) * cos(1 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 2.0 * pi * k / 3.0);
        }
    }
    else if ((D < 0) && (p < 0))
    {
        root[0] = -2.0 * sgn(q) * sqrt(-p / 3.0) * cosh(1 / 3.0 * acosh(-3.0 * abs(q) / (2.0 * p) * sqrt(-3.0 / p)));
    }

    return root;
}