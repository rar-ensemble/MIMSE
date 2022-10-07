/*
 * File: otherf.cpp
 * Project: Quenching
 * File Created: Sunday, 1st July 2020 6:30:50 pm
 * Author: Amruthesh T (amru@seas.upenn.edu)
 * -----
 * Last Modified: Monday, 25th April 2022 10:33:49 pm
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

string convert_int(const int &number)
{
  stringstream ss; //create a stringstream
  ss << number;    //add number to the stream
  return ss.str(); //return a string with the contents of the stream
}

string convert_double(const double &number)
{
  stringstream ss; //create a stringstream
  ss << fixed << setprecision(4) << number;    //add number to the stream
  return ss.str(); //return a string with the contents of the stream
}