//
// Created by lotos on 29.03.2019.
//

#ifndef MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H
#define MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H

constexpr int L = 16;
constexpr int N = 32;
constexpr int number_of_films = 1;
constexpr double J = 1.0; //integral of interaction
constexpr double J2 = - 0.3 * J;
constexpr double delta =  0.0; //0.4;//0.68; //const of anisotropy
constexpr int time_of_relax = 1000;
constexpr int time_of_observation = 3000;
constexpr int static_configuration_count = 100;
constexpr int ITERATIONS_COUNT = 40;
constexpr char path[] = "../results/";
constexpr double H_INIT = 0.0;
constexpr double temperature_start = 0.8 * J;
constexpr double size_of_temperature_step = 0.1;
double stat = 0.0; 
double stat1 = 0.0;
double half_period = 0.0;
double field_period = 0.0;
double t = 0.0;
int mcs;
FILE *metatime;
double *metastable_time = new double [number_of_films];
    int *flags = new int [number_of_films];
    double *magnetization = new double [number_of_films];
    int *configurations = new int [number_of_films];

enum class lattice_types { simple_cubic_lattice, body_centered_cubic_lattice };
const lattice_types system_lattice_type = lattice_types::body_centered_cubic_lattice;

double H = H_INIT;



#endif //MODEL_OF_HEISENBERG_CONFIGURATION_OF_SYSTEM_H
