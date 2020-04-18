//
// Created by hiki on 5/12/19.
//

#ifndef MODEL_OF_HEISENBERG_OUTPUT_DATA_H
#define MODEL_OF_HEISENBERG_OUTPUT_DATA_H

#include "../spin.h"
#include "../configuration_of_system.h"

struct magnetization_output_data
{
    double magnetization  = 0.0;
    double sqr_of_magnetization = 0.0;
    double z_magnetization = 0.0;
    double plane_magnetization = 0.0;

};

struct output_data
{
    double e = 0;
    double e2 = 0;
    magnetization_output_data films_magnetization[number_of_films];
    spin layers_components[N * number_of_films];
};

void operator+=(magnetization_output_data& l_val, const magnetization_output_data&& r_val)
{
    l_val.magnetization += r_val.magnetization;
    l_val.z_magnetization += r_val.z_magnetization;
    l_val.plane_magnetization += r_val.plane_magnetization;
    l_val.sqr_of_magnetization += r_val.sqr_of_magnetization;
}

#endif //MODEL_OF_HEISENBERG_OUTPUT_DATA_H
