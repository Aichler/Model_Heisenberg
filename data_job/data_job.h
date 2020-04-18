//
// Created by hiki on 5/12/19.
//

#ifndef MODEL_OF_HEISENBERG_OUTPUT_JOB_H
#define MODEL_OF_HEISENBERG_OUTPUT_JOB_H

#include "./output_data.h"

class Data_job {

public:
    static void normalize_data(output_data &data, int const &count_of_steps);
    static void update_output_data(spin const *const *const *const sp, output_data &data, int t, int mcs);
    static magnetization_output_data measure_magnetization(spin components_for_magnetization, int size_of_lattice);
    static void data_output(output_data data_for_file[ITERATIONS_COUNT]);
};

#include "./data_job.cpp"
#endif //MODEL_OF_HEISENBERG_OUTPUT_JOB_H
