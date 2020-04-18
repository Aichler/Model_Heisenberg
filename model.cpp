#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include "components_of_Hamiltonian.cpp"
#include <unordered_set>
#include "data_job/data_job.h"

using namespace std;
using namespace chrono;

mt19937_64 generator;
uniform_real_distribution<double> z_rand(-1, 1);
uniform_real_distribution<double> teta_rand(0, 2 * M_PI);

unordered_set<int> rows_with_upper_interaction;
unordered_set<int> rows_with_lower_interaction;

inline spin get_random_spin() noexcept
{
    spin random_spin;

    random_spin.z = z_rand(generator);
    const double teta = teta_rand(generator);
    const double temp = sqrt(1 - random_spin.z * random_spin.z);
    random_spin.x = temp * cos(teta);
    random_spin.y = temp * sin(teta);

    return random_spin;
}


inline void algorithm_of_metropolis()
{
    auto time_of_start = high_resolution_clock::now();
    for (int i = 1; i < number_of_films; i++) {
        rows_with_upper_interaction.insert(i * N);
        rows_with_lower_interaction.insert((i * N) - 1);
    }

    output_data data_for_file[time_of_observation];
    
    for(t = 31; t <= 50 ; t++)
    {
    
        half_period =  50 + t*2;
        stat = 1.0 / (static_configuration_count * 10);
        stat1 = 1.0 / (2 * half_period);
        
        for (int i = 0; i < number_of_films ; ++i) {
            metastable_time[i] = 0.0;
            magnetization[i] = 0.0;
            configurations[i] = 0;
        }
        int m = 0;
    for (int static_configuration = 0; static_configuration < static_configuration_count; static_configuration++)
    {
    

        spin*** lattice = new spin** [L];
        for (size_t i = 0; i < L; i++) {
            lattice[i] = new spin* [L];
            for (size_t j = 0; j < L; j++) {
                lattice[i][j] = new spin[N * number_of_films];
            }
        }
        

        static constexpr bool is_multi_film_struct = number_of_films > 1;

        uniform_real_distribution < double > get_r(0, 1);
        uniform_int_distribution < int > get_i_or_j(0, L - 1);
        uniform_int_distribution < int > get_k(0, N * number_of_films - 1);
            field_period = 0;
            H = -0.55;
            
           // for (int i = 0; i < number_of_films ; ++i) {
           //     flags [i] = 0;
           // }
            
            for (int mcs = 0; mcs < 2 * half_period * 10; mcs++) {
            m = 0;
            
               // mcs = iteration;
                
                 if (field_period == (half_period)) {
                    H = - H;
                    field_period = 0;
                }
                 // std::cout << "H " << H << std::endl;
                field_period++;
                
               // if (H < 0){
                
                //  int z =+ mcs;
                  //std::cout << "Mcs " << z << std::endl;                                  
               // }
                
                
               // for (int i = 0; i < number_of_films ; ++i) {
                
               // if ((flags[i] == 0)) {
               // metastable_time[i] = metastable_time[i] + mcs;
               // std::cout   << ", metastable_time = " << metastable_time[i]  <<  " MSD "<<mcs << std::endl;
               // }
               // }
                
              //  for (int i = 0; i < number_of_films ; ++i) {
                //&& (flags[i] == 0
                //int z =+ mcs;
                //std::cout << "flags[i] " << flags[i] << std::endl;
                //    if ( (H < 0) && (flags[i] != 1)) {
                   //std::cout << "Mcs " << mcs << std::endl;
                        
                   //     metastable_time[i] = metastable_time[i] + mcs;
                  //      flags[i] = 1;
                   //    std::cout   << "H = " << H << ", metastable_time = " << metastable_time[i]  << "  configurations  " << configurations[i] <<  " MSD "<<mcs << std::endl;
                      // std::cout << "mcs = " << metastable_time[i] << std::endl;
                //      configurations[i]++;
                  //  }
               // }
                
                for (int elem_step = 0; elem_step < N * number_of_films * L * L; ++elem_step) {
                    int i = get_i_or_j(generator);
                    int j = get_i_or_j(generator);
                    int k = get_k(generator);

                    const spin selected_spin = lattice[i][j][k];
                    const double E1 = get_energy_of_spin_and_neighbours<is_multi_film_struct, system_lattice_type>
                            (lattice, rows_with_upper_interaction, rows_with_lower_interaction, i, j, k);
                    lattice[i][j][k] = get_random_spin();

                    const double E2 = get_energy_of_spin_and_neighbours<is_multi_film_struct, system_lattice_type>
                            (lattice, rows_with_upper_interaction, rows_with_lower_interaction, i, j, k);

                    const double delta_E = E2 - E1;

                    if (delta_E > 0) {
                        const double W = exp(-delta_E / ((temperature_start)));
                        const double rand_value = get_r(generator);
                        if (rand_value > W) {
                            lattice[i][j][k] = selected_spin;
                        }
                    }
                }

                //std::cout << "mcs = " << mcs << std::endl;
                Data_job::update_output_data(lattice, data_for_file[mcs], t, mcs);
            }

        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                delete[] lattice[i][j];
            }
            delete[] lattice[i];
        }
        delete[] lattice;
    }
      
    Data_job::data_output(data_for_file);

    auto time_of_end = high_resolution_clock::now();

    cout << "I SPEND " << duration_cast < seconds > (time_of_end - time_of_start).count() << " SECONDs OF MY LIFE" << endl;
    }
}

int main() {
    algorithm_of_metropolis();
}
