//
// Created by lotos on 29.03.2019.
//

#ifndef MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H
#define MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H

#include "spin.h"
#include "configuration_of_system.h"
#include <unordered_set>
#include <iostream>

inline double get_energy_of_two_spin(spin const &spin1, spin const &spin2) noexcept
{
    return (-J *  ( ( (1.0 - delta) *  ((spin1.x * spin2.x) + (spin1.y * spin2.y)) ) +  ( (  spin1.z * spin2.z))  )) + (H * spin1.z);
}

inline double get_energy_of_two_spin_whith_rkky_interaction(spin const &spin1, spin const &spin2) noexcept
{
    return -J2 * (((1.0 - delta) * ((spin1.x * spin2.x) + (spin1.y * spin2.y))) + (spin1.z * spin2.z)) + (H* spin1.z);
}

//************************************************************************************************************************\\
// SIMPLE CUBIC LATTICE


template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_with_uppder_interaction
       (spin const *const *const *const sp, int const &i, int const &j, int const &k,
        typename  std::enable_if_t<lattice_type == lattice_types::simple_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    i == 0
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[L - 1][j][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i - 1][j][k]);

    i == L - 1
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[0][j][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i + 1][j][k]);

    j == 0
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][L - 1][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j - 1][k]);

    j == L - 1
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][0][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j + 1][k]);

    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i][j][k - 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);

    return energy;
}

template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_with_lower_interaction
       (spin const *const *const *const sp, int const &i, int const &j, int const &k,
               typename  std::enable_if_t<lattice_type == lattice_types::simple_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    i == 0
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[L - 1][j][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i - 1][j][k]);

    i == L - 1
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[0][j][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i + 1][j][k]);

    j == 0
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][L - 1][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j - 1][k]);

    j == L - 1
        ? energy += get_energy_of_two_spin(sp[i][j][k], sp[i][0][k])
        : energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j + 1][k]);

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i][j][k + 1]);

    return energy;
}

template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_without_rrky_interaction
       (spin const *const *const *const sp, int const &i, int const &j, int const &k,
               typename std::enable_if_t<lattice_type == lattice_types::simple_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    const int i_minus_1 = i == 0 ? L - 1 : i - 1;
    const int i_plus_1  = i == L -1 ? 0 : i + 1;

    const int j_minus_1 = j == 0 ? L -1 : j - 1;
    const int j_plus_1  = j == L - 1? 0 : j + 1;

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j_minus_1][k]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j_plus_1][k]);

    if(k != 0)
    {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k - 1]);
    }

    if(k != (N * number_of_films) - 1)
        {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i][j][k + 1]);
    }

    return energy;
}

//************************************************************************************************************************\\


//************************************************************************************************************************\\
// BODY CENTRED CUBIC LATTICE


template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_without_rrky_interaction
        (spin const *const *const *const sp, int const &i, int const &j, int const &k,
         typename std::enable_if_t<lattice_type == lattice_types::body_centered_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    const int i_minus_1 = i == 0 ? L - 1 : i - 1;
    const int i_plus_1  = i == L -1 ? 0 : i + 1;

    const int j_minus_1 = j == 0 ? L -1 : j - 1;
    const int j_plus_1  = j == L - 1? 0 : j + 1;



    if (k != 0) {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_minus_1][k - 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_minus_1][k - 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_plus_1][k - 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_plus_1][k - 1]);
    }

    if (k != (N * number_of_films) - 1) {
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_minus_1][k + 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_minus_1][k + 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_plus_1][k + 1]);
        energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_plus_1][k + 1]);
    }

    return energy;
}


template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_with_lower_interaction
        (spin const *const *const *const sp, int const &i, int const &j, int const &k,
         typename std::enable_if_t<lattice_type == lattice_types::body_centered_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    const int i_minus_1 = i == 0 ? L - 1 : i - 1;
    const int i_plus_1  = i == L -1 ? 0 : i + 1;

    const int j_minus_1 = j == 0 ? L -1 : j - 1;
    const int j_plus_1  = j == L - 1? 0 : j + 1;

    const int k_minus_1  = k == 0?  ((N * number_of_films) - 1 ) : k - 1;
    const int k_plus_1  = k == ((N * number_of_films) - 1 )? 0 : k + 1;

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_minus_1][k - 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_minus_1][k - 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_plus_1][k - 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_plus_1][k - 1]);

    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_minus_1][j_minus_1][k + 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_plus_1][j_minus_1][k + 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_minus_1][j_plus_1][k + 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_plus_1][j_plus_1][k + 1]);

    return energy;
}

template <lattice_types lattice_type>
double get_energy_of_spin_and_neighbours_with_uppder_interaction
        (spin const *const *const *const sp, int const &i, int const &j, int const &k,
         typename std::enable_if_t<lattice_type == lattice_types::body_centered_cubic_lattice>* = nullptr) noexcept
{
    double energy = 0;

    const int i_minus_1 = i == 0 ? L - 1 : i - 1;
    const int i_plus_1  = i == L -1 ? 0 : i + 1;

    const int j_minus_1 = j == 0 ? L -1 : j - 1;
    const int j_plus_1  = j == L - 1? 0 : j + 1;


    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_minus_1][j_minus_1][k - 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_plus_1][j_minus_1][k - 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_minus_1][j_plus_1][k - 1]);
    energy += get_energy_of_two_spin_whith_rkky_interaction(sp[i][j][k], sp[i_plus_1][j_plus_1][k - 1]);

    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_minus_1][k + 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_minus_1][k + 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_minus_1][j_plus_1][k + 1]);
    energy += get_energy_of_two_spin(sp[i][j][k], sp[i_plus_1][j_plus_1][k + 1]);

    return energy;
}

//************************************************************************************************************************\\





template <bool is_multi_film_struct, lattice_types lattice_type>
double get_energy_of_spin_and_neighbours
        (spin const *const *const *const sp, std::unordered_set<int> const& rows_with_upper_interaction,
         std::unordered_set<int> const& rows_with_lower_interaction,  int const &i, int const &j, int const &k,
         typename std::enable_if_t<is_multi_film_struct>* = nullptr
         ) noexcept
{
    if(rows_with_upper_interaction.find(k) != rows_with_upper_interaction.end())
    {
        return get_energy_of_spin_and_neighbours_with_uppder_interaction<lattice_type>(sp, i, j, k);
    }

    if(rows_with_lower_interaction.find(k) != rows_with_lower_interaction.end())
    {
        return get_energy_of_spin_and_neighbours_with_lower_interaction<lattice_type>(sp, i, j, k);
    }

    return get_energy_of_spin_and_neighbours_without_rrky_interaction<lattice_type>(sp, i, j, k);

}

// TODO: remove sets from signature
template <bool is_multi_film_struct, lattice_types lattice_type>
inline double get_energy_of_spin_and_neighbours
        (spin const *const *const *const sp, std::unordered_set<int> const& rows_with_upper_interaction,
         std::unordered_set<int> const& rows_with_lower_interaction,  int const &i, int const &j, int const &k,
         typename std::enable_if_t<!is_multi_film_struct>* = nullptr
         ) noexcept
{
    return get_energy_of_spin_and_neighbours_without_rrky_interaction<lattice_type>(sp, i, j, k);
}


#endif //MODEL_OF_HEISENBERG_COMPONENTS_OF_HAMILTONIAN_H
