//
// Created by hiki on 5/12/19.
//

#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>

void Data_job::normalize_data(output_data &data, int const &count_of_steps)
{
    data.e  /= static_cast<double>(count_of_steps);
    data.e2 /= static_cast<double>(count_of_steps);
    for(int i = 0; i < number_of_films; i++)
    {
        data.films_magnetization[i].magnetization /= count_of_steps;
        data.films_magnetization[i].z_magnetization /= count_of_steps;
        data.films_magnetization[i].plane_magnetization /= count_of_steps;
        data.films_magnetization[i].sqr_of_magnetization /= count_of_steps;
    }



    for(int i = 0; i < number_of_films * N; i++)
    {
        data.layers_components[i].x /= count_of_steps;
        data.layers_components[i].y /= count_of_steps;
        data.layers_components[i].z /= count_of_steps;
    }

}

magnetization_output_data Data_job::measure_magnetization(spin components_for_magnetization, int size_of_lattice)
{
    magnetization_output_data data;
    double sqr_of_projections_x =  components_for_magnetization.x * components_for_magnetization.x;
    double sqr_of_projections_y =  components_for_magnetization.y * components_for_magnetization.y;
    double sqr_of_projections_z =  components_for_magnetization.z * components_for_magnetization.z;


    data.magnetization = fabs(static_cast<double>(1.0 / (size_of_lattice))
                              * pow( (sqr_of_projections_x + sqr_of_projections_y +  sqr_of_projections_z) , 0.5));

    data.z_magnetization = static_cast<double>(1.0 / (size_of_lattice)) * components_for_magnetization.z;

    data.plane_magnetization = fabs(static_cast<double>(1.0 / (size_of_lattice))
                                    * pow((sqr_of_projections_x + sqr_of_projections_y) , 0.5));

    return data;
}

void Data_job::update_output_data(spin const *const *const *const sp, output_data &data, int t, int mcs)
{
    double  energy = 0.0;
    spin components_for_layer_magnetization[number_of_films];
    spin layers_components[number_of_films * N];

    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int layer = 0; layer < number_of_films; layer++ )
            {
                for (int spin = 0; spin < N; spin++) {
                    const int k = (spin + (layer * N));

                    layers_components[k].x += sp[i][j][k].x;
                    layers_components[k].y += sp[i][j][k].y;
                    layers_components[k].z += sp[i][j][k].z;

                    components_for_layer_magnetization[layer].x += sp[i][j][k].x;
                    components_for_layer_magnetization[layer].y += sp[i][j][k].y;
                    components_for_layer_magnetization[layer].z += sp[i][j][k].z;

                }
            }
        }
    }

    for(int i = 0; i < N * number_of_films; i++)
    {
        data.layers_components[i].x +=  ( layers_components[i].x /  (L * L) );
        data.layers_components[i].y +=  ( layers_components[i].y /  (L * L) );
        data.layers_components[i].z +=  ( layers_components[i].z /  (L * L) );
    }



    for(int i = 0; i < number_of_films; i++)
    {
        data.films_magnetization[i] +=  measure_magnetization(components_for_layer_magnetization[i], N * L * L);

    }
}


void  Data_job::data_output(output_data data_for_file[ITERATIONS_COUNT])
{

    std::cout << "mcs = " << ITERATIONS_COUNT << std::endl; 
    metatime = fopen ("time.dat", "w+");

   // for (int i = 0; i < number_of_films ; ++i) {
  //          metastable_time[i] = 0.0;
   //         magnetization[i] = 0.0;
    //        configurations[i] = 0;
     //   }


        for (int i = 0; i < number_of_films ; ++i) {
                flags [i] = 0;
            }

    for (int step = 0; step < time_of_observation; step++)
    {
        normalize_data(data_for_file[step],  static_configuration_count);
    }

    std::stringstream file_name;
    file_name << "L=" << L << "_N=" << N << "_films_count=" << number_of_films << "_H="<< t<<".dat";
    std::cout << "t = " << t << ", mcs " << mcs << std::endl;

    std::ofstream output_file(file_name.str(), std::ios_base::trunc);
    //output_file << std::fixed << std::setprecision(15) << "steps" << "\t";

    //for(int i = 0; i < 1; i++)
    //{
    //    output_file << std::fixed << std::setprecision(15)  << "M_z_" << i << "\t" ;
    //}

    //for (int i = 0; i < number_of_films * N; i++)
    //{
    //    output_file << std::fixed << std::setprecision(15)  << "X_of_layer" << i << "\t" << "Y_of_layer" << i << "\t" << "Z_of_layer" << i << "\t";
   // }

    output_file << std::endl;


    for (int step = 0; step < 2 * half_period * 10; step++)
    {
//        output_file << std::fixed << std::setprecision(15) << step << "\t";
        for(int i = 0; i < 1; i++)
        {
        
        
        
            output_file << std::fixed << std::setprecision(15)
           // << data_for_file[step].films_magnetization[i].magnetization << "\t"
            << data_for_file[step].films_magnetization[i].z_magnetization << "\t";
            //<< data_for_file[step].films_magnetization[i].plane_magnetization << "\t";

        }

       // for(int i = 0; i < number_of_films * N; i++ )
       // {
        //    output_file << std::fixed << std::setprecision(15)
        //    << data_for_file[step].layers_components[i].x << "\t"
        //    << data_for_file[step].layers_components[i].y << "\t"
        //    << data_for_file[step].layers_components[i].z << "\t";
       // }
        output_file << std::endl;
    }

    output_file.close();
}

//
//#else
//
//for (int h_step = 0; h_step < ITERATIONS_COUNT; h_step++)
//{
//H = H_INIT +  (h_step * size_of_h_step);
//
//output_file << std::fixed << std::setprecision(15) << H << "\t";
//for(int i = 0; i < number_of_films; i++)
//{
//output_file << std::fixed << std::setprecision(15)
//<< data_for_file[h_step].films_magnetization[i].magnetization << "\t"
//<< data_for_file[h_step].films_magnetization[i].z_magnetization << "\t"
//<< data_for_file[h_step].films_magnetization[i].plane_magnetization << "\t";
//
//}
//
//for(int i = 0; i < number_of_films * N; i++ )
//{
//output_file << std::fixed << std::setprecision(15)
//<< data_for_file[h_step].layers_components[i].x << "\t"
//<< data_for_file[h_step].layers_components[i].y << "\t"
//<< data_for_file[h_step].layers_components[i].z << "\t";
//}
//
//output_file << std::endl;
//}
//
//
//#endif