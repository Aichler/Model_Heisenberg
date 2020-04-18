
#pragma argsused
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "randomc.h"
#include "mersenne.cpp"
#include <iostream>
#include <fstream>
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

int L = 16;
int z = 3;
int N = L*L*z, nstat = 1;
long configuration_max = 200;
double T = 4.5115 * 0.8;
double J_s = 1.0;
double J_b = 1.0;

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


double Neigbours( int ***sp, int i, int j, int k )
{


    int left = 0 ,right = 0 ,up = 0 ,down = 0 ,downz = 0 ,upz = 0 ;

	if (i == 0)
	{
        downz = 0;
	}
    else
	{
		downz = sp [i - 1][j][k];
	}

	if (i == z - 1)
    {
        upz = 0;
    }
    else
    {
        upz = sp[i + 1][j][k];
    }

    if (j == 0)
    {
        left = sp[i][L - 1][k];
    }
    else
    {
        left = sp[i][j - 1][k];
    }

    if (j == L - 1)
    {
        right = sp[i][0][k];
    }
    else
    {
	   right = sp[i][j + 1][k];
    }

    if (k == 0)
    {
        down = sp[i][j][L - 1];
    }
    else
    {
      down = sp[i][j][k - 1];
    }

    if (k == L - 1 )
    {
        up = sp[i][j][0];
    }
    else
    {
       up =  sp[i][j][k + 1];
    }

      return (double) left + right + up + down + downz + upz;

}

double Neigbours_layer(int ***sp, int i, int j, int k)
{
     int left = 0 ,right = 0, up = 0, down = 0;

    if (j == 0)
    {
        left = sp[i][L - 1][k];
    }
    else
    {
        left = sp[i][j - 1][k];
    }
    if (j == L - 1)
    {
        right = sp[i][0][k];
    }
    else
    {
       right = sp[i][j + 1][k];
    }
    if (k == 0)
    {
        down = sp[i][j][L - 1];
    }
    else
	{
      down = sp[i][j][k - 1];
    }
    if (k == L - 1 )
    {
        up = sp[i][j][0];
    }
    else
    {
       up =  sp[i][j][k + 1];
    }

      return  (double) left + right + up + down;
}

double Magnetic( int ***sp )
{
    int i, j, k, res = 0;

    for(i = 0;i < z;++i)
    {
        for(j = 0;j < L;++j)
        {
            for(k = 0;k < L;++k)
            {
                res += sp[i][j][k];
            }
        }
    }

    return (double) res;
}

double Magnetic_layer( int ***sp, int z)
{
    int res = 0;

    for(int i = 0; i < L; ++i)
    {
        for ( int j = 0 ; j < L; ++j )
        {
            res += sp[z][i][j];
        }
    }

    return (double) res;
}

int main()
{
    int i, j, a, ***sp, k, nt, ntstat = configuration_max/nstat, T_step;
    int field_period;
    long mcs, configuration, seed;
	double dE, r, w, stat, stat2, stat1;
    char fname[20];
    double field_constant,field_amplitude, metastable_time=0.0;
	int half_period = 0;
	int quantity_periods = 10;
    int periods;
    FILE *OUTn, *disp_Magnetic_result;


	double Magnetic_result=0;
	double	dispersion_magnetic_result;
double pow_Magnetic_result = 0.0;

    double *q_2 = new double [z];
    double *q_4 = new double [z];
    double *q_z = new double [z];
    double *q_eq = new double [z];
    double *m = new double [z];
    double *m_eq = new double [z];
	double *chi = new double [z];
    double *u = new double [z];
	double *tau = new double[z];

	ifstream metatime ("time.dat");
	for (i = 0; i < z; i++)
	{
			  metatime >>tau[i];
		 cout<<tau[i]<<endl;
    }
	sp = new int **[z];
    for(i = 0;i < z; ++i)
    {
        sp[i] = new int *[L];
        for(j = 0; j < L; ++j)
        {
            sp[i][j] = new int [L];
        }
	}

	double ***Magnetic_stat = new double **[z];
		 for(i = 0;i < z; ++i)
	{
		Magnetic_stat[i] = new double *[configuration_max];
		for ( j = 0; j < configuration_max; ++j)
		{
			Magnetic_stat[i][j] = new double [4 * (90 + 5*60) * quantity_periods ];
		}
	}
    seed = time(0);
    CRandomMersenne Mersenne(seed);
	field_constant = -0.55 * J_b;

    for(int t = 1; t <= 20 ; ++t)
    {
		half_period =  15 + t*2;
        stat = 1.0 / (configuration_max * quantity_periods);
        stat1 = 1.0 / (2 * half_period);
        printf("T = %.4lf\n", T);
		printf("t1/2 = %d\n", half_period);


        for (i = 0; i < z; ++i)
        {
            q_eq[i] = 0.0;
			q_2[i] = 0.0;
            q_4[i] = 0.0;
			chi[i] = 0.0;
        }
        periods = 0;
		for(configuration = 0;configuration < configuration_max; ++configuration)
        {
            for (i = 0; i < z; ++i)
            {
				m_eq[i] = 0.0;
            }

			printf("%ld\n", configuration);
            field_constant =  - 0.55 * J_b;
            cout<<field_constant<<endl;
			field_amplitude = 0.55;
            seed = seed + time(0) + (configuration%7)*Mersenne.IRandomX(10,20+2 * half_period * quantity_periods%Mersenne.IRandomX(1,2 * half_period * quantity_periods));
            Mersenne.RandomInit(seed);

            spin*** lattice = new spin** [L];
        for (size_t i = 0; i < L; i++) {
            lattice[i] = new spin* [L];
            for (size_t j = 0; j < L; j++) {
                lattice[i][j] = new spin[N * number_of_films];
            }
        }

        out << static_configuration << endl;
        uniform_real_distribution < double > get_r(0, 1);
        uniform_int_distribution < int > get_i_or_j(0, L - 1);
        uniform_int_distribution < int > get_k(0, N * number_of_films - 1);


            //for(i=0;i<z;++i)
			//{
            //    for(j=0;j<L;++j)
             //   {
             //       for(k=0;k<L;++k)
              //      {
              //          sp[i][j][k] = 1;
				//	}
              //  }
            //}

            field_period = 0;
			for(mcs = 1;mcs <= 2 * half_period * quantity_periods ;++mcs)
			{
				if(field_period == (half_period))
                {
                    field_constant = - field_constant;
					field_period = 0;
                }

				field_period++;
                for (i = 0; i < z ; ++i)
				{
					m[i] = Magnetic_layer(sp,i)/(L*L);
					m_eq[i] = m_eq[i] + m[i] * stat1;
					Magnetic_stat[i][configuration][mcs] = m[i];


                }

                if (mcs % (half_period * 2)== 0)
				{
                    periods++;
                    cout<<mcs<<endl<<endl;
                    for (i = 0; i < z ; ++i)
                    {
                    q_z[i] = abs(m_eq[i]);
                    m_eq[i] = 0.0;

                    if (periods >= 3)
                   {




                    q_eq[i] = q_eq[i] + q_z[i] * stat;
                    q_2[i] = q_2[i] + q_z[i] * q_z[i] * stat;
                    q_4[i] = q_4[i] + q_z[i] * q_z[i] * q_z[i] * q_z[i] * stat;
                    chi[i] = chi[i] + q_z[i] * q_z[i] * stat;
                   }
                    q_z[i] = 0.0;
					}
                }

                for(a = 0;a < N;++a)
                {
					int i = get_i_or_j(generator);
                    int j = get_i_or_j(generator);
                    int k = get_k(generator);

                    //dE = (get_energy_of_spin_and_neighbours(lattice,rows_with_upper_interaction, rows_with_lower_interaction,i,j,k);
                    //  lattice[i][j][k] = get_random_spin();
                    //if (i == 0)
                    //{
                     //  dE = 2.0 * sp[i][j][k] * (Neigbours_layer(sp,i,j,k) * J_s + field_constant + sp[i + 1][j][k] * J_b);
                    // }

                    //if (i == z - 1)
                    //{
                    //       dE = 2.0 * sp[i][j][k] * (Neigbours_layer(sp,i,j,k) * J_s + field_constant + sp[i - 1][j][k] * J_b);
                    // }

                    // if(dE <= 0.0)
                    //{
                    //    sp[i][j][k] = -sp[i][j][k];
                    // }

                    const double E1 = get_energy_of_spin_and_neighbours(lattice, rows_with_upper_interaction, rows_with_lower_interaction, i, j, k);
                    lattice[i][j][k] = get_random_spin();

                    const double E2 = get_energy_of_spin_and_neighbours(lattice, rows_with_upper_interaction, rows_with_lower_interaction, i, j, k);

                    const double delta_E = E2 - E1;

                    const double W = exp(-dE / ((T)));

                    if (delta_E > 0) {
                        const double rand_value = get_r(generator);
                        if (rand_value > W) {
                            lattice[i][j][k] = selected_spin;
                        }
                    }
                }
            }
		}

       for (i = 0; i < z; ++i)
		{

			 sprintf(fname,"Magnetization_%d_theta=%.4lf.dat",i,(double)half_period/tau[i]);

				OUTn = fopen(fname,"w+");
				disp_Magnetic_result = fopen("disp_Magnetic_result.dat", "w+");



			 for ( j = 1; j <= 2 * half_period * quantity_periods; ++j)
		{
			for ( k = 0; k < configuration_max; ++k)
			{
				cout<< Magnetic_stat[i][k][j];
				Magnetic_result += Magnetic_stat[i][k][j];
				pow_Magnetic_result += pow(Magnetic_stat[i][k][j] , 2);
			}
		dispersion_magnetic_result = (double)pow_Magnetic_result/configuration_max - ((double) Magnetic_result /configuration_max)*((double) Magnetic_result /configuration_max);
            Magnetic_result = (double) Magnetic_result /configuration_max;

		fprintf(OUTn,"%d\t%.4lf\n", j, Magnetic_result);
fprintf(disp_Magnetic_result,"%d\t%.4lf\n", j, Magnetic_result);
        Magnetic_result = 0.0;
	dispersion_magnetic_result = 0.0;

        }
        fclose(OUTn);
fclose(disp_Magnetic_result);

        }






        for (i = 0; i < z; ++i)
        {
			sprintf(fname,"Q_%d.dat",i);
            if (t == 0)
            {
				OUTn = fopen(fname,"w+");
            }
            else
            {
				OUTn = fopen(fname,"a+");
           }
            chi[i] = (chi[i] - q_eq[i] * q_eq[i]) * L * L;
            u[i] = 1.0 - ((q_4[i]) / (3 * q_2[i] * q_2[i]));
			fprintf(OUTn,"%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", (double)half_period/tau[i],q_eq[i],sqrt(q_2[i]-q_eq[i]*q_eq[i])/configuration_max,chi[i],u[i]);
            fclose(OUTn);
        }
    }

    delete(sp);
    return 0;
}

//---------------------------------------------------------------------------
