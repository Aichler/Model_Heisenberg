
//---------------------------------------------------------------------------

#pragma argsused
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "randomc.h"
#include "mersenne.cpp"
#include <iostream>

#include <iomanip>
#include <fstream>



#include "stdio.h"


using namespace std;
int L = 16;
int z = 6;
int N = L*L*z;
long mcs_max = 340, configuration_max = /*100*/ 100;
double T = 4.5115 * 0.8;
//double T = 3.8701;
//double T = 2.269;
double J_s = 1.0;
double J_b = 1.0;



double e1 = 1.0 / exp(1.0);



double Neigbours( int ***sp, int i, int j, int k ) {
    int left = 0 ,right = 0 ,up = 0 ,down = 0 ,downz = 0 ,upz = 0;

    if (i == 0) {
        downz = 0;
    } else {
        downz = sp [i - 1][j][k];
    }

    if (i == z - 1) {
        upz = 0;
    } else {
        upz = sp[i+1][j][k];
    }

    if (j == 0) {
        left = sp[i][L-1][k];
    } else {
        left = sp[i][j-1][k];
    }

    if (j == L - 1) {
        right = sp[i][0][k];
    } else {
       right =sp[i][j+1][k];
    }

    if (k == 0) {
        up = sp[i][j][L-1];
    } else {
      up =sp[i][j][k-1];
    }

    if (k == L - 1 ) {
        down = sp[i][j][0];
    } else {
       down =  sp[i][j][k+1];
    }

return (double) left + right + up + down + downz + upz;

}

double Neigbours_layer(int ***sp, int i, int j, int k) {
    int left = 0 ,right = 0, up = 0, down = 0;

    if (j == 0) {
        left = sp[i][L - 1][k];
    } else {
        left = sp[i][j - 1][k];
    }

    if (j == L - 1) {
        right = sp[i][0][k];
    } else {
       right = sp[i][j + 1][k];
    }

    if (k == 0) {
        down = sp[i][j][L - 1];
    } else {
      down = sp[i][j][k - 1];
    }

    if (k == L - 1 ) {
        up = sp[i][j][0];
    } else {
       up =  sp[i][j][k + 1];
    }

return  (double) left + right + up + down;

}

double Magnetic( int ***sp ) {
    int i, j, k, res=0;

    for (i = 0; i < z; ++i) {
        for (j = 0; j < L; ++j) {
            for( k = 0; k < L; ++k) {
                res = res + sp[i][j][k];
            }
        }
    }

return (double)res;

}

double Magnetic_layer( int ***sp, int z) {
    int res = 0;
    for (int i = 0; i < L; ++i) {
        for ( int j = 0 ; j < L; ++j ) {
            res+= sp[z][i][j];
        }
    }

return (double) res;

}

int main() {
    int i, j, a, ***sp, k, nt;
    int field_period;
    int mcs, configuration, seed;
    double dE, r, w[3], w1, stat, stat2;
    double *M, m, *X;
    double field_constant,field_amplitude;
    double Magnetic_result_surf;
    double Magnetic_result_bulk;
    double dispersion_magnetic_result_bulk, dispersion_magnetic_result_surf;
    double pow_magnetic_result_bulk, pow_magnetic_result_surf;
    double m_surf, m_bulk;
    int half_period = 0;
    int conf_stat;
    double Magnetic_result;

    FILE *FMAGNET, *field, *metatime;
    FILE *fm_bulk, *fm_surface, *fmtest, *D_fm_bulk, *D_fm_surf;
    double *metastable_time = new double [1];
    int *flags = new int [1];
    double *magnetization = new double [1];
    int *configurations = new int [1];
    char fname[353];
    FILE *OUTn;

	
    sp = new int **[z];
    for (i = 0; i < z; ++i) {
        sp[i] = new int *[L];
        for (j=0; j<L; ++j) {
            sp[i][j] = new int [L];
        }
    }

    double **Magnetic_stat_surf = new double *[configuration_max];
    //double Magnetic_stat_surf[100][10000]; // Ã¯Ã°Ã¨Ã¬Ã¥Ã° Ã±Ã²Ã Ã²Ã¨Ã·Ã¥Ã±ÃªÃ®Ã© Ã¯Ã Ã¬Ã¿Ã²Ã¨

   for (int i = 0; i < configuration_max; i++) {
        Magnetic_stat_surf[i] = new double [mcs_max*13];//!_!_!_!_!_!_!++++++++++++++!+!+!+!+!
   }

    double **Magnetic_stat_bulk = new double *[configuration_max];

    for (int i = 0; i < configuration_max; i++) {
        Magnetic_stat_bulk[i] = new double [mcs_max*13];//_!_!_!_!_!_!_!_+++++++++++!_!_!_!_!_!_
    }

    double ***Magnetic_stat = new double **[z];

    for (i = 0; i < z; ++i) {
        Magnetic_stat[i] = new double *[configuration_max];
        for ( j = 0; j < configuration_max; ++j) {
            Magnetic_stat[i][j] = new double [4 * half_period * 5];
        }
    }

    seed = time(0);
    CRandomMersenne Mersenne(seed);

    field = fopen("field.dat","w+");
	metatime = fopen ("time.dat", "w+");
    fm_bulk = fopen ("m_bulk.dat", "w+");
    fm_surface = fopen ("m_surface.dat", "w+");
    fmtest = fopen ("m_conf.dat", "w+");
	D_fm_bulk = fopen ("D_fm_bulk.dat", "w+");
D_fm_surf = fopen ("D_fm_surf.dat", "w+");


    for (int t = 1; t <= 50 ; ++t) {
        half_period = 14 + t*2;
        stat = 1.0 / (configuration_max * (2 * half_period * 5) * 2 * half_period);
        //printf("T = %.4lf\n", T);
        //printf("t1/2 = %d\n", half_period);
        int h = 2 * half_period * 10 + 2;

        for (i = 0; i < 1 ; ++i) {
            metastable_time[i] = 0.0;
            magnetization[i] = 0.0;
            configurations[i] = 0;
        }
        
        int n = h;//÷èñëî ñòðîê
		int m = 1;//÷èñëî ñòîëáöîâ íà åäèíèöó áîëüøå ÷èñëà ïðîáåëîâ
		double **x;
   double *mass = new double[h];
   
		x = new double*[n];
		for (int i = 0; i<n; i++) x[i] = new double[m];
        
        std::stringstream file_name;
        file_name << "L=16_N=6_films_count=2_H="<< t <<".dat";
        std::ifstream out("L=16_N=6_films_count=2_H=" + std::to_string(t) + ".dat");
       
	if (out.is_open())
	{
 std::cout << " Go 1 " << std::endl;
		//Åñëè îòêðûòèå ôàéëà ïðîøëî óñïåøíî
 
		//Âíà÷àëå ïîñ÷èòàåì ñêîëüêî ÷èñåë â ôàéëå
		int count = 0;// ÷èñëî ÷èñåë â ôàéëå
		double temp;//Âðåìåííàÿ ïåðåìåííàÿ
 std::cout << h << std::endl;
		while (!out.eof())// ïðîáåãàåì ïîêà íå âñòðåòèì êîíåö ôàéëà eof
		{
   //std::cout << temp << std::endl;
			out >> temp;//â ïóñòîòó ñ÷èòûâàåì èç ôàéëà ÷èñëà
			count++;// óâåëè÷èâàåì ñ÷åò÷èê ÷èñëà ÷èñåë
		}
 std::cout << "Go 2" << std::endl;
		//×èñëî ÷èñåë ïîñ÷èòàíî, òåïåðü íàì íóæíî ïîíÿòü ñêîëüêî
		//÷èñåë â îäíîé ñòðîêå
		//Äëÿ ýòîãî ïîñ÷èòàåì ÷èñëî ïðîáåëîâ äî çíàêà ïåðåâîäà íà íîâóþ ñòðîêó 
 
		//Âíà÷àëå ïåðåâåäåì êàðåòêó â ïîòîêå â íà÷àëî ôàéëà
		out.seekg(0, std::ios::beg);
		out.clear();
 std::cout << "Go 3" << std::endl;
		//×èñëî ïðîáåëîâ â ïåðâîé ñòðî÷êå âíà÷àëå ðàâíî 0
		int count_space = 0;
		char symbol;
		//while (!in.eof())//íà âñÿêèé ñëó÷àé öèêë îãðàíè÷èâàåì êîíöîì ôàéëà
		//{
			//òåïåðü íàì íóæíî ñ÷èòûâàòü íå ÷èñëà, à ïîñèìâîëüíî ñ÷èòûâàòü äàííûå
		//	in.get(symbol);//ñ÷èòàëè òåêóùèé ñèìâîë
		//	if (symbol == ' ') count_space++;//Åñëè ýòî ïðîáåë, òî ÷èñëî ïðîáåëîâ óâåëè÷èâàåì
		//	if (symbol == '\n') break;//Åñëè äîøëè äî êîíöà ñòðîêè, òî âûõîäèì èç öèêëà
	//	}
		//cout << count_space << endl;
 std::cout << "Go 4" << std::endl;
		//Îïÿòü ïåðåõîäèì â ïîòîêå â íà÷àëî ôàéëà
		out.seekg(0, std::ios::beg);
		out.clear();
 
		//Òåïåðü ìû çíàåì ñêîëüêî ÷èñåë â ôàéëå è ñêîëüêî ïðîáåëîâ â ïåðâîé ñòðîêå.
		//Òåïåðü ìîæåì ñ÷èòàòü ìàòðèöó.
 
		
 
 std::cout << "Go 5" << std::endl;
		//Ñ÷èòàåì ìàòðèöó èç ôàéëà
		for (int i = 0; i < m; i++){
			for (int j = 0; j < n; j++){
				out >> x[i][j];
        //std::cout << x[i][j] << std::endl;
 }
 }
 
 
 

 
		out.close();//ïîä êîíåö çàêðîåì ôàéëà
	}
	else
	{
		//Åñëè îòêðûòèå ôàéëà ïðîøëî íå óñïåøíî
		std::cout << "Ôàéë íå íàéäåí.";
	}
        
        
        std::cout << "Go 6" << std::endl;
        
     
        for (configuration = 0; configuration < 100; ++configuration) {
            
            field_constant = -0.55;
            field_amplitude = 0.55;
           

            for (i = 0; i < 1 ; ++i) {
                flags [i] = 0;
            }
            //std::cout << "x   " << x[1][configuration] << std::endl;

            for (int iteration = 0; iteration <= h ; ++iteration) {
            //std::cout << "x   " << x[0][iteration] << std::endl;
            
            
               for (i = 0; i < 1; ++i) {
                //for(int j = 0 ; j < h ; j ++ ){
                //std::cout << "x   " << x[i][j] << std::endl;
                    if ( ( x[0][iteration] <= 0.0) && (flags[i] == 0)) {
                    
                        configurations[i]++;
                        metastable_time[i] = metastable_time[i] + iteration;
                        flags[i] = 1;
                         std::cout <<  "  MCS " << metastable_time[i]  << " configurations  " << configurations[i]  << "  LKO  " << iteration  << std::endl;
                    }
                  //  }
                }

                

                fprintf(field,"%ld\t%.4lf\n", mcs, field_constant/field_amplitude);

                
                }
            }
            
            double tau = 0.0;
            
            
            for (i = 0; i < 1 ; ++i) {
            tau = metastable_time[i]/configurations[i];
            std::cout <<  "  MCS " << metastable_time[i]  << " configurations  " << configurations[i]  <<  " tau "<< tau <<  std::endl;
        }
        
        
        for (i = 0; i < 1; ++i)
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
			fprintf(OUTn,"%.4lf\t%.4lf\n", half_period/tau);
            fclose(OUTn);
        }
            
             //std::cout <<  "  MCS " << metastable_time[i]  << " configurations  " << configurations[i]  << std::endl;
        }
}

     /*   Magnetic_result_surf = 0.0;
        Magnetic_result_bulk = 0.0;
	dispersion_magnetic_result_bulk = 0.0;
	dispersion_magnetic_result_surf = 0.0;
	pow_magnetic_result_bulk = 0.0;
	pow_magnetic_result_surf = 0.0;
        for ( i = 1; i <=mcs_max; ++i) {

            for (int j = 0; j < configuration_max; ++j)
            {
                Magnetic_result_surf += Magnetic_stat_surf[j][i];//!!!!!!!!!!!!_!_!_!_!Ã³Ã¡Ã°Ã Ã²Ã¼ Ã­Ã³Ã«Ã¼
                Magnetic_result_bulk += Magnetic_stat_bulk[j][i];//!!!!!!!!!!!_!_!_!_!_
		pow_magnetic_result_surf += pow(Magnetic_stat_surf[j][i], 2);
		pow_magnetic_result_bulk += pow(Magnetic_stat_bulk[j][i] , 2);
            }

	Magnetic_result_bulk = (double)Magnetic_result_bulk/configuration_max;
	dispersion_magnetic_result_bulk =  (double)pow_magnetic_result_bulk/ 
configuration_max - ((double)Magnetic_result_bulk/configuration_max)*((double)Magnetic_result_bulk/configuration_max);

            Magnetic_result_surf = (double)Magnetic_result_surf/configuration_max;
dispersion_magnetic_result_surf =  (double)pow_magnetic_result_surf/ 
configuration_max - ((double)Magnetic_result_surf/configuration_max)*((double)Magnetic_result_surf/configuration_max);

		
             // Ã¬Ã®Ã¦Ã¥Ã² Ã¡Ã»Ã²Ã¼ Ã§Ã¤Ã¥Ã±Ã¼ Ã¤Ã®Ã«Ã¦Ã­Ã® Ã¡Ã»Ã²Ã¼ Magnetic_result_bulk
	    fprintf(D_fm_bulk,"%d\t%.8lf\n", i, dispersion_magnetic_result_bulk);
		fprintf(D_fm_surf,"%d\t%.8lf\n", i, dispersion_magnetic_result_surf);
            fprintf(fm_bulk,"%d\t%.8lf\n", i, Magnetic_result_bulk);
            fprintf(fm_surface,"%d\t%.8lf\n", i, Magnetic_result_surf);
pow_magnetic_result_surf = 0.0;
pow_magnetic_result_bulk = 0.0;
            Magnetic_result_surf = 0.0;
            Magnetic_result_bulk = 0.0;
        }
std::cout << "Yse" << endl;*/

        

//---------------------------------------------------------------------------
