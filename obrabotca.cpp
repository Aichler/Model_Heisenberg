#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include "stdio.h"





 void algorithm_of_metropolis()
{
FILE *fo;
int L = 16;
int N = 6;
int number_of_films = 1;

double one, two;



 int time_of_relax = 1000;
 int time_of_observation = 1500;
 int static_configuration_count = 100;
 int ITERATIONS_COUNT = 40;



double** sm;

 double half_period = 0.0;
double field_period = 0.0;
int mcs = 0;
FILE *metatime;
double *metastable_time = new double [number_of_films];
int *flags = new int [number_of_films];
double *magnetization = new double [number_of_films];
int *configurations = new int [number_of_films];




      
    

    
    double H = -0.55;
    constexpr int t = 1;
     half_period =  14 + t*2;
     int z = 2 * half_period * 10;
     
     int n = z;//число строк
		int m = 1;//число столбцов на единицу больше числа пробелов
		double **x;
		x = new double*[n];
		for (int i = 0; i<n; i++) x[i] = new double[m];
        //cout <<
    //double *sm = new double[2][z];
        
                //if ((fo = fopen("L=16_N=6_films_count=2_H=1.dat", "r")) == NULL)
                //{
                //printf("Oshibka otkr fila.");
                //}
               // else
                //{
               // for (int i = 0; i < 2 * half_period * 10; i++)
               // {
               // fscanf(fo, "%d %d", &one, &two);
               // sm[0][i] = one;
               // sm[1][i] = two;
               // std::cout << "FFFFFFFFFFFFFFFF     " << two << std::endl;
                //}
               // }
                
                
                
                std::ifstream in("L=16_N=6_films_count=2_H=2.dat");
 
	if (in.is_open())
	{
 std::cout << "Ti che bly " << std::endl;
		//Если открытие файла прошло успешно
 
		//Вначале посчитаем сколько чисел в файле
		int count = 0;// число чисел в файле
		double temp;//Временная переменная
 std::cout << "Ti che bly 4" << std::endl;
		while (!in.eof())// пробегаем пока не встретим конец файла eof
		{
   //std::cout << temp << std::endl;
			in >> temp;//в пустоту считываем из файла числа
			count++;// увеличиваем счетчик числа чисел
		}
 std::cout << "Ti che bly 5" << std::endl;
		//Число чисел посчитано, теперь нам нужно понять сколько
		//чисел в одной строке
		//Для этого посчитаем число пробелов до знака перевода на новую строку 
 
		//Вначале переведем каретку в потоке в начало файла
		in.seekg(0, std::ios::beg);
		in.clear();
 std::cout << "Ti che bly 6" << std::endl;
		//Число пробелов в первой строчке вначале равно 0
		int count_space = 0;
		char symbol;
		//while (!in.eof())//на всякий случай цикл ограничиваем концом файла
		//{
			//теперь нам нужно считывать не числа, а посимвольно считывать данные
		//	in.get(symbol);//считали текущий символ
		//	if (symbol == ' ') count_space++;//Если это пробел, то число пробелов увеличиваем
		//	if (symbol == '\n') break;//Если дошли до конца строки, то выходим из цикла
	//	}
		//cout << count_space << endl;
 std::cout << "Ti che bly 7" << std::endl;
		//Опять переходим в потоке в начало файла
		in.seekg(0, std::ios::beg);
		in.clear();
 
		//Теперь мы знаем сколько чисел в файле и сколько пробелов в первой строке.
		//Теперь можем считать матрицу.
 
		
 
 std::cout << "Ti che bly 8" << std::endl;
		//Считаем матрицу из файла
		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++){
				in >> x[i][j];
        //std::cout << x[i][j] << std::endl;
 }
 }

 
		in.close();//под конец закроем файла
	}
	else
	{
		//Если открытие файла прошло не успешно
		std::cout << "Файл не найден.";
	}
                   
            
                
                
         for (int static_configuration = 0; static_configuration < static_configuration_count; static_configuration++)
    {
  
            field_period = 0;
            H = -0.55;
            
            for (int i = 0; i < number_of_films ; ++i) {
                flags [i] = 0;
            }
            
            for (mcs = 0; mcs < 2 * half_period * 10; mcs++) {
            int m = 0;
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
                
                
                
                
                for (int i = 0; i < number_of_films ; ++i) {
                //&& (flags[i] == 0
                //int z =+ mcs;
                //std::cout << "flags[i] " << flags[i] << std::endl;
                for(int j = 0; j <2 * half_period * 10; j++ ){
                    if ( (x[i][j] < 0) && (flags[i] == 0)) {
                //std::cout << "Mcs " << mcs << std::endl;
                std::cout << x[i][j] << std::endl;
                        
                        metastable_time[i] = metastable_time[i] + mcs;
                        flags[i] = 1;
                       std::cout   << "H = " << H << ", metastable_time = " << metastable_time[i]  << "  configurations  " << configurations[i] <<  " MSD "<< mcs << std::endl;
                       std::cout << "mcs = " << metastable_time[i] << std::endl;
                     configurations[i]++;
  }
        
    }
    }        
        
      }  
        
   }      
    
}
int main() {
std::cout << "Ti che bly " << std::endl;
    algorithm_of_metropolis();
}