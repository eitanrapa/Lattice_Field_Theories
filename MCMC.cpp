#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

/*
Functions
*/

double magnetization(int **state, int matrix_size)
{
    double sum = 0;
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            sum = sum + state[i][j];
        }
    }
    return sum/pow(matrix_size,2.0);
}

double energy_state(double h, double j, int pos, int neig1, int neig2,
int neig3, int neig4)
{
    double energy_j = 0.0;
    energy_j = energy_j + pos*neig1;
    energy_j = energy_j + pos*neig2;
    energy_j = energy_j + pos*neig3;
    energy_j = energy_j + pos*neig4;
    return -j*energy_j - h*pos;
}

vector<double> MCMC(int **state, 
double h, double j, double beta, int matrix_size, int **later_state)
{
    vector<double> values;
    int** current_state = state;
    values.push_back(magnetization(current_state, matrix_size));
    int counter = 0;
    for (int i = 1; i < 10000*pow(matrix_size,2.0); i++)
    {       
        int RandIndex_i = rand() % matrix_size;
        int RandIndex_j = rand() % matrix_size;
        int pos = current_state[RandIndex_i][RandIndex_j];
        int neig1 = 0;
        int neig2 = 0;
        int neig3 = 0;
        int neig4 = 0; 
        if (RandIndex_i - 1 > -1)
        {
            neig1 = current_state[RandIndex_i - 1][RandIndex_j];
        }
        if (RandIndex_i + 1 < matrix_size)
        {
            neig2 = current_state[RandIndex_i + 1][RandIndex_j];
        }
        if (RandIndex_j - 1 > -1)
        {
            neig3 = current_state[RandIndex_i][RandIndex_j - 1];
        }
        if (RandIndex_j + 1 > matrix_size)
        {
            neig4 = current_state[RandIndex_i][RandIndex_j + 1];
        }
        double energy1 = beta/abs(beta)*
        energy_state(h, j, -1*pos, neig1, neig2, neig3,neig4);
        double energy2 = beta/abs(beta)*
        energy_state(h, j, pos, neig1, neig2, neig3,neig4);
        if (energy1 < energy2)
        {
            current_state[RandIndex_i][RandIndex_j] = -1*pos;
        }
        else
        {
            double random = ((double) rand() / (RAND_MAX));
            if (random < exp(-beta*(energy1 - energy2)))
            {
                current_state[RandIndex_i][RandIndex_j] = -1*pos;
            }
        }
        if (counter % (int) pow(matrix_size,2.0) == 0)
        {
            values.push_back(magnetization(current_state, matrix_size));
        }
        if (counter == 9000*pow(matrix_size,2.0))
        {
            for(int n = 0; n < matrix_size; n++)
            {
                for(int m = 0; m < matrix_size; m++)
                {
                    later_state[n][m] = current_state[n][m];
                }
            }
        } 
        counter = counter + 1;
    }
    return values;
}

void magnetization_density_average(int **state, double h, double j, 
double beta, int matrix_size, double *returnArray, int **later_state)
{
    vector<double> values = MCMC(state, h, j, beta, matrix_size, later_state);
    double burn_time = 500;
    vector<double> burn_values = 
    vector<double>(values.begin() + burn_time, values.end());
    double sum_of_elems = 
    accumulate(burn_values.begin(), burn_values.end(), 0.0);
    double average = sum_of_elems/burn_values.size();
    double std = 0.0;
    for (int i = 0; i < (signed) burn_values.size(); i++)
    {   
        std = std + pow(burn_values[i] - average, 2.0);
    }
    std = sqrt((1.0/burn_values.size())*std)*(1/sqrt(burn_values.size()));
    returnArray[0] = average;
    returnArray[1] = std;
}

void magnetization_density_average_sq(int **state, double h, double j, 
double beta, int matrix_size, double *returnArray, int **later_state)
{
    vector<double> values = MCMC(state, h, j, beta, matrix_size, later_state);
    double burn_time = 500;
    vector<double> burn_values = 
    vector<double>(values.begin() + burn_time, values.end());
    transform(burn_values.begin(), burn_values.end(),
    burn_values.begin(), [](double f)->double { return f * f; });
    double sum_of_elems = 
    accumulate(burn_values.begin(), burn_values.end(), 0.0);
    double average = sum_of_elems/burn_values.size();
    returnArray[0] = average;
    returnArray[1] = 0;
}

void susceptibility_average(int **state, double h, double j, 
double beta, int matrix_size, double *returnArray, int **later_state)
{
    vector<double> values = MCMC(state, h, j, beta, matrix_size, later_state);
    double burn_time = 500;
    vector<double> burn_values = 
    vector<double>(values.begin() + burn_time, values.end());
    
    double sum_of_elems_1 = 
    accumulate(burn_values.begin(), burn_values.end(), 0.0);
    double average_1 = sum_of_elems_1/burn_values.size();
    
    transform(burn_values.begin(), burn_values.end(),
    burn_values.begin(), [](double f)->double { return f * f; });
    double sum_of_elems_2 = 
    accumulate(burn_values.begin(), burn_values.end(), 0.0);
    double average_2 = sum_of_elems_2/burn_values.size();
    returnArray[0] = pow(matrix_size,2.0)*beta*(average_2 - pow(average_1,2.0));
    returnArray[1] = 0;
}

void binder_cumulant(int **state, double h, double j, 
double beta, int matrix_size, double *returnArray, int **later_state)
{
    vector<double> values = MCMC(state, h, j, beta, matrix_size, later_state);
    double burn_time = 500;
    vector<double> burn_values = 
    vector<double>(values.begin() + burn_time, values.end());
    vector<double> burn_values_2 = 
    vector<double>(values.begin() + burn_time, values.end());
    
    transform(burn_values.begin(), burn_values.end(),
    burn_values.begin(), [](double f)->double { return f * f; });
    transform(burn_values_2.begin(), burn_values_2.end(),
    burn_values_2.begin(), [](double f)->double { return f * f * f * f; });
    
    double sum_of_elems_sq = 
    accumulate(burn_values.begin(), burn_values.end(), 0.0);
    double average_sq = sum_of_elems_sq/burn_values.size();
    
    double sum_of_elems_4th = 
    accumulate(burn_values_2.begin(), burn_values_2.end(), 0.0);
    double average_4th = sum_of_elems_4th/burn_values_2.size();
    
    returnArray[0] = 1 - (average_4th)/(3*pow(average_sq,2.0));
    returnArray[1] = 0;
}
/*
Main function
*/

int main() 
{

    srand(time(NULL));
    
    /*For <m> as a function of beta*/
    /*
    double h = 1.0;
    double j = 0.0;
    int beta_size = 21;
    double beta[beta_size] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,
    5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10};
    int matrix_size = 10;
    //Averages and Standard dvs. 
    double values[beta_size][2];
    
    int **later_state;
    later_state= new int *[matrix_size];
    for(int i = 0; i < matrix_size; i++)
    {
        later_state[i] = new int[matrix_size];
    }
    for(int i = 0; i < matrix_size; i++)
    {
        for(int j = 0; j < matrix_size; j++)
        {
            later_state[i][j] = 0;
        }
    }
    
    for (int i = 0; i < beta_size; i++)
    {
        int **initial_state;
        initial_state = new int *[matrix_size];
        for(int i = 0; i < matrix_size; i++)
        {
            initial_state[i] = new int[matrix_size];
        }
        for(int i = 0; i < matrix_size; i++)
        {
            for(int j = 0; j < matrix_size; j++)
            {
                int ranchoice = (rand() > RAND_MAX/2) ? -1 : 1;
                initial_state[i][j] = ranchoice; 
            }
        }
        magnetization_density_average(initial_state, 
        h, j, beta[i], matrix_size, values[i],later_state);
    }
    
    ofstream myfile ("Check_in_data.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < beta_size; i++) 
        {
            myfile << to_string(beta[i]) + " " + to_string(values[i][0]) +
            " " + to_string(values[i][1]);
            myfile << "\n";
        }
    myfile.close();
    }
    */
   
    /*Warm up*/
    
    
    double h = 0.0;
    double j = 1.0;
    double beta= -32.0;
    double values[2];
    int matrix_size = 100;
    
    int **later_state;
    later_state= new int *[matrix_size];
    for(int i = 0; i < matrix_size; i++)
    {
        later_state[i] = new int[matrix_size];
    }
    for(int i = 0; i < matrix_size; i++)
    {
        for(int j = 0; j < matrix_size; j++)
        {
            later_state[i][j] = 0;
        }
    }
    
    int **initial_state;
    initial_state = new int *[matrix_size];
    for(int i = 0; i < matrix_size; i++)
    {
       initial_state[i] = new int[matrix_size];
    }
    for(int i = 0; i < matrix_size; i++)
    {
        for(int j = 0; j < matrix_size; j++)
        {
            int ranchoice = (rand() > RAND_MAX/2) ? -1 : 1;
            initial_state[i][j] = ranchoice; 
        }
    }
    
    magnetization_density_average(initial_state, 
        h, j, beta, matrix_size, values, later_state);
    
    ofstream myfile ("Check_in_data.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < matrix_size; i++) 
        {
            for (int k = 0; k < matrix_size; k++) 
            {
                myfile << to_string(later_state[i][k])  + " ";
            }
            myfile << "\n";
        }
    myfile.close();
    }
    
    
    /*Locating the transition*/
    /*
    double h = 0.0;
    double j = 1.0;
    int beta_size = 21;
    
    double beta[beta_size] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,
    5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10};
    
    double beta[beta_size] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
    0.55,0.6,0.65};
    
    int matrix_size = 100;
    //Averages and Standard dvs. 
    double values[beta_size][2];
    
    int **later_state;
    later_state= new int *[matrix_size];
    for(int i = 0; i < matrix_size; i++)
    {
        later_state[i] = new int[matrix_size];
    }
    for(int i = 0; i < matrix_size; i++)
    {
        for(int j = 0; j < matrix_size; j++)
        {
            later_state[i][j] = 0;
        }
    }
    
    for (int i = 0; i < beta_size; i++)
    {
        int **initial_state;
        initial_state = new int *[matrix_size];
        for(int i = 0; i < matrix_size; i++)
        {
            initial_state[i] = new int[matrix_size];
        }
        for(int i = 0; i < matrix_size; i++)
        {
            for(int j = 0; j < matrix_size; j++)
            {
                int ranchoice = (rand() > RAND_MAX/2) ? -1 : 1;
                initial_state[i][j] = ranchoice; 
            }
        }
        magnetization_density_average_sq(initial_state, 
        h, j, beta[i], matrix_size, values[i],later_state);
        
        susceptibility_average(initial_state, h, j, beta[i], matrix_size,
        values[i], later_state);
        
        binder_cumulant(initial_state, 
        h, j, beta[i], matrix_size, values[i],later_state);
    }
    
    ofstream myfile ("Transition_data.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < beta_size; i++) 
        {
            myfile << to_string(beta[i]) + " " + to_string(values[i][0]) +
            " " + to_string(values[i][1]);
            myfile << "\n";
        }
    myfile.close();
    }
    */
 
    return 0;
}

