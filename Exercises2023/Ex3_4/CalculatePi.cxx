#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <random>


double calculate_pi(double radius, int N_random){

    // Throw darts using the random pdfs in x and y directions
    double x, y; // The positions of the randomly thrown darts.
    unsigned int seed = 0;
    double pi; // The value of pi to be calculated.
    int count = 0; // The number of dart throws that land in the circle
    std::mt19937 mtEngine{seed};
    std::uniform_real_distribution<double> uniformPDF_x{0.0, 2.0*radius};
    std::uniform_real_distribution<double> uniformPDF_y{0.0, 2.0*radius};

    for (int i=0; i<N_random; i++){
        x = uniformPDF_x(mtEngine);
        y = uniformPDF_y(mtEngine);

        //std::cout << "(x,y) value " << x << " " << y << std::endl;

        if (sqrt(pow(x-radius, 2) + pow(y-radius, 2)) < radius){
            count += 1;
        }
    }

    pi = 4.0 * count /  N_random;
    return pi;
}


int main(int argc, char *argv[]){

    double pi;
    pi = calculate_pi(std::stod(argv[1]), std::stod(argv[2]));

    std::cout << "Value of PI is - " << pi << std::endl;

}