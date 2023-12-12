#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <random>


double calculate_pi(double radius, int N_random){

    // Throw darts using the random pdfs in x and y directions
    double x, y; // The positions of the randomly thrown darts.
    unsigned int seed = 1;
    double pi; // The value of pi to be calculated.
    int count = 0; // The number of dart throws that land in the circle
    std::mt19937 mtEngine{seed};
    std::uniform_real_distribution<double> uniformPDF_x{0.0, 2.0*radius};
    std::uniform_real_distribution<double> uniformPDF_y{0.0, 2.0*radius};

    for (int i=0; i<N_random; i++){
        // Randomly generate (x,y) coordinates of each dart throw.
        x = uniformPDF_x(mtEngine); // Randomly generate x value inside box.
        y = uniformPDF_y(mtEngine); // Randomly generate y value inside box.

        // Check (x,y) pair to see if they lie within "radius" of circle.
        if (sqrt(pow(x-radius, 2) + pow(y-radius, 2)) < radius){ // If within radius
            count += 1; // Count the dart as a strike within circle
        }
    }

    // Probability that the darts land in circle is simply the fraction of area occupied by the circle in the embedding square.
    // For any radius r, the box side is 2r.
    // This means the probability is PI/4.
    // Thus PI is 4 x the probability of darts landing on circle.
    pi = 4.0 * ((float) count /  (float) N_random);
    return pi;
}


int main(int argc, char *argv[]){

    double pi;
    pi = calculate_pi(std::stod(argv[1]), std::stod(argv[2]));

    // Use the printf command with .10f to display a floating point value of pi with 10 decimal places.
    std::printf("The value of PI is - %.10f \n", pi);

}