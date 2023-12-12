// Name: Sahyadri Krishna
// Email: sk322@st-andrews.ac.uk

#include <iostream>
#include <fstream>
#include <math.h>
#include "FiniteFunctions.h"
#include <random>

int main(){

    // Load the Mystery Data
    std::string file_loc = "Outputs/data/MysteryData24111.txt";
    std::ifstream input_file(file_loc); // Input file object
    std::string line; // Holds each line of data
    std::vector<double> x; // The vector that holds input data

    if (!input_file.is_open()){
        std::cout << "Issue opening file!" << std::endl;
    }
    else{
        while (std::getline(input_file, line)){ // Access each available line without knowing size apriori.
            x.push_back(std::stod(line)); // Append the line at the back of x.
        }
    }

    // Create an instance of FiniteFunction
    std::vector<double> x2;
    double xmin = -5.0;
    double xmax = 5.0;
    std::string outName = "FiniteFunction"; // Name for the resultant plot

    FiniteFunction finite(xmin, xmax, outName);
    x2 = finite.metro_sample(10000); // Sample data from the inverse-square function.
    finite.plotData(x2, 40, false); // Plot the binned sampled data.
    finite.plotData(x, 40, true); // Plot the mystery data.
    finite.plotFunction(); // Plot the function.


    // Create an instance of NormalDistribution.
    // Save plot as a file "NormalDsitribution.png"
    double mean = 2.0;
    double sigma = 2.0;
    outName = "NormalDistribution";

    NormalDistribution normal(mean, sigma, xmin, xmax, outName); // Instance of NormalDistribution
    
    x2 = normal.metro_sample(10000); // Sampled data for same distribution
    normal.plotData(x2, 40, false); // Plot binned sample data
    normal.plotData(x, 40, true); // Plot binned mystery data.
    normal.plotFunction(); // Plot the normal distribution for the given mean and sigma.
    normal.printPars();


    // Create an instance of CauchyDistribution.
    // Save plot as file "CauchyDistribution.png"
    // TODO: Add a check to see if gamma > 0.
    double gamma = 1.5;
    double x_init = 2.0;
    outName = "CauchyDistribution";
    
    CauchyDistribution cDist(gamma, x_init, xmin, xmax, outName);

    x2 = cDist.metro_sample(10000); // Sampled data for same distribution
    cDist.plotData(x2, 40, false); // Plot binned sample data
    cDist.plotData(x, 40, true); // Plot binned mystery data
    cDist.plotFunction(); // Plot the Cauchy distribution
    cDist.printPars();


    // Create an instance of Negative Crystal Ball distribution
    double n = 2.0;
    double alpha = 1.0;
    mean = 2.0;
    sigma = 1.5;
    outName = "CrystalBallDistribution";
    
    CrystalBall cBall(mean, sigma, n, alpha, xmin, xmax, outName);

    x2 = cBall.metro_sample(10000); // Sampled data for same distribution
    cBall.plotData(x2, 40, false); // Plot the binned sampled data
    cBall.plotData(x, 40, true); // Plot the mystery data
    cBall.plotFunction(); // Plot the Negative Crystal Ball function.
    cBall.printPars();


    return 0;
}