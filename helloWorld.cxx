#include <iostream>
#include <math.h>

// Create a function to calculate the magnitude of a vector
double calcMag(double x, double y){
    double mag;
    mag = sqrt(x*x + y*y);
    return mag;
}


int main()
{
    std::cout << "Hello World!" << std::endl;

    double x = 2.3;
    double y = 4.5;
    double mag, mag2;

    // Want to calculate the magnitude of the vector
    mag = sqrt(x*x + y*y);
    std::cout << "Magnitude of the vectors is (non-function)= " << mag << std::endl;

    // Now use the code to get the magnitude
    mag2 = calcMag(x, y);
    std::cout << "Magnitued of the vectors is (function) = " << mag2 << std::endl;

    return 0;

}