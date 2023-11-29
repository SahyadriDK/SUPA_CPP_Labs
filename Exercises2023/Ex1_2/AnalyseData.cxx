#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "CustomFunctions.h"

int main(){
    
    // First, read the file
    std::string filename = "Exercises2023/Ex1_2/input2D_float.txt";
    std::vector<std::string> v;
    std::vector<std::string> x;
    std::vector<std::string> y;
    std::vector<double> mag;
    std::vector<double> expArr;
    double exp;

    // Load the file
    v = read_file(filename, x, y);

    // Print the file
    int n = 3; // No. of lines to print
    //print_file(n, v);

    // Calculate the magnitudes of the vector
    mag = calculate_magnitude(x, y);
    for (int i=0; i<mag.size(); i++){
        std::cout << mag[i] << std::endl;
    }

    // Print out the best fit line
    FitLine(x, y);

    // Print out x^y
    calcExpArr(x, y, expArr, 0);
    for (int i=0; i<expArr.size(); i++){
        std::cout << expArr[i] << std::endl;
    }

    return 0;

}