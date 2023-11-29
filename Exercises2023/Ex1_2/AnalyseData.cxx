// Author: Sahyadri Krishna
// Email: sk322@st-andrews.ac.uk

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "CustomFunctions.h"

int main(){

    int option, retry; // Variables that reflect user choice.
    bool go = true; // Set to false if code to be exited.

    // First, read the file
    std::string filename = "Exercises2023/Ex1_2/input2D_float.txt";
    std::vector<std::string> v; // Vector that holds individual read lines.
    std::vector<std::string> x; // Vector that holds x values from data file.
    std::vector<std::string> y; // Vector that holds y values from data file.
    std::vector<double> mag; // Vector that holds calculated magnitude values.
    std::vector<double> expArr; // Vector that holds calculated exponent values.
    double exp; // Holds exponent value temporarily.
    std::string new_file; // Variable that stores name of newly written output file.

    // Load the file
    v = read_file(filename, x, y);

    // While loop that allows multiple user inputs.ß
    while (go){

        std::cout << "Which of the following functions would you like to use?" << std::endl;
        std::cout << "Print data (enter 1), Calculate Magnitudes (enter 2), Fit a line (enter 3), Calculate Exponents (enter 4)" << std::endl;
        std::cin >> option;

        switch(option){

            case 1: {
                // Print out dataß
                int n;
                std::cout << "" << std::endl;
                std::cout << "How many lines of data do you want printed?"  << std::endl;
                std::cin >> n;
                std::cout << "" << std::endl;
                print_file(n, v); // Prints n lines of data vector v.ß
                break;

            }

            case 2: {
                // Calculate the magnitudes of the vector
                mag = calculate_magnitude(x, y); // Calculates magnitude vector from x and y vectors.ßß
                new_file = "Exercises2023/Ex1_2/magnitudes.txt";
                save_file(mag, new_file); // Saves magnitude vector to file.
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            case 3: {
                // Print out the best fit line
                std::string fit;
                fit = FitLine(x,y); // Fits line via least squares fitting and calculates reduced chisquare.
                new_file = "Exercises2023/Ex1_2/best_fit.txt";
                save_file(fit, new_file); // Saves information to text file.
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            case 4: {
                // Calculate x^y
                calcExpArr(x, y, expArr); // Recursively calculates the value of x^y for each x and y.
                new_file = "Exercises2023/Ex1_2/exponents.txt";
                save_file(expArr, new_file); // Saves exponent values to text file.
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            default:{
                // Offers alternative chance if wrong option presented
                std::cout << "" << std::endl;
                std::cout << "Wrong choice. Try again!" << std::endl;
                break;
            }
        }

        std::cout << "" << std::endl;
        std::cout << "Output complete! Would you like to use another function?" << std::endl;
        std::cout << "Press 1 to exit. Press anything else to continue" << std::endl;
        std::cin >> retry;

        // Quits code if retry set to 1. Else continues executing code under the while loop.
        if (retry == 1){
            std::cout << "" << std::endl;
            std::cout << "Exiting code!" << std::endl;
            go = false;
            break;
        }
    }

    return 0;

}