#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "CustomFunctions.h"

int main(){

    int option, retry;
    bool go = true;

    // First, read the file
    std::string filename = "Exercises2023/Ex1_2/input2D_float.txt";
    std::vector<std::string> v;
    std::vector<std::string> x;
    std::vector<std::string> y;
    std::vector<double> mag;
    std::vector<double> expArr;
    double exp;
    std::string new_file;

    // Load the file
    v = read_file(filename, x, y);


    while (go){

        std::cout << "Which of the following functions would you like to use?" << std::endl;
        std::cout << "Print data (enter 1), Calculate Magnitudes (enter 2), Fit a line (enter 3), Calculate Exponents (enter 4)" << std::endl;
        std::cin >> option;

        switch(option){

            case 1: {
                // Print out data
                int n;
                std::cout << "" << std::endl;
                std::cout << "How many lines of data do you want printed?"  << std::endl;
                std::cin >> n;
                std::cout << "" << std::endl;
                print_file(n, v);
                break;

            }

            case 2: {
                // Calculate the magnitudes of the vector
                mag = calculate_magnitude(x, y);
                new_file = "Exercises2023/Ex1_2/magnitudes.txt";
                save_file(mag, new_file);
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            case 3: {
                // Print out the best fit line
                std::string fit;
                fit = FitLine(x,y);
                new_file = "Exercises2023/Ex1_2/best_fit.txt";
                save_file(fit, new_file);
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            case 4: {
                // Calculate x^y
                calcExpArr(x, y, expArr);
                new_file = "Exercises2023/Ex1_2/exponents.txt";
                save_file(expArr, new_file);
                std::cout << "" << std::endl;
                std::cout << "Output saved to file" << std::endl;
                break;
            }

            default:{
                std::cout << "" << std::endl;
                std::cout << "Wrong choice. Try again!" << std::endl;
                break;
            }
        }

        std::cout << "" << std::endl;
        std::cout << "Output complete! Would you like to use another function?" << std::endl;
        std::cout << "Press 1 to exit. Press anything else to continue" << std::endl;
        std::cin >> retry;

        if (retry == 1){
            std::cout << "" << std::endl;
            std::cout << "Exiting code!" << std::endl;
            go = false;
            break;
        }
    }

    return 0;

}