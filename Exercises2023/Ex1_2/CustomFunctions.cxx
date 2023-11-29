#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

// Define a function to open a file
std::vector<std::string> read_file(std::string filename, std::vector<std::string>& x, std::vector<std::string>& y){

    std::ifstream myFile(filename);
    std::vector<std::string> v;
    int nLine = 0;

    // Iterate over each filename
    std::string line;
    std::cout << "Printing lines in file" << std::endl;
    std::cout << " " << std::endl;

    // Add each line to a vector called v.
    if (!myFile.is_open()){
        std::cout << "Issue opening file!" << std::endl;
    }
    else{
        while (std::getline(myFile, line)){
            v.push_back(line);
            nLine++;
        }
    }

    // Split each line into individual x and y values. Add the x and y values into individual vectors.
    int pos;
    for (int i=0; i<v.size(); i++){
        if (i==0) continue;
        pos = v[i].rfind(",");
        x.push_back(v[i].substr(0, pos));
        y.push_back(v[i].substr(pos+1, v[i].size()));
    }

    return v;
}


// Define a function to print the lines
void print_file(int n, std::vector<std::string> v){

    // Use a loop to print each line, if requested number of lines doesnt exceed those available in file.
    for (int i=0; i<n; i++){
        if (n > v.size()){
            std::cout << "The number of lines to be printed is larger than the number of lines in the file. Exiting!" << std::endl;
            break;
        }
        else{
            std::cout << v[i] << std::endl;
        }
    }
}


// Define a function to calculate magnitudes given x and y vectors
std::vector<double> calculate_magnitude(std::vector<std::string> x, std::vector<std::string> y){

    std::vector<double> mag; // The vector that will hold magnitudes.
    
    // Loop over all items in x and y and calculate the individual magnitudes. Store them in vector.
    for (int i=0; i<x.size(); i++){
        mag.push_back(sqrt(std::stof(x[i]) * std::stof(x[i]) + std::stof(y[i])*std::stof(y[i])));
    }

    return mag;
}

// Want to calculate the best-fit line given a set of (x,y)
// Return: m (slope) and c (intercept)
void FitLine(std::vector<std::string>x, std::vector<std::string>y){
    
    int N = x.size(); // Length of the vectors
    double m, c;

    // Fit the line using Least Squares Method.
    double sum_xy, sum_x, sum_y, sum_xx;
    for (int i=0; i<x.size(); i++){
        if (i==0){
            sum_xy = std::stod(x[i]) * std::stod(y[i]);
            sum_xx = std::stod(x[i]) * std::stod(x[i]);
            sum_x = std::stod(x[i]);
            sum_y = std::stod(y[i]);
        }
        else{
            sum_xy += (std::stod(x[i]) * std::stod(y[i]));
            sum_xx += (std::stod(x[i]) * std::stod(x[i]));
            sum_x += std::stod(x[i]);
            sum_y += std::stod(y[i]);
        }
    }

    m = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_xx) - (sum_x * sum_x));
    c = ((sum_xx * sum_y) - (sum_xy * sum_x)) / ((N * sum_xx) - (sum_x * sum_x));

    // Print out the equation
    std::cout << "The least square fit line is y = " << m << "x + " << c << std::endl;

    // Read the error file
    std::string err_file = "Exercises2023/Ex1_2/error2D_float.txt";
    std::vector<std::string> err_x, err_y;
    std::vector<std::string> errLine;
    errLine = read_file(err_file, err_x, err_y);

    // Use the errors to caclulate the chisquared
    double chi_sq;
    int Ndf = err_x.size() - 2; // Degrees of freedom
    for (int i=0; i<err_x.size(); i++){
        if (i==0){
            chi_sq = ((std::stod(y[i]) - (m*std::stod(x[i])+c)) *  (std::stod(y[i]) - (m*std::stod(x[i])+c))) / (std::stod(err_y[i]) * std::stod(err_y[i]));
        }
        else{
            chi_sq += (((std::stod(y[i]) - (m*std::stod(x[i])+c)) *  (std::stod(y[i]) - (m*std::stod(x[i])+c))) / (std::stod(err_y[i]) * std::stod(err_y[i])));
        }
    }
    chi_sq /= Ndf; // Division by DoF.

    // Print out the the chi_sq
    std::cout << "The reduced chisquared for the fit is " << chi_sq << std::endl;


}


// Want to define a function to calculate exponent recursively.
double calcExp(double x, double y){

    // Round y to an integer
    int y_int;
    y_int = (int) y;

    // Check if the y value is less than 1
    if (y < 1){
        return 1;
    }


    // Recursively calculate x^y.
    if (y_int==1){
        return x;
    }
    
    return x * calcExp(x, y-1); // This will return actual exponent
}

// The previous function calculates exponents for each individual element of a vector
// Use recursion to access each element of x and y vectors
int calcExpArr(std::vector<std::string> x, std::vector<std::string> y, std::vector<double>& exp_arr, int i){

    if (i==x.size()){
        return 0;
    }

    double val;
    val = calcExp(std::stod(x[i]), std::stod(y[i]));
    std::cout << val << std::endl;
    exp_arr.push_back(val); // Store value in array
    
    return calcExpArr(x, y, exp_arr, i+1); // Repeat calculation till you reach the end

}