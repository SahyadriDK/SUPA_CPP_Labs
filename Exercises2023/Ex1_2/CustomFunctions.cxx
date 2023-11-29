// Author: Sahyadri Krishna
// Email: sk322@st-andrews.ac.uk

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>


// Define a function to open a file
std::vector<std::string> read_file(std::string filename, std::vector<std::string>& x, std::vector<std::string>& y){

    std::ifstream myFile(filename); // Open the file at location "filename" using ifstream.
    std::vector<std::string> v; // Vector that will hold each individual line.
    int nLine = 0; // Line counter.

    // Iterate over each filename
    std::string line;
    std::cout << "Printing lines in file" << std::endl;
    std::cout << " " << std::endl;

    // Add each line to a vector called v.
    if (!myFile.is_open()){ // Check if no issues opening file
        std::cout << "Issue opening file!" << std::endl;
    }
    else{
        while (std::getline(myFile, line)){ // Access each available line without knowing size apriori.
            v.push_back(line); // Append the line at the back of v.
            nLine++; // Add 1 to nLine.
        }
    }

    // Split each line into individual x and y values. Add the x and y values into individual vectors.
    int pos;
    for (int i=0; i<v.size(); i++){
        if (i==0) continue;
        pos = v[i].rfind(","); // Find the exact index where the ',' delimiter occurs.
        x.push_back(v[i].substr(0, pos)); // Find that part of the string that is before ',' and add it to vector x.
        y.push_back(v[i].substr(pos+1, v[i].size())); // Find that part of the string that is after ',' and add it to y.
    }

    return v;
}


// Define a function to print the lines
void print_file(int n, std::vector<std::string> v){

    // Use a loop to print each line, if requested number of lines doesnt exceed those available in file.
    for (int i=0; i<n; i++){
        if (n > v.size()){ // Check to see if number of lines doesnt exceed those available.
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
        // NOTE: x and y values stored earlier as strings. Coverting to floats/doubles.
        // Add magnitude to the end of mag vector.
        mag.push_back(sqrt(std::stof(x[i]) * std::stof(x[i]) + std::stof(y[i])*std::stof(y[i])));
    }

    return mag;
}


// Want to calculate the best-fit line given a set of (x,y)
// Return: m (slope) and c (intercept)
std::string FitLine(std::vector<std::string>x, std::vector<std::string>y){
    
    int N = x.size(); // Length of the vectors
    double m, c; // Variables for slope and intercept respectively.s

    // Fit the line using Least Squares Method.
    double sum_xy, sum_x, sum_y, sum_xx;
    // Use the for loop to iteratively calculate the sum of each component of formula.
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

    // Use previously calculated sums in the formulas given in assignement sheet.
    m = ((N * sum_xy) - (sum_x * sum_y)) / ((N * sum_xx) - (sum_x * sum_x));
    c = ((sum_xx * sum_y) - (sum_xy * sum_x)) / ((N * sum_xx) - (sum_x * sum_x));

    // Read the error file
    std::string err_file = "Exercises2023/Ex1_2/error2D_float.txt"; // Load error file.
    std::vector<std::string> err_x, err_y; // variables that hold x and y error bars.
    std::vector<std::string> errLine; // Each line of the error file.
    errLine = read_file(err_file, err_x, err_y); // Load error bars.

    // Use the errors to caclulate the chisquared
    double chi_sq;
    int Ndf = err_x.size() - 2; // Degrees of freedom
    // Calculate chisquare in a loop
    for (int i=0; i<err_x.size(); i++){
        if (i==0){
            chi_sq = ((std::stod(y[i]) - (m*std::stod(x[i])+c)) *  (std::stod(y[i]) - (m*std::stod(x[i])+c))) / (std::stod(err_y[i]) * std::stod(err_y[i]));
        }
        else{
            chi_sq += (((std::stod(y[i]) - (m*std::stod(x[i])+c)) *  (std::stod(y[i]) - (m*std::stod(x[i])+c))) / (std::stod(err_y[i]) * std::stod(err_y[i])));
        }
    }
    chi_sq /= Ndf; // Division by DoF. Gives reduced chisquare.

    // Print out the the chi_sq
    std::string out_string;
    out_string = "The least square fit line is y = " + std::to_string(m) + "x + " + std::to_string(c) + ". The reduced chisquared for the fit is " + std::to_string(chi_sq);
    
    return out_string;

}


// Want to define a function to calculate exponent recursively.
double calcExp(double x, double y){

    // Round y to an integer
    int y_int;
    y_int = (int) y; // Forces double y to be an integer (typecast). This acts like a floor function.

    // Check if the y value is less than 1
    if (y < 1){
        return 1;
    }

    // Recursively calculate x^y.
    if (y_int==1){
        return x;
    }
    
    return x * calcExp(x, y-1); // Repeat multiplication. This will return actual exponent.
}


// The previous function calculates exponents for each individual element of a vector
// Use recursion to access each element of x and y vectors
int calcExpArr(std::vector<std::string> x, std::vector<std::string> y, std::vector<double>& exp_arr, int i){

    if (i==x.size()){
        return 0;
    }

    double val;
    val = calcExp(std::stod(x[i]), std::stod(y[i])); // For each pair of x and y, recursively calculate x^y.
    std::cout << val << std::endl;
    exp_arr.push_back(val); // Store value in array.
    
    return calcExpArr(x, y, exp_arr, i+1); // Repeat calculation till you reach the end.

}


// Create a function to store the outputs into a text file.
// Function 1: Input uses two vectors.
void save_file(std::vector<std::string> x, std::vector<std::string> y, std::string filename){
    
    std::ofstream myOutput; // Load ofstream variable.
    myOutput.open(filename); // Open said file.

    // Save output using loop to ofstream variable.
    for (int i=0; i<x.size(); i++){
        myOutput << x[i] << "," << y[i] << std::endl;
    }
    myOutput.close(); // Close writing of file

}


// Function 2: Input uses a single vector.
void save_file(std::vector<double> x, std::string filename){

    std::ofstream myOutput; // Load ofstream variable
    myOutput.open(filename); // Open said file

    for (int i=0; i<x.size(); i++){
        myOutput << x[i] << std::endl;
    }
    myOutput.close(); // Close writing of file

}


// Function 3: Specifically only for y=mx+c.
void save_file(std::string data_string, std::string filename){
    
    std::ofstream myOutput;
    myOutput.open(filename);

    myOutput << data_string << std::endl; // Simply save string to ofstream variable.
    myOutput.close();
}