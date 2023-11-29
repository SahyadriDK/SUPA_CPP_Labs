#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

std::vector<std::string> read_file(std::string, std::vector<std::string>&, std::vector<std::string>&);
void print_file(int, std::vector<std::string>);
std::vector<double> calculate_magnitude(std::vector<std::string>, std::vector<std::string>);
void FitLine(std::vector<std::string>, std::vector<std::string>);
double calcExp(double, double);
int calcExpArr(std::vector<std::string>, std::vector<std::string>, std::vector<double>&, int);