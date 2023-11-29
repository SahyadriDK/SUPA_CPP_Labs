#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

std::vector<std::string> read_file(std::string, std::vector<std::string>&, std::vector<std::string>&);
void print_file(int, std::vector<std::string>);
std::vector<double> calculate_magnitude(std::vector<std::string>, std::vector<std::string>);
std::string FitLine(std::vector<std::string>, std::vector<std::string>);
double calcExp(double, double);
int calcExpArr(std::vector<std::string>, std::vector<std::string>, std::vector<double>&, int = 0);
void save_file(std::vector<std::string>, std::vector<std::string>, std::string);
void save_file(std::vector<double>, std::string);
void save_file(std::string, std::string);