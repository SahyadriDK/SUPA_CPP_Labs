#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  virtual double integral(int Ndiv = 1000); 
  virtual std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  std::vector<double> metro_sample(int n_random);
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  std::string m_OutName;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  virtual double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  
private:
  double invxsquared(double x); //The default functional form
};


class NormalDistribution : public FiniteFunction{

public:
  NormalDistribution();
  NormalDistribution(double mean, double sigma);
  NormalDistribution(double mean, double sigma, double range_min, double range_max, std::string outfile);
  double callFunction(double x);
  std::vector< std::pair<double,double> > scanFunction(int Ndiv);
  double integral(int Ndiv=1000);
  void printPars();

protected:
  double m_Mean;
  double m_Sigma;
  double integrate(int Ndiv); 


private:
  double gauss(double x, double mean, double sigma);
};


class CauchyDistribution : public FiniteFunction{

public:
  CauchyDistribution(); // assume gamma=1 and x0 = 0
  CauchyDistribution(double gamma, double x_init);
  CauchyDistribution(double gamma, double x_init, double xmin, double xmax, std::string outfile);
  double callFunction(double x);
  std::vector< std::pair<double,double> > scanFunction(int Ndiv);
  double integral(int Ndiv=1000);
  void printPars();

protected:
  double m_Gamma;
  double m_Xinit;
  double integrate(int Ndiv);

private:
  double cauchy(double x, double gamma, double x_init);
};


class CrystalBall : public FiniteFunction{

public:
  CrystalBall(); // assume gamma=1 and x0 = 0
  CrystalBall(double mean, double sigma, double n, double alpha);
  CrystalBall(double mean, double sigma, double n, double alpha, double xmin, double xmax, std::string outfile);
  double callFunction(double x);
  std::vector< std::pair<double,double> > scanFunction(int Ndiv);
  double integral(int Ndiv=1000);
  void printPars();

protected:
  double m_Mean;
  double m_Sigma;
  double m_n;
  double m_Alpha;
  double integrate(int Ndiv);

private:
  double negCrysBall(double x, double mean, double sigma, double n, double alpha);


};