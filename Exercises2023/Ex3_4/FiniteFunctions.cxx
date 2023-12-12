#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <random>

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  m_OutName = outfile;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  // Use Trapezoid rule to calculate the integral.
  double dx = (m_RMax - m_RMin) / (double) Ndiv;
  double start = m_RMin;
  for (int i=0; i<Ndiv; i++){
    if (i==0){
      m_Integral = 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    else{
      m_Integral += 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    start += dx;
  }
  return m_Integral;  
}

double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
#################################
// Metropolis-Hastings sampling
#################################
*/
std::vector<double> FiniteFunction::metro_sample(int n_random){
  std::vector<double> sampled_data; // The output vector of sampled data points
  double f_x, f_y, x_init, y, A, T; // Various parameters used in algorithm

  // Use the Metropolis-Hastings algorithm to generate random numbers.
  // First, initialise random number generator
  unsigned int seed = 42;
  std::mt19937 mtEngine{seed};
  std::uniform_real_distribution<double> uniformPDF{m_RMin, m_RMax};
  std::uniform_real_distribution<double> uniformPDF2{0.0, 1.0};
  // Create an initial value of x, called x_init
  x_init = uniformPDF(mtEngine);

  // Loop over n_random and generate random numbers using the algorithm
  for (int i=0; i<n_random; i++){
    std::normal_distribution<double> normalPDF{x_init, 1.5}; // Sigma arbitrary.
    y = normalPDF(mtEngine);
    f_y = this->callFunction(y); // PROBLEM IS HERE
    f_x = this->callFunction(x_init); // PROBLEM IS HERE
    if (f_y/f_x < 1){
      A = f_y/f_x;
    }
    else{
      A = 1;
    }
    T = uniformPDF2(mtEngine);
    if (T < A){
      x_init = y;
    }

    // Append the newly generated random numbers to vector
    sampled_data.push_back(x_init);
    
  }

  return sampled_data;
  //this->callFunction(1.0);
}



/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}




/*
###########################################################################
##### Creating class contents for child class "NormalDistribution" ########
###########################################################################
*/

NormalDistribution::NormalDistribution():FiniteFunction(){
  m_Mean = 0.0; // Default mean attribute
  m_Sigma = 1.0; // Default sigma attribute
}

NormalDistribution::NormalDistribution(double mean, double sigma):FiniteFunction(){
  m_Mean = mean;
  m_Sigma = sigma;
}

NormalDistribution::NormalDistribution(double mean, double sigma, double range_min, double range_max, std::string outfile):FiniteFunction(range_min, range_max, outfile){
  m_Mean = mean;
  m_Sigma = sigma;
}

double NormalDistribution::gauss(double x, double mean, double sigma) {return (1.0/(sqrt(2*M_PI)*sigma)) * exp(-pow((x-mean), 2)/pow(sigma, 2));};
double NormalDistribution::callFunction(double x) {return this->gauss(x, m_Mean, m_Sigma);};


// Modified "Integral" to include m_Sigma and m_Mean
double NormalDistribution::integrate(int Ndiv){ //private
  // Use Trapezoid rule to calculate the integral.
  double dx = (m_RMax - m_RMin) / (double) Ndiv;
  double start = m_RMin;
  for (int i=0; i<Ndiv; i++){
    if (i==0){
      m_Integral = 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    else{
      m_Integral += 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    start += dx;
  }
  return m_Integral;  
}

double NormalDistribution::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

//Modified "Scan" to have callFunction include m_Sigma and m_Mean
std::vector< std::pair<double,double> > NormalDistribution::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

// Define function to print out parameters (TEMP)
void NormalDistribution::printPars(){
  std::cout << "Mean is " << m_Mean << std::endl;
  std::cout << "Sigma is " << m_Sigma << std::endl;
  std::cout << "Outfile name: " << m_OutName << std::endl;

}



/*
###########################################################################
##### Creating class contents for child class "CauchyDistribution" ########
###########################################################################
*/

CauchyDistribution::CauchyDistribution():FiniteFunction(){
  m_Gamma = 1.0; // Default mean attribute
  m_Xinit = 0.0; // Default sigma attribute
}

CauchyDistribution::CauchyDistribution(double gamma, double x_init):FiniteFunction(){
  m_Gamma = gamma;
  m_Xinit = x_init;
}

CauchyDistribution::CauchyDistribution(double gamma, double x_init, double range_min, double range_max, std::string outfile):FiniteFunction(range_min, range_max, outfile){
  m_Gamma = gamma;
  m_Xinit = x_init;
}

double CauchyDistribution::cauchy(double x, double gamma, double x_init) {return (1/(M_PI*gamma) / (1 + pow((x-x_init)/gamma, 2)));};
double CauchyDistribution::callFunction(double x) {return this->cauchy(x, m_Gamma, m_Xinit);};


// Modified "Integral" to include m_Sigma and m_Mean
double CauchyDistribution::integrate(int Ndiv){ //private
  // Use Trapezoid rule to calculate the integral.
  double dx = (m_RMax - m_RMin) / (double) Ndiv;
  double start = m_RMin;
  for (int i=0; i<Ndiv; i++){
    if (i==0){
      m_Integral = 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    else{
      m_Integral += 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    start += dx;
  }
  return m_Integral;  
}

double CauchyDistribution::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

//Modified "Scan" to have callFunction include m_Sigma and m_Mean
std::vector< std::pair<double,double> > CauchyDistribution::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

// Define function to print out parameters (TEMP)
void CauchyDistribution::printPars(){
  std::cout << "Gamma is " << m_Gamma << std::endl;
  std::cout << "x0 is " << m_Xinit << std::endl;
  std::cout << "Outfile name: " << m_OutName << std::endl;
}




/*
###########################################################################
##### Creating class contents for child class "CrystalBall" ########
###########################################################################
*/

CrystalBall::CrystalBall():FiniteFunction(){
  m_Mean = 0.0; // Default mean attribute
  m_Sigma = 1.0; // Default sigma attribute
  m_n = 2.0; // Default n (>1)
  m_Alpha = 1.0; // Default alpha (>0)
}

CrystalBall::CrystalBall(double mean, double sigma, double n, double alpha):FiniteFunction(){
  m_Mean = mean;
  m_Sigma = sigma;
  m_n = n;
  m_Alpha = alpha;
}

CrystalBall::CrystalBall(double mean, double sigma, double n, double alpha, double range_min, double range_max, std::string outfile):FiniteFunction(range_min, range_max, outfile){
  m_Mean = mean;
  m_Sigma = sigma;
  m_n = n;
  m_Alpha = alpha;
}

double CrystalBall::negCrysBall(double x, double mean, double sigma, double n, double alpha) {

  double A, B, C, D, N;
  A = pow(n/alpha, n) * exp(-pow(alpha,2)/2);
  B = (n/alpha) - alpha;
  C = ((n/alpha)/(n-1)) * exp(-pow(alpha, 2)/2);
  D = sqrt(0.5*M_PI) * (1 + erf(alpha/sqrt(2)));
  N = 1 / (sigma * (C + D));

  if ((x-mean)/sigma > -alpha){
    return exp(-pow(x-mean, 2.0)/(2*pow(sigma, 2.0)));
  }
  else if ((x-mean)/sigma <= -alpha){
    return A * pow(B - ((x-mean)/sigma), -n);
  }
};

double CrystalBall::callFunction(double x) {
  return this->negCrysBall(x, m_Mean, m_Sigma, m_n, m_Alpha);
};


// Modified "Integral" to include m_Sigma and m_Mean
double CrystalBall::integrate(int Ndiv){ //private
  // Use Trapezoid rule to calculate the integral.
  double dx = (m_RMax - m_RMin) / (double) Ndiv;
  double start = m_RMin;
  for (int i=0; i<Ndiv; i++){
    if (i==0){
      m_Integral = 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    else{
      m_Integral += 0.5 * (this->callFunction(start+dx) + this->callFunction(start)) * dx;
    }
    start += dx;
  }
  return m_Integral;  
}

double CrystalBall::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}

//Modified "Scan" to have callFunction include m_Sigma and m_Mean
std::vector< std::pair<double,double> > CrystalBall::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

// Define function to print out parameters (TEMP)
void CrystalBall::printPars(){
  std::cout << "Mean is " << m_Mean << std::endl;
  std::cout << "Sigma is " << m_Sigma << std::endl;
  std::cout << "n is " << m_n << std::endl;
  std::cout << "Alpha is " << m_Alpha << std::endl;
  std::cout << "Outfile name: " << m_OutName << std::endl;
}