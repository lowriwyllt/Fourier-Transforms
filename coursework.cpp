// Include the libaries needed
#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include <string>
// Define constants etc - T is time period
using namespace std;
const double pi = acos(-1);
const double T = 2 * pi;
ofstream myfile_out;
ifstream myfile_in;
// Make a complex class for real and imaginary parts and member functions to work out the 
modulus and print complex numbers
class complex
{
public:
 double Re;
 double Im;
 complex(double, double, char);
 double modDbl();
 void printcomplex();
};
// Define the real and imaginary parts of a complex number (if only one number given assumes it is 
real part)
// polar : alpha =r, beta = theta, id 'p' where r * exp(i *theta)
// rectangular : alpha =Re, beta = Im, id 'r' where Re + Im *i
complex::complex(double alphaDbl, double betaDbl = 0.0, char idCh = 'r')
{
 if (idCh == 'p') // For polar complex numbers
 {
 Re = alphaDbl * cos(betaDbl);
 Im = alphaDbl * sin(betaDbl);
 }
 else if (idCh == 'r') // For rectangular complex numbers
 {
 Re = alphaDbl;
 Im = betaDbl;
 }
 else // Error message
 {
 cout << "constructor arguments are (r, theta, \'p\' or Re, Im, \'r\')" << endl;
 }
}
//Define member function to find the modulus of a complex number
double complex::modDbl()
{
 return sqrt(
 pow(Re, 2) +
 pow(Im, 2)
 );
}
//Define function modvector so the mod can work on a whole vector of complex numbers 
//Where F is the vector of complex which the mod of needs to be found
vector<double> modVecDbl(vector<complex> FVecCmplx)
{
 vector <double> modVecDbl = {};
 // Applies mod function to all elements in a vector
 for (int iInt = 0; iInt < FVecCmplx.size(); iInt++) 
 {
 modVecDbl.push_back(
 FVecCmplx[iInt].modDbl()
 );
 }
 return modVecDbl;
}
// Define member function to print a complex number
void complex::printcomplex()
{
 printf("%.2f + %.2f i", Re, Im);
}
//Define function printvectorcomplex so the printcomplex can work on a whole vector of complex 
numbers 
//where omega is a list of frequencies (can be list of time) and F is the vector of complex numbers to 
be printed 
void printvectorcomplex(vector<double> omegaVecDbl, vector<complex> FVecCmplx)
{
 cout << "n , w , Re + Im i , Mod \n";
 // Print the index of an element, the frequency and apply printcomplex function to an element 
and repeat over every element in a vector
 for (int iInt = 0; iInt < FVecCmplx.size(); iInt++) 
 {
 printf("%d , %.2f , ", iInt, omegaVecDbl[iInt]);
 FVecCmplx[iInt].printcomplex();
 printf(" , %.2f \n", FVecCmplx[iInt].modDbl());
 }
}
// Define a function add which adds two complex numbers
// Where add_one and add_two are both complex numbers
complex addCmplx(complex add_oneCmplx, complex add_twoCmplx)
{
 return complex(
 add_oneCmplx.Re + add_twoCmplx.Re,
 add_oneCmplx.Im + add_twoCmplx.Im,
 'r'
 );
}
// Define a function multiply which multiplies two complex numbers
// Where multi_one and multi_two are both complex numbers
complex multiplyCmplx(complex multi_oneCmplx, complex multi_twoCmplx)
{
 return complex(
 (multi_oneCmplx.Re * multi_twoCmplx.Re) - (multi_oneCmplx.Im * multi_twoCmplx.Im),
 (multi_oneCmplx.Im * multi_twoCmplx.Re) + (multi_oneCmplx.Re * multi_twoCmplx.Im),
 'r'
 );
}
// Define the function h1 as a complex number
// Where t is the time
complex h1Cmplx(double tDbl)
{
 return addCmplx(
 complex(1, tDbl, 'p'),
 complex(1, 5 * tDbl, 'p')
 );
}
// Define the function h2 as a complex number
//Where t is the time
complex h2Cmplx(double tDbl)
{
 return complex(exp(pow((tDbl - pi), 2) / 2));
}
//Define tuple which returns the vector of sample of a function and vector of double sample times
// Where f is the function and N is the number of samples needed
tuple<vector<double>, vector<complex>> sampleTup(complex(*fCmplx)(double tDbl), int NInt)
{
 vector<complex> fsampleVecCmplx;
 vector<double> tsampleVecDbl;
 double NDbl = (double)NInt;
 // Put sample value of of function f into fsample and the time (k*T/N) to tsample
 for (double kDbl = 0; kDbl < NInt; kDbl++) 
 {
 fsampleVecCmplx.push_back(
 fCmplx(kDbl * (T / NDbl))
 );
 tsampleVecDbl.push_back(
 kDbl * (T / NDbl)
 );
 }
 return make_tuple(tsampleVecDbl, fsampleVecCmplx);
}
// Prints time, complex number and mod from vector to file
//Where name is the name of the file, t is the time and f is the complex vector
void printfile(string filenameStr, vector<double> tVecDbl, vector<complex> fVecCmplx)
{
 myfile_out.open(filenameStr);
 myfile_out << "t" << "," << "Re" << "," << "Im" << "," << "Mod" << "\n";
 // for loop prints time, real, imaginary and modulus of complex number for every elemenet in 
cmplx vector
 for (size_t i = 0; i < fVecCmplx.size(); i++) 
 {
 myfile_out
 << tVecDbl[i] << ","
 << fVecCmplx[i].Re << ","
 << fVecCmplx[i].Im << ","
 << fVecCmplx[i].modDbl() << "\n";
 }
 myfile_out.close();
}
// Prints sample of a function to file
// Where name is the name of the file, f is the function and N is the number of samples needed
void printfile(string filenameStr, complex(*fCmplx)(double tDbl), int NInt)
{
 vector<complex> fsampleVecCmlx;
 vector<double> tsampleVecDbl;
 tie(tsampleVecDbl, fsampleVecCmlx) = sampleTup(fCmplx, NInt);
 printfile(filenameStr, tsampleVecDbl, fsampleVecCmlx);
}
//Define discrete Fourier transform using values
// Where h is vector of complex numbers to fourier transform and N is the number of samples
tuple<vector<double>, vector<complex>> DFTTup(vector <complex> hVecCmplx)
{
 vector<complex> HVecCmplx = {};
 vector<double> omegaVecDbl = {};
 double sizeDbl = hVecCmplx.size();
 // for loop sets the running total to zero then puts the running total into the Hcomplex vector and 
a value of omega into omega vector
 for (double nDbl = 0.0; nDbl < hVecCmplx.size(); nDbl++)
 {
 complex running_totalCmplx(0, 0, 'r');
 // for loop adds a multiplication of a element of a vector h[k] and exp(-2*i*pi*n*k/N) to a 
runnning total for every value of k from 0 to N-1
 for (double kDbl = 0.0; kDbl < hVecCmplx.size(); kDbl++)
 {
 running_totalCmplx = addCmplx(
 running_totalCmplx,
 multiplyCmplx(
 hVecCmplx[kDbl],
 complex(1, (-2. * pi * nDbl * kDbl / sizeDbl),
 'p'))
 );
 }
 HVecCmplx.push_back(running_totalCmplx);
 omegaVecDbl.push_back(nDbl * 2. * pi / sizeDbl);
 }
 return make_tuple(omegaVecDbl, HVecCmplx);
}
//Define tuple for discrete Fourier transform using a function returning double vector omega and 
vector complex H
//Where f is a function and N is the number of samples
tuple<vector<double>, vector<complex>> DFTTup(complex(*fCmplx)(double tDbl), int NInt)
{
 vector<complex> fsampleVecCmplx;
 vector<double> tsampleVecDbl;
 vector<complex> HVecCmplx;
 vector<double> omegaVecDbl;
 tie(tsampleVecDbl, fsampleVecCmplx) = sampleTup(fCmplx, NInt);
 tie(omegaVecDbl, HVecCmplx) = DFTTup(fsampleVecCmplx);
 return make_tuple(omegaVecDbl, HVecCmplx);
}
//Define inverse Fourier transform using values
// Where H is a complex vector and N is the number of samples used
tuple<vector<double>, vector<complex>>IFTTup(vector<complex> HVecCmplx)
{
 vector<complex> h_VecCmplx = {};
 vector<double> timeVecDbl = {};
 double sizeDbl = HVecCmplx.size();
 // for loop which sets the running total to zero and put the value of the running_total * (1/N) into 
hcomplex vector and value of time into time vector
 for (double kDbl = 0.0; kDbl < HVecCmplx.size(); kDbl++)
 {
 complex running_totalCmplx(0, 0, 'r');
 // for loop which goes through every element of n < N, only skipping n==s, and adds a
multiplication of H[n] element with exp(2*pi*n*k/N) to running total
 for (double nDbl = 0.0; nDbl < HVecCmplx.size(); nDbl++)
 {
 running_totalCmplx = addCmplx(
 running_totalCmplx,
 multiplyCmplx(
 HVecCmplx[nDbl],
 complex(
 1,
 (2 * pi * nDbl * kDbl / sizeDbl),
 'p')
 )
 );
 }
 h_VecCmplx.push_back(multiplyCmplx(
 complex(1 / sizeDbl),
 running_totalCmplx)
 );
 timeVecDbl.push_back(kDbl * T / sizeDbl);
 }
 return make_tuple(timeVecDbl, h_VecCmplx);
}
//Define a function which sets the value of the an element in a vector of complex numbers to 0
// where H is a vector of complex and element is the integer to be set to zero
vector<complex> setzeroVecCmplx(vector<complex> HVecCmplx, int elementInt)
{
 HVecCmplx[elementInt].Re = 0.0;
 HVecCmplx[elementInt].Re = 0.0;
 return HVecCmplx;
}
//Define a tuple which takes in a file and gives out n, time and complex values
// Where name is the name of the file where to get the values from
tuple<vector<int>, vector<double>, vector<complex>> getinfoTup(string nameStr)
{
 myfile_in.open(nameStr);
 vector<int> nFVecInt;
 vector<double> timeFVecDbl;
 vector<complex> FVecCmplx;
 // if statment checks the file is open
 if (myfile_in.is_open())
 {
 //while statments carries on if we havent reached the end of the file
 // while makes a string for each value needed, uses getline to get each value, assigns each value 
to the correct type and puts them into the correct vector
 while (!myfile_in.eof())
 {
 string nStr;
 string timeStr;
 string realStr;
 string imaginaryStr;
 getline(myfile_in, nStr, ',');
 getline(myfile_in, timeStr, ',');
 getline(myfile_in, realStr, ',');
 getline(myfile_in, imaginaryStr, '\n');
 //if nString is empty then we are reading an empty line, skip it
 if (nStr.empty())
 {
 continue;
 }
 int nInt = stoi(nStr);
 double timeDbl = stod(timeStr);
 double realFDbl = stod(realStr);
 double imaginaryFDbl = stod(imaginaryStr);
 complex FCmplx(realFDbl, imaginaryFDbl, 'r');
 nFVecInt.push_back(nInt);
 timeFVecDbl.push_back(timeDbl);
 FVecCmplx.push_back(FCmplx);
 }
 myfile_in.close();
 }
 return make_tuple(nFVecInt, timeFVecDbl, FVecCmplx);
}
// Define function finds greatest values and returns a vector only with the greatest values
// where nInt is how many greatest values is needed and F is a vector of complex number to find the 
greatest in where j is the integer
vector<complex> maxFVecCmplx(int nInt, vector<complex> FVecCmplx)
{
 vector<complex> maxF(FVecCmplx.size(), complex(0, 0, 'r'));
 vector<double> FmodVecDbl = modVecDbl(FVecCmplx);
 int maxjInt;
 // for loop sets maxn to zero then changes the value of the tempmod[maxn] to zero, so we don't 
count this again, then puts the complex value into maxF vector
 for (int i = 0; i < nInt; i++) {
 maxjInt = 0;
 // for loop checks every element in a vector of modulus to find the greatest and its index
 for (int jInt = 1; jInt < FmodVecDbl.size(); jInt++) {
 // if the maxn element is greater or equal to the n element don't do anyhting maxn is greatest
 if (FmodVecDbl[maxjInt] >= FmodVecDbl[jInt]) {
 }
 // otherwise (if maxn element is less than element n) then change maxn to be n this is the 
greatest value
 else {
 maxjInt = jInt;
 }
 }
 FmodVecDbl[maxjInt] = 0.0;
 maxF[maxjInt] = FVecCmplx[maxjInt];
 }
 return maxF;
}
int main() 
{
 // Q3 b - write h1 and h2 to file
 //N is how many samples
 int NInt = 100;
 // Printing h1 and h2 to a file
 printfile("h1sample.csv", &h1Cmplx, NInt);
 printfile("h2sample.csv", &h2Cmplx, NInt);
 // Q3 e - print H1 and H2
 // Define H1 and H2 and frequencies (omega)
 vector<double> omega_H1VecDbl; vector<complex> H1VecCmplx;
 vector<double> omega_H2VecDbl; vector<complex> H2VecCmplx;
 tie(omega_H1VecDbl, H1VecCmplx) = DFTTup(&h1Cmplx, NInt);
 tie(omega_H2VecDbl, H2VecCmplx) = DFTTup(&h2Cmplx, NInt);
 //printing H1 and H2
 cout << "H1 values \n";
 printvectorcomplex(omega_H1VecDbl, H1VecCmplx);
 cout << "\nH2 values \n";
 printvectorcomplex(omega_H2VecDbl, H2VecCmplx);
 // Q3 f - inverse fourier transform to get h'1 and h'2 (h_1 and h_1)
 vector<double> time_h_1VecDbl; vector<complex> h_1VecCmplx;
 vector<double> time_h_2VecDbl; vector<complex> h_2VecCmplx;
 tie(time_h_1VecDbl, h_1VecCmplx) = IFTTup(setzeroVecCmplx(H1VecCmplx, 1));
 tie(time_h_2VecDbl, h_2VecCmplx) = IFTTup(setzeroVecCmplx(H2VecCmplx, 0));
 // Q3 g - write h'1 and h'2 (h_1 and h_2) to file
 printfile("h'1sample.csv", time_h_1VecDbl, h_1VecCmplx);
 printfile("h'2sample.csv", time_h_2VecDbl, h_2VecCmplx);
 // Q3 i - get file in with h3 data where N = 200
 NInt = 200;
 vector<int> nh3VecInt; 
 vector<double> timeh3VecDbl; 
 vector<complex> h3VecCmplx;
 tie(nh3VecInt, timeh3VecDbl, h3VecCmplx) = getinfoTup("h3.txt");
 
 // Q3 j - inverse Fourier transform h3 to get H3
 vector<complex> H3VecCmplx;
 vector<double> omega_H3VecDbl;
 tie(omega_H3VecDbl, H3VecCmplx) = DFTTup(h3VecCmplx);
 // Q3 k - find 4 greatest values of mod H3 apply fourier transfor to these
 //find the index of the 4 greatest amplitude of mod H3
 int maxneeded = 4;
 vector<complex> maxH3VecCmplx = maxFVecCmplx(maxneeded, H3VecCmplx);
 // apply fourier transform to 4 greatest values
 vector<complex> h_3VecCmplx; 
 vector<double>time_h_3VecDbl;
 tie(time_h_3VecDbl, h_3VecCmplx) = IFTTup(maxH3VecCmplx);
 //Q3 l - write outcome for h'3 in file
 printfile("h'3sample.csv", time_h_3VecDbl, h_3VecCmplx);
 
 return 0;
}
