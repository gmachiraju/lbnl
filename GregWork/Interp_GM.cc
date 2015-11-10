//Functions referenced from Numerical Recipes  -- Assembled and written by Gautam Machiraju

#include "nr3.h"
#include "interp_1d.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
using namespace std;

/*
**=========================================
**           Begin main method
**========================================= 
*/

// Adapted from the Numerical Recipes forum

int main() {

    int n2 = 0; 
	std::string line2;
    std::ifstream myfile2("RatioFile2"); //raw data points
    while (std::getline(myfile2, line2))
        ++n2; //getting number of lines in data file

	int n3 = 0;
	std::string line3;
    std::ifstream myfile3("RatioFile3"); //binned values
    while (std::getline(myfile3, line3))
        ++n3; //getting number of lines in data file

	// Spline_interp needs vecDoubs as input
    VecDoub xvd2(n2);
    VecDoub yvd2(n2);
   	VecDoub xvd3(n3);
    VecDoub yvd3(n3);

    std::ifstream x2("x_RatioFile2");
    std::istream_iterator<double> startx2(x2), endx2;
    std::vector<double> xv2(startx2, endx2);
    //std::cout << "Read " << xv2.size() << " numbers" << std::endl;

    std::ifstream y2("y_RatioFile2");
    std::istream_iterator<double> starty2(y2), endy2;
    std::vector<double> yv2(starty2, endy2);
    //std::cout << "Read " << yv2.size() << " numbers" << std::endl;

    std::ifstream x3("x_RatioFile3");
    std::istream_iterator<double> startx3(x3), endx3;
    std::vector<double> xv3(startx3, endx3);
    //std::cout << "Read " << xv3.size() << " numbers" << std::endl;

    std::ifstream y3("y_RatioFile3");
    std::istream_iterator<double> starty3(y3), endy3;
    std::vector<double> yv3(starty3, endy3);
    //std::cout << "Read " << yv3.size() << " numbers" << std::endl;

    int i = 0 ; 
    while ( i < n2 ) 
    {
        xvd2[i] = xv2[i] ;
        yvd2[i] = yv2[i] ;
        i++ ;
    }


    int j = 0 ;
    while ( j < n3 ) 
    {
        xvd3[j] = xv3[j] ;
        yvd3[j] = yv3[j] ;
        j++ ;   
    }

    // Function call
    Spline_interp si = Spline_interp(xvd2, yvd2);

    // Will print n values of y and interpolated value of y for xmin <= x <= xmax (RatioFile3 - binned data)
    char *Err_RatFile = "ErrorFile";          
  	ofstream ERRORFILE;                             
  	ERRORFILE.open(Err_RatFile);
 

    cout << "    x       y     interpolated     bin error    " << endl;
    cout << "------------------------------------------------" << endl;
    cout << fixed;
    for (int i = 0; i < n3; i++) {
        double xx3 = xv3[i];
        double yy3 = yv3[i];
        double binErr = abs(si.interp(xx3) - yy3);
        ERRORFILE << xx3 << " " << binErr  << "\n";
        cout << xx3 << " " << setw(10) << yy3 << " " << setw(10) <<  si.interp(xx3) << " " << setw(10) << binErr << endl;
    }

    cout << endl << "Your data set has been interpolated!" << endl << endl;

    //This automatic call to gnuplot doesn't work, so call manually.
	// if (!system("gnuplot gnuplot_GM.gpt")) {
	//     return 1;
 //    } else {
	//     return 0;
	// }

    return 0;
}

