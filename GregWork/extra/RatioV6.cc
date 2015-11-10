/*The purpose of this code is to calculate a parameter which indicates the oscillations 
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

struct DataPoint
{
  double qval;
  double Inten;
  double Err;
};


/*START OF MAIN PROGRAM*/

int main(int argc, char *argv[])
{
  //This assumes the analysis is conducted on a particle size smaller than the ribosome
  double MaxParticleSize = 400;

  //The maximum the intensity can change for a particle  400 angstroms in size is 10^12. So I scale all data sets to 10^12 and remove all data points which are below 1.
  double IzeroSetPoint = 1*10e12;
  

// Input Checking
  if (argc !=3 )
    {
      cerr << "Usage: " << argv[0] << " ScatteringFile1 ScatteringFile2\n";
      exit(EXIT_FAILURE);
    }

  ifstream Scat1;
  ifstream Scat2;
  
  Scat1.open(argv[1]);
  Scat2.open(argv[2]);

  if(!Scat1.is_open())
    cerr << "Could not open your first scattering file\n";

  if(!Scat2.is_open())
    cerr << "Could not open your second scattering file\n";
//End of Input Checking

//Reading Files, Dividing them by one another and outputting to a temp file  
  
  DataPoint DataScat1;
  DataPoint DataScat2;
  double AveDat1 = 0.0;
  double AveDat2 = 0.0;
  int n = 0;
  
  double MaxI1 = 0.0;
  double MaxI2 = 0.0;

  while(Scat1 >> DataScat1.qval >> DataScat1.Inten >> DataScat1.Err)
    { 
      Scat2 >> DataScat2.qval >> DataScat2.Inten >> DataScat2.Err;
      if(DataScat2.qval != DataScat1.qval)
	cerr << "The Two data files are not properly interpolated\n";
      if(DataScat1.Inten > MaxI1) MaxI1 = DataScat1.Inten;
      if(DataScat2.Inten > MaxI2) MaxI2 = DataScat2.Inten;
      AveDat1 += DataScat1.Inten;
      AveDat2 += DataScat2.Inten;
      n++;
    }
  AveDat1 = AveDat1 / n;
  AveDat2 = AveDat2 / n;

  double ScaleFactor1 = IzeroSetPoint/MaxI1;
  double ScaleFactor2 = IzeroSetPoint/MaxI2;

  Scat1.close();
  Scat2.close();

  ifstream Scat1_1;
  ifstream Scat2_1;
  
  Scat1_1.open(argv[1]);
  Scat2_1.open(argv[2]);
  double Volit= 0;
  double ratio,ratErr;

  char *RatFile1 = "RatioFile1";
  ofstream RATFILE1;
  RATFILE1.open(RatFile1);
  double center1 = 0;

  while(Scat1_1 >> DataScat1.qval >> DataScat1.Inten >> DataScat1.Err)
    { 
     Scat2_1 >> DataScat2.qval >> DataScat2.Inten >> DataScat2.Err;
     DataScat1.Inten = DataScat1.Inten*ScaleFactor1;
     DataScat2.Inten = DataScat2.Inten*ScaleFactor2;
     DataScat1.Err = DataScat1.Err*ScaleFactor1;
     DataScat2.Err = DataScat2.Err*ScaleFactor2;
    
//Here I make sure that the minimum value is physically possible given the 10^12 decay discussed above.
     if(DataScat1.Inten < 1) DataScat1.Inten = 1;
     if(DataScat2.Inten < 1) DataScat2.Inten = 1;  

     ratio = (DataScat1.Inten)/(DataScat2.Inten);
     /*     if(ratio < 1)
       ratio = 1/ratio;
     */    
     //      cout << MaxI2 <<"\n";
     ratErr = ratio*sqrt(pow((DataScat1.Err/DataScat1.Inten),2)+(pow((DataScat2.Err/DataScat2.Inten),2)));
     center1 += ratio;
     RATFILE1 << DataScat1.qval << " " << ratio << " " << ratErr << "\n";
    } 

  center1 = center1/n;

  RATFILE1.close();
  ifstream RATFILE2;
  RATFILE2.open(RatFile1);
  double step1;

  char *RatFile3 = "RatioFile2";
  ofstream RATFILE3;
  RATFILE3.open(RatFile3);

  int i,j,k;
  DataPoint RatioPoint[n];
  
  for(i = 0; i < n; i++)
    {   
      RATFILE2 >> RatioPoint[i].qval >> RatioPoint[i].Inten >> RatioPoint[i].Err;
    }

  for(i = 0; i < n; i++)
    {
      RatioPoint[i].Inten = RatioPoint[i].Inten/center1;
      RatioPoint[i].Err = RatioPoint[i].Err/center1;
      RATFILE3 << RatioPoint[i].qval << " " << RatioPoint[i].Inten << " " <<  RatioPoint[i].Err << "\n";
    }

  double qmin = RatioPoint[0].qval;
  double qmax = RatioPoint[n-1].qval;

  

  double deltaq = 3.14159/MaxParticleSize;
  int NumBox =  int ((qmax/deltaq) + 0.5);
  
  double StandDevDiff[NumBox];
  double AveBox[NumBox];
  double BinnedQval[NumBox];
  double Counter[NumBox];
  double BoxErr[NumBox];
  int BoxNum;

  //  cout << qmin/deltaq << "\n";

  for(k =0; k< NumBox;k++)
    {
      AveBox[k] = 0;
      BinnedQval[k] = 0;
      Counter[k] = 0;
      BoxErr[k] = 0;
    }
  
  for(i = 0; i < n; i++)
    {
      BoxNum = int (RatioPoint[i].qval/deltaq);
      AveBox[BoxNum] += RatioPoint[i].Inten;
      Counter[BoxNum]++;
//      cout << BoxNum << " " << AveBox[BoxNum] << "\n";
  }
 
  for(i = 0; i < NumBox; i++)
    {
      if(Counter[i] > 0)
	{
	  BinnedQval[i] = (i*deltaq)+(0.5*deltaq); 
	  AveBox[i] = AveBox[i]/Counter[i];
	  //      cout << BinnedQval[i] << " " << i << "\n";
	}
    }
  
 for(i = 0; i < n; i++)
    {
      BoxNum = int ((RatioPoint[i].qval/deltaq) + 0.5);
      BoxErr[BoxNum] += (AveBox[BoxNum] - RatioPoint[i].Inten) * (AveBox[BoxNum] - RatioPoint[i].Inten); 
    }

 for(i = 0; i < NumBox; i++)
   {
     if(Counter[i] > 0)
       {
	 BoxErr[i] = sqrt(BoxErr[i]/Counter[i]);
       }
   }

  char *RatFile4 = "RatioFile3";
  ofstream RATFILE4;
  RATFILE4.open(RatFile4);

  for(k = 0; k < NumBox; k++)
    {
      if(Counter[k] > 0 && BinnedQval[k] > qmin && BinnedQval[k] < qmax)
	{
	  RATFILE4 <<  BinnedQval[k] << " " << AveBox[k] << " " << BoxErr[k] << "\n";
	}
    }

  double Hit = 0;
  double step2;

  for(k = 1; k < NumBox-1; k++)
    {
      if(Counter[k] > 0 && BinnedQval[k] > qmin && BinnedQval[k] < qmax)
	{
	  step1 = (abs(AveBox[k] - AveBox[k+1]));
	  step2 = (AveBox[k] + AveBox[k+1])/2;
	  /*      cout << Hit <<"\n"; */
	  Hit += step1/step2;
	}
    }

  cout << 100*(Hit/NumBox);
  return 0;
}

