/*
 * =====================================================================================
 *
 *       Filename:  colorCalculation.hpp
 *
 *    Description:  :w
 *
 *
 *        Version:  1.0
 *        Created:  14.09.2012 14:19:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *	
 *
 *         Author:  Michal Zychowski (), michal@fotonowy.pl
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <memory>
#include <float.h>
//#include <assert.h>
//#include "spline.hpp"

using namespace std;

#ifndef COLORCALCULATION_H
#define COLORCALCULATION_H

#define LERP(a,b,c)     (((b) - (a)) * (c) + (a))

struct XYZColor
{
	double X;
	double Y;
	double Z;
};

struct YuvColor
{
	double Y;
	double u;
	double v;
};


struct RGBColor
{
	unsigned char Red;
	unsigned char Green;
	unsigned char Blue;
};

struct LabColor
{
	double L;
	double a;
	double b;
};

struct Colors
{
	RGBColor RGB;
	XYZColor XYZ;
	LabColor Lab;
	YuvColor Yuv;
};
struct CD
{
	double c;
	double d;
};

typedef struct UVT {
        double  u;
        double  v;
        double  t;
} UVT;

struct TCSobj{
	double *Table;
	gsl_interp_accel *acc;
	gsl_spline *line;
};



class colorCalculation
{
	
	public:
	colorCalculation(vector <double> xScale)
	{
		initialObservers();
		initialTCS();
		
		for(unsigned int i=0;i<xScale.size();i++)
		{
			xscale.push_back(xScale[i]);
		}

		
		
	}
	colorCalculation()
	{
		initialObservers();
	}
	


	double dE1976(vector <double> oldspectrum,vector <double> newspectrum);
	double dE2000(vector <double> oldspectrum,vector <double> newspectrum);

	RGBColor colorFromWave(double wavelength);
	Colors calculateColors(vector <double>spectrum);
	double calculateColorTemperature(vector <double>spectrum);
	double calculateRobbertosnCCT(vector<double> spectrum);
	CD calculateCD(YuvColor Yuv);
	
	
	XYZColor calculateColorsXYZ(vector <double> spectrum);
	RGBColor calculateColorsRGB(XYZColor xyz);
	LabColor calculateColorsLab(XYZColor xyz);
	YuvColor calculateColorsYuv(XYZColor xyz);

	XYZColor getColorsXYZ();
	RGBColor getColorRGB();
	LabColor getColorLab();


	void initialTCS();

	void setXScale(vector <double> xScale);
	void getWhiteReference();
	void getWhiteReference(vector <double>spectrum);
	double blackBody(double temperature,double wavelenght);
	
	double calculateRa(vector <double>spectrum,double CTT);
	double calculateRa(vector <double>spectrum);

	~colorCalculation();  
	private:
	
	static const double h = 6.62606957e-34;	//planc const
	static const double c = 299792458;		//light speed
	static const double k_b= 1.3806488e-23;	//boltzman const
	
	vector <double> observerScale;
	vector <double> observerRed;
	vector <double> observerGreen;
	vector <double> observerBlue;
	vector <double> xscale;

	double redObserver(double wavelenght);
	double greenObserver(double wavelength);
	double blueObserver(double wavelength);

	void initialObservers();
	
	XYZColor XYZ;
	RGBColor RGB;
	LabColor Lab;
	
	double linearInterpolate(vector <double>xValues,vector <double>yValues ,double xValue); 

	XYZColor XYZWhiteReference;
	
	Colors colors;
	
	
	gsl_interp_accel *redacc;
	gsl_interp_accel *greenacc;
	gsl_interp_accel *blueacc;
        
	gsl_spline *redspline;
	gsl_spline *greenspline;
	gsl_spline *bluespline;

	
	double *x;
	double *yred;
	double *ygreen;
	double *yblue;

	double f(double t);
	double *TCSScale;
	TCSobj TCS[14];
	bool noFileFlag;

};

#endif
