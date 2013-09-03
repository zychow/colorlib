/*
 * =====================================================================================
 *
 *       Filename:  colorCalculation.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  13.09.2012 23:49:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michal Zychowski (), michal@fotonowy.pl
 *   Organization:  
 *
 * =====================================================================================
 */



#include "colorCalculation.hpp"

using namespace std;

colorCalculation::~colorCalculation()
{
         gsl_spline_free (redspline  );
         gsl_spline_free (greenspline);
         gsl_spline_free (bluespline );
         
	 gsl_interp_accel_free (redacc  );
         gsl_interp_accel_free (greenacc);
         gsl_interp_accel_free (blueacc );
	
	 for (int i=0;i<14;i++){
		gsl_spline_free(TCS[i].line);
		gsl_interp_accel_free(TCS[i].acc);
		free(TCS[i].Table);
	 }

	 free (TCSScale);
	 free (x);
	 free (yred);
	 free (ygreen);
	 free (yblue);

}

double colorCalculation::blackBody(double temperature,double wavelength)
{
	double wl = wavelength/1e9;
	
	double radiance = ((2*h*c*c)/(wl*wl*wl*wl*wl)) *1/(exp(h*c/(wl*k_b*temperature))-1);

	return radiance;
}

double colorCalculation::calculateColorTemperature(vector <double>spectrum)		// wikipedia cct
{
	XYZColor colors = calculateColorsXYZ(spectrum);

	//wikipedia color temperature approsximation 2856K to 6504K

	double xe = 0.3366;
	double ye = 0.1735;
	double A0 = -949.86315;
	double A1 = 6253.80338;
	double t1 = 0.92159;
	double A2 = 28.70599;
	double t2 = 0.20039;
	double A3 = 0.0004;
	double t3 = 0.07125;	
	
	
	double n = (colors.X - xe)/(colors.Y - ye);			

	double cct = A0+A1*exp(-n/t1)+A2*exp(-n/t2)+A3*exp(-n/t3);


	if (cct>1000&&cct<9000){
		return cct;
	}
	else{
		return -666.666;
	}
}

double colorCalculation::calculateRobbertosnCCT(vector <double>spectrum)	// roberson method from http://www.brucelindbloom.com/Eqn_XYZ_to_T.html
{
	XYZColor colors = calculateColorsXYZ(spectrum);


  double us, vs, p, di, dm;
        int i;

	
	double rt[31] = {       /* reciprocal temperature (K) */
	         DBL_MIN,  10.0e-6,  20.0e-6,  30.0e-6,  40.0e-6,  50.0e-6,
	         60.0e-6,  70.0e-6,  80.0e-6,  90.0e-6, 100.0e-6, 125.0e-6,
	        150.0e-6, 175.0e-6, 200.0e-6, 225.0e-6, 250.0e-6, 275.0e-6,
	        300.0e-6, 325.0e-6, 350.0e-6, 375.0e-6, 400.0e-6, 425.0e-6,
	        450.0e-6, 475.0e-6, 500.0e-6, 525.0e-6, 550.0e-6, 575.0e-6,
	        600.0e-6
	};

	UVT uvt[31] = {
	        {0.18006, 0.26352, -0.24341},
	        {0.18066, 0.26589, -0.25479},
	        {0.18133, 0.26846, -0.26876},
	        {0.18208, 0.27119, -0.28539},
	        {0.18293, 0.27407, -0.30470},
	        {0.18388, 0.27709, -0.32675},
	        {0.18494, 0.28021, -0.35156},
	        {0.18611, 0.28342, -0.37915},
	        {0.18740, 0.28668, -0.40955},
	        {0.18880, 0.28997, -0.44278},
	        {0.19032, 0.29326, -0.47888},
	        {0.19462, 0.30141, -0.58204},
	        {0.19962, 0.30921, -0.70471},
	        {0.20525, 0.31647, -0.84901},
	        {0.21142, 0.32312, -1.0182},
	        {0.21807, 0.32909, -1.2168},
	        {0.22511, 0.33439, -1.4512},
	        {0.23247, 0.33904, -1.7298},
	        {0.24010, 0.34308, -2.0637},
	        {0.24792, 0.34655, -2.4681},	/* Note: 0.24792 is a corrected value for the error found in W&S as 0.24702 */
	        {0.25591, 0.34951, -2.9641},
	        {0.26400, 0.35200, -3.5814},
	        {0.27218, 0.35407, -4.3633},
	        {0.28039, 0.35577, -5.3762},
	        {0.28863, 0.35714, -6.7262},
	        {0.29685, 0.35823, -8.5955},
	        {0.30505, 0.35907, -11.324},
	        {0.31320, 0.35968, -15.628},
	        {0.32129, 0.36011, -23.325},
	        {0.32931, 0.36038, -40.770},
	        {0.33724, 0.36051, -116.45}
	};


        if ((colors.X < 1.0e-20) && (colors.Y < 1.0e-20) && (colors.Z < 1.0e-20))
                return(-666.666);     /* protect against possible divide-by-zero failure */
        us = (4.0 * colors.X) / (colors.X + 15.0 * colors.Y + 3.0 * colors.Z);
        vs = (6.0 * colors.Y) / (colors.X + 15.0 * colors.Y + 3.0 * colors.Z);
        dm = 0.0;
        for (i = 0; i < 31; i++) {
                di = (vs - uvt[i].v) - uvt[i].t * (us - uvt[i].u);
                if ((i > 0) && (((di < 0.0) && (dm >= 0.0)) || ((di >= 0.0) && (dm < 0.0))))
                        break;  /* found lines bounding (us, vs) : i-1 and i */
                dm = di;
        }
        if (i == 31)
                return(-666.666);     /* bad XYZ input, color temp would be less than minimum of 1666.7 degrees, or too far towards blue */
        di = di / sqrt(1.0 + uvt[i    ].t * uvt[i    ].t);
        dm = dm / sqrt(1.0 + uvt[i - 1].t * uvt[i - 1].t);
        p = dm / (dm - di);     /* p = interpolation parameter, 0.0 : i-1, 1.0 : i */
        p = 1.0 / (LERP(rt[i - 1], rt[i], p));
//        *temp = p;
        return(p);      /* success */	

}

void colorCalculation::initialTCS()
{
	FILE *tcsFile = NULL;
	if(!(tcsFile = fopen("./TCS.dat","r"))){
		noFileFlag = true;
		return;
	}
//	printf("%d\n",tcsFile);
	int i=0;

	if(!(TCSScale = (double *)malloc(94*sizeof(double)))){
		noFileFlag= true;
		printf("no memory");
		return;
	}

	for (i=0;i<14;i++){
		if(!(TCS[i].Table = (double *)malloc(94*sizeof(double)))){
			noFileFlag = true;
			printf("no memory");
			return;
		}
	}

	i=0;
//	printf("reading:\n");
	double tmp;
	
	for (i=0;i<94;++i){
		fscanf(tcsFile,"%lf",&tmp);
		TCSScale[i]=tmp;
		for (int j=0;j<14;++j){
			fscanf(tcsFile,"%lf",&tmp);
			TCS[j].Table[i]=tmp;
		}
	}

	for (int i=0;i<14;i++){
		TCS[i].acc = gsl_interp_accel_alloc();
		TCS[i].line  = gsl_spline_alloc(gsl_interp_cspline, 94);
        	gsl_spline_init(TCS[i].line,TCSScale, TCS[i].Table, 94);
	}
	noFileFlag= false;

}


void colorCalculation::initialObservers()
{

      		ifstream myfile;
    		string line;
		stringstream tempstring;
		double temp;
		
		myfile.open ("observerScale.dat");
    		if (myfile.is_open())
		{
			while(myfile.good())
			{
				temp = 0;
				line.clear();
				tempstring.clear();
				getline(myfile,line);
				
				tempstring<<line;
				tempstring>>temp;
//				cout<<temp<<endl;
				if(temp>0){
					observerScale.push_back(temp);
				}
			}

		}else{
			cout<<"no file"<<endl;
		}
		myfile.close();




    		myfile.open ("observerRed.dat");
    		if (myfile.is_open())
		{
			while(myfile.good())
			{
				temp = 0;
				line.clear();
				tempstring.clear();
				getline(myfile,line);
				tempstring<<line;
				tempstring>>temp;
//				cout<<temp<<endl;
				if(observerRed.size()<observerScale.size())
					observerRed.push_back(temp);
			}

		}	
		myfile.close();
    		myfile.open ("observerGreen.dat");
    		if (myfile.is_open())
		{
			while(myfile.good())
			{
				temp = 0;
				line.clear();
				tempstring.clear();
				getline(myfile,line);
				tempstring<<line;
				tempstring>>temp;
//				cout<<temp<<endl;
				if(observerGreen.size()<observerScale.size())
					observerGreen.push_back(temp);
			}

		}	
		myfile.close();

		
    		myfile.open ("observerBlue.dat");
    		if (myfile.is_open())
		{
			while(myfile.good())
			{
				temp = 0;
				line.clear();
				tempstring.clear();
				getline(myfile,line);
				tempstring<<line;
				tempstring>>temp;
//				cout<<temp<<endl;
				if(observerBlue.size()<observerScale.size())
					observerBlue.push_back(temp);
			}

		}	
		myfile.close();



	x 	= (double*)malloc(sizeof(double)*observerScale.size());
//	yred 	= new double(observerScale.size());
//	ygreen 	= new double(observerScale.size());
//	yblue 	= new double(observerScale.size());
	yred 	= (double*)malloc(sizeof(double)*observerScale.size());
	ygreen 	= (double*)malloc(sizeof(double)*observerScale.size());
	yblue 	= (double*)malloc(sizeof(double)*observerScale.size());

  	for(unsigned int i=0;i<observerScale.size();i++)
	{
		x[i] = observerScale[i];
		yred[i] = observerRed[i];
		ygreen[i] = observerGreen[i];
		yblue[i] = observerBlue[i];
//		cout<<x[i]<<"\t"<<yred[i]<<"\t"<<ygreen[i]<<"\t"<<yblue[i]<<endl;

	}

	redacc = gsl_interp_accel_alloc();
        greenacc= gsl_interp_accel_alloc();
        blueacc= gsl_interp_accel_alloc();

	redspline  = gsl_spline_alloc(gsl_interp_cspline, observerScale.size());
        greenspline= gsl_spline_alloc(gsl_interp_cspline, observerScale.size());
        bluespline = gsl_spline_alloc(gsl_interp_cspline, observerScale.size());

        gsl_spline_init(redspline, x, yred, observerScale.size());
        gsl_spline_init(greenspline, x, ygreen, observerScale.size());
        gsl_spline_init(bluespline, x, yblue, observerScale.size());

//	getWhiteReference();

}

void colorCalculation::getWhiteReference()
{
	string substring;
	stringstream temp;
	double value;
	vector <double>whiteRef;
	vector <double>xscale;
	fstream myfile;
	myfile.open("whiteref.dat",fstream::in);
	
	while(!myfile.eof())
	{
		getline(myfile,substring,'\t');
		temp<<substring;
		temp>>value;
		xscale.push_back(value);
		substring.clear();
		temp.clear();
		getline(myfile,substring,'\r');
		temp<<substring;
		temp>>value;
		cout<<value<<endl;
		whiteRef.push_back(value);
		substring.clear();
		temp.clear();
		cout<<substring<<endl;
		substring.clear();
	}
	
	XYZWhiteReference = calculateColorsXYZ(whiteRef);

	myfile.close();
}


void colorCalculation::setXScale(vector <double>xScale)
{
	xscale.clear();

	
		for(unsigned int i=0;i<xScale.size();i++)
		{
			xscale.push_back(xScale[i]);
		}

}

Colors colorCalculation::calculateColors(vector <double>spectrum)
{

	colors.XYZ =calculateColorsXYZ(spectrum);
	colors.RGB=RGB;
	colors.Lab=Lab;
	colors.Yuv= calculateColorsYuv(colors.XYZ);
	return colors;
	
}
double colorCalculation::calculateRa(vector <double>spectrum)
{
	double cct = calculateRobbertosnCCT(spectrum);


	return calculateRa(spectrum,cct);
}
YuvColor colorCalculation::calculateColorsYuv(XYZColor xyz)
{
	YuvColor Yuv;
	Yuv.Y=xyz.Y;
	Yuv.u= 4*xyz.X/(xyz.X+15*xyz.Y+3*xyz.Z);
	Yuv.v=6*xyz.Y/(xyz.X+15*xyz.Y+3*xyz.Z);
	return Yuv;
	
}

CD colorCalculation::calculateCD(YuvColor Yuv)
{
	CD cd;
	cd.c = (4 - Yuv.u - 10*Yuv.v) / Yuv.v;
	cd.d = (1.708*Yuv.v + 0.404 - 1.481*Yuv.u) / Yuv.v;
	return cd;

}

double colorCalculation::calculateRa(vector <double> spectrum,double CCT)
{	

	int spectrumSize = spectrum.size();

	vector<vector <double> > refSpects;
	vector<vector <double> > spectSpects;
	vector<double> blackB;
	double bbvalue;
	CD CDr[14];
	CD CDk[14];
	YuvColor kip[14];
	double Wstarr[14];
	double Ustarr[14];
	double Vstarr[14];
	
	double Wstark[14];
	double Ustark[14];
	double Vstark[14];

	double Ckm,Dkm;		// C,D of pure spectrum
	double Crm,Drm;		// C,D of pure black body spectrum
	double colorDistance;	// distance from color and black body
	double Yknormal,Yrnormal; //y to nromalization reference and true spectrum
	double deltaE[14];
	double R[14];

	Colors refColors[14];		//colors of all reference samples
	Colors spectColors[14];		//colors of spects samples
	Colors mainSpectrumColor;	
	Colors referenceSpectrumColor;

	refSpects.resize(spectrumSize);
	spectSpects.resize(spectrumSize);
	refSpects.resize(14);
	spectSpects.resize(14);
	for (unsigned int i=0;i<refSpects.size();++i){
		refSpects[i].resize(spectrumSize);
		spectSpects[i].resize(spectrumSize);
	}

	if(noFileFlag){
		printf ("nofile\n");
		return -666.666;
	}

	for (int i=0;i<spectrumSize;++i){
		bbvalue=blackBody(CCT,xscale[i]);
		for (int j=0;j<14;++j){
			if((xscale[i]<825) &&(xscale[i]>360)){
				spectSpects[j][i] = spectrum[i]*gsl_spline_eval(TCS[j].line,xscale[i],TCS[j].acc);
				refSpects[j][i] = bbvalue*gsl_spline_eval(TCS[j].line,xscale[i],TCS[j].acc);
			}else{
				spectSpects[j][i] = 0.0;
				refSpects[j][i] = 0.0;
			}
	//		blackB.push_back(blackBody(xscale[i],CCT));
		}
		blackB.push_back(bbvalue);
	}

	printf("colors:\n");
	for (int j=0;j<14;++j){
		refColors[j] = calculateColors(refSpects[j]);
		spectColors[j]= calculateColors(spectSpects[j]);
//		printf("refcolors %f %f %f\n",refColors[j].XYZ.X,refColors[j].XYZ.Y,refColors[j].XYZ.Z);
//		printf("spectColors %f %f %f\n",spectColors[j].XYZ.X,spectColors[j].XYZ.Y,spectColors[j].XYZ.Z);
	}
	printf("\n");

	mainSpectrumColor = calculateColors(spectrum);
	referenceSpectrumColor = calculateColors(blackB);
	Yknormal = 100/mainSpectrumColor.XYZ.Y;
	Yrnormal = 100/referenceSpectrumColor.XYZ.Y;
	{
		CD tmp = calculateCD(mainSpectrumColor.Yuv);
		Ckm = tmp.c;
		Dkm = tmp.d;
		tmp = calculateCD(referenceSpectrumColor.Yuv);
		Crm = tmp.c;
		Drm = tmp.d;
	}
	double tmpKw;
	for (int i=0;i<14;i++){
		CDk[i].c = (4- spectColors[i].Yuv.u - 10*spectColors[i].Yuv.v)/spectColors[i].Yuv.v;
		CDk[i].d = (1.708*spectColors[i].Yuv.v +0.404 - 1.481*spectColors[i].Yuv.u) / spectColors[i].Yuv.v;
		kip[i].u = (10.872 + 0.404 * (Crm / (Ckm*CDk[i].c) - 4*(Drm/(Drm*CDk[i].d))))/(16.518 + 1.481 * (Crm/(Ckm*CDk[i].c)) - Drm/(Dkm*CDk[i].d));
		kip[i].v = 5.520/(16.518 + 1.481 * (Crm/(Ckm*CDk[i].c)) - Drm/(Dkm*CDk[i].d));
	
		Wstarr[i] = pow(25*(refColors[i].XYZ.Y*Yrnormal),0.3333) - 17;
		Ustarr[i] = 13*Wstarr[i]*(refColors[i].Yuv.u - referenceSpectrumColor.Yuv.u);
		Vstarr[i] = 13*Wstarr[i]*(refColors[i].Yuv.v - referenceSpectrumColor.Yuv.v);
	
		Wstark[i] = pow(25*(spectColors[i].XYZ.Y*Yknormal),0.3333) - 17;
		Ustark[i] = 13*Wstark[i]*(spectColors[i].Yuv.u - mainSpectrumColor.Yuv.u);
		Vstark[i] = 13*Wstark[i]*(spectColors[i].Yuv.v - mainSpectrumColor.Yuv.v);
		tmpKw = (Ustarr[i]-Ustark[i])*(Ustarr[i]-Ustark[i]) + (Vstarr[i]-Vstark[i])*(Vstarr[i]-Vstark[i]) +  (Wstarr[i]-Wstark[i])*(Wstarr[i]-Wstark[i]);
		
		deltaE[i] = sqrt(tmpKw>=0?tmpKw:0 );
	
		R[i] = 100-4.6 *deltaE[i];
	}
	double Ra;
	for (int i=0;i<8;++i){
		Ra +=R[i];
	}
	Ra/=8;
	
	return Ra;

//	DC = sqrt((uk-ur).^2 + (vk-vr).^2);
}


XYZColor colorCalculation::calculateColorsXYZ(vector <double> spectrum)
{
	double x=0;
	double y=0;
	double z=0;
	int observeriterator=0;
	
	while(xscale[0]>observerScale[observeriterator])
	{
		observeriterator++;
	}
//	cout<< observeriterator<<"\t"<<observerScale[observeriterator] <<endl;

	//////////////////calculate red observer responce	//working as long as you have xscale resolution beter than 5nm
	for(unsigned int i =0 ;i<spectrum.size()&&xscale[i]<750;i++)
	{
		if(xscale[i]>observerScale[observeriterator])
		{
			observeriterator++;
		}		
		
		x+=spectrum[i]*observerRed[observeriterator];
		y+=spectrum[i]*observerGreen[observeriterator];
		z+=spectrum[i]*observerBlue[observeriterator];
//		xsum+=observerRed[observeriterator];
//		ysum+=observerRed[observeriterator];
//		zsum+=observerRed[observeriterator];

	}	
	//////////////////////////////////////////////////////////////
	

	XYZColor xyz;
	if(x==0&&y==0&&z==0){
		xyz.X=0;
		xyz.Y=0;
		xyz.Z=0;
		return xyz;
	}


	xyz.X= x/(x+y+z);
	xyz.Y= y/(x+y+z);
	xyz.Z= z/(x+y+z);
//	cout<<xyz.X<<"\t"<<xyz.Y<<"\t"<<xyz.Z<<endl;

	return xyz;
}



double colorCalculation::redObserver(double wavelength)
{
	if (redspline&&redacc)
		return gsl_spline_eval(redspline,wavelength,redacc);
	else 
		return 0.0;
}


double colorCalculation::greenObserver(double wavelength)
{
	if (greenspline&&greenacc)
		return gsl_spline_eval(greenspline,wavelength,greenacc);
	else 
		return 0.0;
}

double colorCalculation::blueObserver(double wavelength)
{
	if (bluespline&&blueacc)
		return gsl_spline_eval(bluespline,wavelength,blueacc);
	else 
		return 0.0;

}

double colorCalculation::dE1976(vector <double>oldspectrum,vector<double>newspectrum)
{
	XYZColor xyz=calculateColorsXYZ(oldspectrum);
	LabColor LabOld= calculateColorsLab(xyz);
	LabColor LabNew=calculateColorsLab(calculateColorsXYZ(newspectrum));
	double temp = (LabNew.L-LabOld.L)*(LabNew.L-LabOld.L)+(LabNew.a-LabOld.a)*(LabNew.a-LabOld.a)+(LabNew.b-LabOld.b)*(LabNew.b-LabOld.b);
	return sqrt(temp);	

}

LabColor colorCalculation::calculateColorsLab(XYZColor xyz)
{
		LabColor color;
		color.L=116*(f(xyz.Y/XYZWhiteReference.Y)-16);
		color.a=500*(f(xyz.X/XYZWhiteReference.X)-f(xyz.Y/XYZWhiteReference.Y));
		color.b=200*(f(xyz.Y/XYZWhiteReference.Y)-f(xyz.Z/XYZWhiteReference.Z));
		return color;
}


double colorCalculation::dE2000(vector <double>oldspectrum,vector<double>newspectrum)
{
	XYZColor xyzold=calculateColorsXYZ(oldspectrum);
//	LabColor LabOld= calculateColorsLab(xyz);
	XYZColor xyznew=calculateColorsXYZ(newspectrum);

	double L1,L2,a1,a2,b1,b2,dL,medL,C1,C2,dC,medC,a1prim,a2prim,h1,h2,dh,dH,medH,T,Sl,Sc,Sh,Rt,dE; //Cab
	
	if (XYZWhiteReference.Y>0&&XYZWhiteReference.Z>0)
	{
		L1=116*f(xyzold.Y/XYZWhiteReference.Y)-16;
		L2=116*f(xyznew.Y/XYZWhiteReference.Y)-16;
		a1=500*(f(xyznew.X/XYZWhiteReference.X)-f(xyznew.Y/XYZWhiteReference.Y));
		a2=500*(f(xyzold.X/XYZWhiteReference.X)-f(xyzold.Y/XYZWhiteReference.Y));
		b1=200*(f(xyznew.Y/XYZWhiteReference.Y)-f(xyznew.Z/XYZWhiteReference.Z));
		b2=200*(f(xyzold.Y/XYZWhiteReference.Y)-f(xyzold.Z/XYZWhiteReference.Z));
		dL=L2-L1;

		medL =(L2+L1)/2;
		C1 	= sqrt (a1*a1+b1*b1);
		C2	= sqrt (a2*a2+b2*b2);
//		Cab	= C1-C2;
		medC	= (C1+C2)/2;
		dC	= C1-C2>0?C1-C2:C2-C1; 			//check it!!
		
		a1prim= a1+a1/2*(1 - sqrt(pow(medC,7)/(pow(medC,7)*pow(25,7))));
		a2prim= a2+a2/2*(1 - sqrt(pow(medC,7)/(pow(medC,7)*pow(25,7))));
		
	
		h1=atan(b1/a1prim);
		h2=atan(b2/a2prim);
		
		if (h1<0)
		{
			h1+=2*M_PI;
		}
		if (h2<0)
		{
			h2+=2*M_PI;
		}


		if (fabs(h1-h2)<=M_PI)
		{
			dh= h2-h1;
		}
		else if(fabs(h1-h2)>M_PI&&h2<=h1)
		{
			dh = h2-h1 +2*M_PI;
		}
		else if(fabs(h1-h2)>M_PI&&h2>h1)
		{
			dh = h2-h1 -2*M_PI;
		}

		dH= 2*sqrt(C1*C2)*sin(dh/2);

		if(fabs(h1-h2)>M_PI)
		{
			medH = (h2+h1+2*M_PI)/2;	
		}
		else
		{
			medH= (h1+h2)/2;
		}

		T= 1-0.17*cos(medH-M_PI/6)+0.24*cos(2*medH)+0.32*cos(3*medH+M_PI/30)-0.20*cos(4*medH-21*M_PI/60);

		Sl = 1+((0.015*(medL-50)*(medL-50))/sqrt(20+(medL-50)*(medL-50)));
		Sc=1+0.045*medC;
		Sh=1+0.015*medC*T;
	
		Rt = -2*sqrt(pow(medC,7)/(pow(medC,7)*pow(25,7)))*sin(M_PI/3*exp(-1*(((dH*(180/M_PI)-275)/25)*((dH*(180/M_PI)-275)/25))));
		dE = sqrt( (dL/Sl)*(dL/Sl) + (dC/Sc)*(dC/Sc)+(dH/Sh)*(dH/Sh) +Rt*(dC/Sc)*dH/Sh);

		return dE;
	}

	return -666.666;

}

double colorCalculation::f(double t)
{
	double temp=0;

	if (t>(6/29)*(6/29)*(6/29))
	{
		temp = pow(t,(1/3));
	}
	else
	{
		temp = 1/3*(29/6)*(29/6)*t+(4/29);
	}
	return temp;
}

double linearInterpolate(vector <double>xValues,vector <double>yValues ,double xValue)
{
	double solve;

	if(xValues.size()!=yValues.size())
	{
		cout<<"bad size x y"<<endl;
		return 0;
	}



	if(xValues[0]<xValues[1])			//if grownig
	{
		unsigned int i=0;
		while(xValues[i]<xValue)
		{
			i++;
			if(i>=xValues.size()-2)
			{
				cout<<"out of range";

			}
		}

		double xproporc = xValues[i]-xValues[i-1];
		double yproporc = yValues[i]-yValues[i-1]; 	
		solve = yValues[i-1]+xproporc*yproporc;



	}
	else if(xValues[0]>xValues[1])			//oposite
	{
		unsigned int i=0;
		while(xValues[i]>xValue)
		{
			i++;
			if(i>=xValues.size()-2)
			{
				cout<<"out of range";

			}
		}





	}
	else
	{
		cout<<"bad xsize"<<endl;		//bad xValues vector
		return 0;
	}

	return solve;

}


/*  
Public Function lstar(Y, Yn)
    ratio = Y / Yn
    If ratio > 0.008856 Then
        lstar = 116 * (ratio ^ (1 / 3)) - 16
    Else
        lstar = 903.3 * ratio
    End If
End Function

Public Function astar(X, Y, Xn, Yn)
    xratio = X / Xn
    yratio = Y / Yn
    third = 1 / 3
    
    If xratio > 0.008856 Then
        fx = xratio ^ third
    Else
        fx = 7.787 * xratio + (16 / 116)
    End If

    If yratio > 0.008856 Then
        fy = yratio ^ third
    Else
        fy = 7.787 * yratio + (16 / 116)
    End If
    
    astar = 500 * (fx - fy)
    
End Function

Public Function bstar(Y, Z, Yn, Zn)
    
    yratio = Y / Yn
    zratio = Z / Zn
    third = 1 / 3
    
    If zratio > 0.008856 Then
        fz = zratio ^ third
    Else
        fz = 7.787 * zratio + (16 / 116)
    End If

    If yratio > 0.008856 Then
        fy = yratio ^ third
    Else
        fy = 7.787 * yratio + (16 / 116)
    End If
    
    bstar = 200 * (fy - fz)
    
End Function

Public Function deltaE2000(l1, a1, b1, l2, a2, b2)

    ' the usual chroma calculation
    cab1 = (a1 ^ 2 + b1 ^ 2) ^ 0.5
    cab2 = (a2 ^ 2 + b2 ^ 2) ^ 0.5

    meanc = (cab1 + cab2) / 2#
    
    ' this bit of nonsense is required because VB seems to choke on the big numbers
    logmeanc = Application.Log10(meanc)
    log25 = Application.Log10(25)
    cratio = (10 ^ (logmeanc * 7 - Application.Log10((10 ^ (logmeanc * 7) + 10 ^ (log25 * 7))))) ^ 0.5
    g = 0.5 * (1 - cratio)

    ' adjust a* values
    a1 = (1 + g) * a1
    a2 = (1 + g) * a2

    cab1 = (a1 ^ 2 + b1 ^ 2#) ^ 0.5
    cab2 = (a2 ^ 2 + b2 ^ 2#) ^ 0.5

    deltac = cab1 - cab2
    
    ' calculate delta L*
    deltal = l1 - l2

    ';calculate hue angles
    h1 = Application.Degrees(Application.Atan2(a1, b1))
    h2 = Application.Degrees(Application.Atan2(a2, b2))
    If h1 < 0 Then
        h1 = h1 + 360
        End If
    If h2 < 0 Then
        h2 = h2 + 360
    End If
    
    ' calculate delta h and mean h
    deltah = h1 - h2
    meanh = h1 + h2

    If deltah > 180 Then
        deltah = deltah - 360
        meanh = meanh - 360
    End If
    
    If deltah < -180 Then
        deltah = deltah + 360
        meanh = meanh - 360
    End If

    delth = 2 * ((cab1 * cab2) ^ 0.5) * Sin(Application.Radians(deltah / 2#))

    deltal = l1 - l2
    deltaE = (deltal ^ 2 + deltac ^ 2 + delth ^ 2) ^ 0.5

    deltaha = (deltaE ^ 2 - deltal ^ 2 - deltac ^ 2)
    If deltaha < 0 Then
        deltaha = 0
    End If
    
    delth = deltaha ^ 0.5
    delth = 2 * ((cab1 * cab2) ^ 0.5) * Sin(Application.Radians(deltah / 2))
    
    meanl = (l1 + l2) / 2#
    meanc = (cab1 + cab2) / 2#
    meanh = meanh / 2#


    ' Computes the weighting functions.

    sl = 1 + ((0.015 * (meanl - 50) ^ 2) / (20 + (meanl - 50) ^ 2) ^ 0.5)
    sc = 1 + 0.045 * meanc

    ' find the hue angle goodness
    ' split these up because I don't know how to split lines in vb
    t1 = 1 - 0.17 * Cos(Application.Radians(meanh - 30))
    t2 = 0.24 * Cos(Application.Radians(2 * meanh))
    t3 = 0.32 * Cos(Application.Radians(3 * meanh + 6))
    t4 = -0.2 * Cos(Application.Radians(4 * meanh - 63))
    t = t1 + t2 + t3 + t4
    
    sh = 1 + 0.015 * meanc * t
    sr = sc * sh
    
    ' same problem as above with the big numbers
    logmeanc = Application.Log10(meanc)
    log25 = Application.Log10(25)
    rc_numerator = 10 ^ (logmeanc * 7)
    rc_denominator = rc_numerator + 10 ^ (log25 * 7)
    rc = 2 * (rc_numerator / rc_denominator) ^ 0.5

    dthet = 30 * Exp(-1 * ((meanh - 275) / 25) ^ 2)
    rt = -Sin(2 * Application.Radians(dthet)) * rc

    ' Now computes deltaV
    vl = (deltal / sl) ^ 2
    vc = (deltac / sc) ^ 2
    vh = (delth / sh) ^ 2

    deltaE2000 = (vl + vc + vh + (rt * deltac * delth / sr)) ^ 0.5

;*/


/*  function  deltae2k, lab2, lab1, k

; delta E 2000 code created to match the spreadsheet from Klaus Witt
; DWyble 06Aug02

; If you have some time, make this work for 3xN lab vectors
; As it stands, you have to pass in two 3x1 lab triplets and an optional triplet of k's

    if (n_params() lt 3) then begin
		k = [1, 1, 1]
    endif

    meanL = (lab1(0)+lab2(0))/2.
    
    Cab1 = (lab1(1)^2.+lab1(2)^2.)^.5
    Cab2 = (lab2(1)^2.+lab2(2)^2.)^.5
    meanC = (Cab1 + Cab2)/2.
    
    G = 0.5*(1-((meanC^7.)/(meanC^7. + 25^7.))^0.5)
        
    aprime1=lab1(1)*(1+G)
    aprime2=lab2(1)*(1+G)
    
    CabPrime1 = (aprime1^2.+lab1(2)^2.)^.5
    CabPrime2 = (aprime2^2.+lab2(2)^2.)^.5
    meanCprime = (CabPrime1+CabPrime2)/2
    
    hprime1 = atan(lab1(2), aprime1)*180/!pi
    hprime2 = atan(lab2(2), aprime2)*180/!pi
    if hprime1 LT 0 then hprime1 = hprime1 + 360
    if hprime2 LT 0 then hprime2 = hprime2 + 360

    if ABS(hprime1-hprime2) GT 180 then $
	meanH = (hprime1+hprime2+360)/2. $
	else meanH = (hprime1+hprime2)/2.

    T = 1.-0.17*cos((meanH-30)*!pi/180) + $
		0.24*COS(2*meanH*!pi/180) + $
		0.32*COS((3.*meanH+6)*!pi/180) - $
		0.2*COS((4.*meanH-63)*!pi/180)

    if abs(hprime2 - hprime1) LE 180 then delhp = hprime2 - hprime1 $
    else begin
	if hprime2 EQ min([hprime1,hprime2]) then $
		delhp = hprime2 - hprime1+ 360 $
	else delhp = hprime2 - hprime1 - 360
    endelse
    dellp = lab2(0) - lab1(0)
    delcp = CabPrime2 - CabPrime1
    delCapHp=2*(CabPrime2*CabPrime1)^0.5*sin((!pi/180.)*delhp/2)

    Sl = 1.+(0.015*(meanL-50.)^2.)/((20+(meanL-50.)^2.)^0.5)
    Sc = 1.+0.045*meanCprime
    Sh = 1.+0.015*meanCprime*T

    delTheta = 30.*exp(-1.*((meanH-275.)/25.)^2)
    Rc = ((meanCprime^7.)/(meanCprime^7.+25.^7.))^0.5
    Rt = -Rc*2.*SIN(2.*delTheta*!pi/180.)

    lterm = (dellp/(k[0]*Sl))^2.
	cterm = (delcp/(k[1]*Sc))^2.
	hterm = (delCapHp/(k[2]*Sh))^2.
	rotate = Rt*(delcp/(k[1]*Sc))*(delCapHp/(k[2]*Sh))
	
    de2000 = ( lterm + cterm + hterm + rotate) ^ .5

    return, de2000
end
 *
*/


