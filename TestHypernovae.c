/* ***************************************************************
   * NEUCOSMA -- NEutrinos from Cosmic Accelerators              *
   * (c) 2010 NEUCOSMA team                          		 *
   *************************************************************** */

// Based on the GRB model from TestGRB.c, we here try to create a model for the collisions of Hypernovae. Even though there are photohadronic interactions inside the shock of a Hypernova during the Sedov-phase, pp-collisions should actually be the more relevant source of neutrinos. Hence, the standard photohadronic interactions used for our GRB calculations are here actually complemented by code for pp-collisions.

// Last modified: Mar. 12, 2014 by P. Baerwald

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "nco_boost.h"
#include "myio.h"
#include "nco_utils.h"
#include "nco_photo.h"
#include "nco_decays.h"
#include "nco_synchr.h"
#include "nco_steady.h"

#include <gsl/gsl_integration.h>


// THE PARAMETERS BELOW STILL NEED TO BE UPDATED FOR HYPERNOVAE!

//***************************************************************************************************************************
//***************************************************************************************************************************

// Calculation ranges protons/electrons and photons [GeV] - in the SRF!
static double EPMIN = 1;
static double EPMAX = 1e15;
static double EPHMIN; // The minimal and maximal photon energy are calculated in the function main() as the energy needs to be transformed from observer to SRF.
static double EPHMAX;

static double particlemass;
static double particlecharge;
static int particle;


// The parameters of the hypernova shell. Still need to figure out which the best ones are.

static double B; // Magnetic field in the SRF - in [G]; is calculated automatically in the program!
static double lorentz = 1.25; // Lorentz or Gamma factor - in [1]
static double viewingangle = 0.52;  // Viewing angle assumed to be 1/Gamma
static double redshift = 2.0;
static double luminosity = 1e38; // Luminosity at the source (but NOT in SRF) - in [erg/s]
static double duration = 10.0; // Duration (T90) at the observer - in [s]
static double tvar = 0.1; // Observed variability time scale (or dynamical time scale) - in [s]
static double grbradius;
static double shellthickness;
static double accelerationefficiency = 0.1; // Acceleration efficiency of the source. Can take values between 0 and 1. (We use 0.1 as standard value.)
static double fe = 0.01; // Energy in electrons as fraction of energy in protons.

static double photonbreak = 14.76e-6; // Parameter for break in photon spectrum - in SRF and in [GeV].
static double firstindex = 1.0; // Lower index of photon spectrum.
static double secondindex = 2.0; // Upper index of the photon spectrum.



//***************************************************************************************************************************
//***************************************************************************************************************************

static double maxprotonenergy = 1e9; // Arbitrary value for testing purposes!!!! NEEDS TO BE CHANGED!!!

static double Cgamma; // Normalization factor for the photon spectrum.
static double Cp; // Normalization factor for the proton spectrum.

static double intenergygamma; // Integrated energy in the photon spectrum.
static double intenergyp; // Integrated energy in the proton spectrum.
static double intlgammagamma; // Integrated (inverse) interaction length for pair creation from photons.
static double numofevents; // Number of events obtained by convolution of effective area, muon neutrino spectrum, and exposure time.

static double oscneubare[301];
static double oscneue[301];
static double oscneubarmu[301];
static double oscneumu[301];
static double oscneubartau[301];
static double oscneutau[301];


// Interpolating functions for later use.
nco_ip s_proton_steady;  // steady state
nco_ip s_neutron_steady;
nco_ip s_piplus_steady;
nco_ip s_piminus_steady;
nco_ip s_kplus_steady;
nco_ip s_muplusr_steady;
nco_ip s_muplusl_steady;
nco_ip s_muminusr_steady;
nco_ip s_muminusl_steady;
nco_ip s_piplus_in;
nco_ip s_piminus_in;
nco_ip s_kplus_in;
nco_ip s_neutron_in;
nco_ip s_proton_in;
nco_ip s_muplusr_in;
nco_ip s_muplusl_in;
nco_ip s_muminusr_in;
nco_ip s_muminusl_in;

nco_ip s_pcool;
nco_ip s_pesc;
nco_ip s_ncool;
nco_ip s_nesc;


// returns proton spectrum in 1/GeV 1/cm^3; energy in SRF in GeV; normalization 1
double grb_protons(double energy)
{
	if(energy>EPMIN) return 1.0/(energy*energy)*exp(-pow(energy/maxprotonenergy,2.0)); // Exponential cut-off at maxprotonenergy. The value of maxprotonenergy will be calculated in the programm.
		else return 0.0;
}

// Returns photon spectrum in 1/GeV 1/cm^3; energy in SRF in GeV; normalization 1. The parameters will be set during the run of the programm.
double grb_photons(double energy)
{
	if(energy<EPHMIN) return 0.0;
		else if(energy<=photonbreak) return pow(photonbreak/energy,firstindex); // Break at around 1keV in source frame (~250keV in observer's frame). Compare to Band paper - ApJ 418 281 (1993).
			else if(energy<=EPHMAX) return pow(photonbreak/energy,secondindex);
				else return 0.0;
}


// Wrapper for cooling and escape.
double w_pcool(double x)
{
	if(x>EPMAX || x<EPMIN) return 0.0;
		else return pow(10,ncoIP(&s_pcool,log10(x)));
}

double w_pesc(double x)
{
	if(x>EPMAX || x<EPMIN) return 0.0;
		else return pow(10,ncoIP(&s_pesc,log10(x)));
}

double w_ncool(double x)
{
	if(x>EPMAX || x<EPMIN) return 0.0;
		else return pow(10,ncoIP(&s_ncool,log10(x)));
}

double w_nesc(double x)
{
	if(x>EPMAX || x<EPMIN) return 0.0;
		else return pow(10,ncoIP(&s_nesc,log10(x)));
}


// Wrapper for energy loss.
double Eloss(double x) // x is related to energy. energy = pow(10,x).
{
	double phcooling=0.0;
	if(particle==NCO_NEUTRON) { return w_ncool(x);}
	else if(particle==NCO_PROTON) phcooling=w_pcool(x);  
	return phcooling+ncoComputeSynchrLossRate(particlemass, particlecharge, B, x)+(1.0+redshift)/(lorentz*tvar); // Version WITH adiabatic losses.
//	return phcooling+ncoComputeSynchrLossRate(particlemass, particlecharge, B, x); // Version without adiabatic losses.
}


// Wrapper for particle escape.
double Eescape(double x)
{
	double phescape=0.0;
	if(particle==NCO_PROTON) phescape=w_pesc(x);  
	if(particle==NCO_NEUTRON) phescape=w_nesc(x);
	return phescape+ncoComputeDecayEscapeRate(particle,x);
}


// Destroy all interpolating functions.
void ncoDestroySpectra()
{  
	ncoIPFree(&s_neutron_steady);
	ncoIPFree(&s_piplus_steady);
	ncoIPFree(&s_piminus_steady);
	ncoIPFree(&s_kplus_steady);
	ncoIPFree(&s_muplusr_steady);
	ncoIPFree(&s_muplusl_steady);
	ncoIPFree(&s_muminusr_steady);
	ncoIPFree(&s_muminusl_steady);
	ncoIPFree(&s_pcool);
	ncoIPFree(&s_pesc);
	ncoIPFree(&s_ncool);
	ncoIPFree(&s_nesc);
	ncoIPFree(&s_neutron_in);
	ncoIPFree(&s_proton_in);
	ncoIPFree(&s_piplus_in);
	ncoIPFree(&s_piminus_in);
	ncoIPFree(&s_kplus_in);
	ncoIPFree(&s_muplusr_in);
	ncoIPFree(&s_muplusl_in);
	ncoIPFree(&s_muminusr_in);
	ncoIPFree(&s_muminusl_in);
}


// The following functions are for the proton-proton interactions, analogous to the photohadronic interactions.

double inelcrosssection(double protonenergy) 	// The inelastic interaction cross section in [mb]=[10^-27 cm^2] as function of the high-energy proton's energy (in [GeV]).
{						// Based on eq. (73) of Kelner, Aharonian & Bugayov, astro-ph/0606058.
	double thresholdenergy = 0.938+2.0*0.135+0.135*0.135/(2.0*0.938); // Estimated threshold energy in pp-collisions for pi^0 production.
	
	if(protonenergy > thresholdenergy) {return (34.3+1.88*log(protonenergy*1e-3)+0.25*log(protonenergy*1e-3)*log(protonenergy*1e-3))*(1.0-pow((thresholdenergy/protonenergy),4.0))*(1.0-pow((thresholdenergy/protonenergy),4.0));} // Factor 1e-3 needed for conversion from [GeV] to [TeV]
	else {return 0.0;}
}

double Fppforpions(double x, double protonenergy) 	// Scaling function for the pion production with the fraction of pion energy compared to proton energy x=E/E_p
{							// and the proton energy E_p (in [GeV]).
	// Different intermediate variables also defined by Kelner, Aharonian & Bugayov, astro-ph/0606058.
	double a = 3.67+0.83*log(protonenergy*1e-3)+0.075*log(protonenergy*1e-3)*log(protonenergy*1e-3);
	double Bpi = a+0.25;
	double alpha = 0.98/sqrt(a);
	double r = 2.6/sqrt(a);
	
	if (x*protonenergy > 0.135) {return 4.0*alpha*Bpi*pow(x,alpha-1.0)*pow((1.0-pow(x,alpha))/(1.0+r*pow(x,alpha)*(1.0-pow(x,alpha))),4.0)*(1.0/(1.0-pow(x,alpha))+(r*(1.0-2.0*pow(x,alpha))/(1.0+r*pow(x,alpha)*(1.0-pow(x,alpha)))))*(1.0-0.135/(x*protonenergy));}
	else return 0.0;
}


// The following function is the "main" part of the calculation, namely the calculation of the spectra. The function calculates the particle interactions step by step starting at the protons and photons at the source. Currently only the oscillated neutrinos are printed to a file. Other spectra can also be obtained, however these functions are currently commented out. Moreover, the program also has the functionality to automatically generate (sorted) filenames, but then the functions needs the parameters which are currently commented out.

int ComputeNeutrinoSpectraFromScratch(double B, double lorentz, double viewingangle, double redshift, double Cp, double Cgamma, double grbradius, double shellthickness)
{
	double jetradius = 1e8*0.09; // Second factor is tangent of jet's half-opening angle. 0.03 for an assumed value of 1.5° or 0.09 for 5°. Jet's radius in [cm].
	double bulkvelocity = 3e10*sqrt(1.0-1.0/(lorentz*lorentz));
	double kineticenergydensity = luminosity/(M_PI*jetradius*jetradius*bulkvelocity);
	
	B = sqrt(8.0*M_PI*kineticenergydensity);
  
	double ncold = 0.9*624.15*luminosity/(lorentz*0.938*M_PI*jetradius*jetradius*bulkvelocity);
	
	double synchlimit = 2.00982e11*sqrt(accelerationefficiency/B);
	double esclimit = 1e-9*(lorentz*jetradius)*accelerationefficiency*3e8*(B*1e-4); // (eV->GeV)*(size of region in [m])*(acceleration efficiency)*(magnetic field in [T])

	if(synchlimit<esclimit) {maxprotonenergy = synchlimit;}
		else {maxprotonenergy = esclimit;}
	printf("E_p,max = %g\n",maxprotonenergy);
	
	
	double FGAMMA(double y, void* params)
	{
		return pow(10,y)*grb_photons(pow(10,y))*pow(10,y)*log(10.0); // New version for logarithmic integration
	}
	
	double FP(double y, void* params) // New version for logarithmic integration
	{
		return pow(10,y)*grb_protons(pow(10,y))*pow(10,y)*log(10.0); // First factor is energy, second is the proton spectrum, the third and fourth factors are needed for transition from linear to logarithmic scale.
	}
	
	const double RELERROR=1e-3;
	double spectraerrorest;
	int res, resgammagamma, res2;
	gsl_function fun,fun2;
	fun.function=&FGAMMA;
	fun.params=NULL;
	fun2.function=&FP;
	fun2.params=NULL;
	
	gsl_integration_workspace* WSspectra=gsl_integration_workspace_alloc(5000);
	
//************	Calculation of the normalization factor Cgamma - including the recursive adaption of the maximal photon energy that can escape! ******************************
	
	// Integration normally from 0.2 keV up to 30 MeV (upper limit of GBM on Fermi). The range 1 keV to 10 MeV is as in our IC40 reanalysis.
	double lowerintlimit = log10(0.2e-6*(1.0+redshift)/lorentz);
	double upperintlimit = log10(3e-2*(1.0+redshift)/lorentz); 
	double oldupperintlimit, taugammagamma;
	//Ephcutoff = 30e-3; // First test energy for cut-off: 30 MeV in SRF (before boosting).
	double Ephcutoff = 511e-6; // Second test energy for cut-off: 511 keV in SRF (before boosting).
	double step;
	int i,j;
	
	double FGAMMAGAMMA(double y) // Function needed for calculation of pair creation from photons. First part is photon spectrum, second part is "responce function" for the pair creation, and third part (last two factors) is the correction for logarithmic integration.
	{
		return grb_photons(pow(10,y))*(1.0-(511e-6*511e-6*511e-6*511e-6)/(pow(10,y)*pow(10,y)*Ephcutoff*Ephcutoff))*pow(10,y)*log(10.0); // New version for logarithmic integration
	}
	
	
	res= gsl_integration_qag(&fun,lowerintlimit,upperintlimit,0.0,RELERROR,5000,GSL_INTEG_GAUSS41,WSspectra,&intenergygamma,&spectraerrorest);
	if(res!=0)
	{
		printf("Int. error %i occured in normalization of photon spectrum\n",res);
		intenergygamma = 0.0;
	}

	Cgamma = fe*624.15*luminosity/(lorentz*intenergygamma*M_PI*jetradius*jetradius*bulkvelocity);
	
	oldupperintlimit = upperintlimit; // Define an "old" value for the upper limit for later use.
	
//		printf("%e %e %e %e %e\n",fractioninphotons,Eiso,Viso,lorentz,intenergygamma);
	
	step = -1.0;
	NCO_INTEGR_METHOD=NCO_GSL;
	
	for(j=0;j<2000;j++)
	{
	intlgammagamma = ncoIntegrate2(FGAMMAGAMMA,log10((511e-6*511e-6)/Ephcutoff),6.0,511); // Lower bound from pair creation threshold. Upper bound should be infinity, but set to 1 PeV in the SRF (!). Remember that the photon spectrum is only defined up to EPHMAX.
	
	taugammagamma = (lorentz*shellthickness*1e5)*(3.0/16.0*0.665e-24*Cgamma*intlgammagamma);
	
//		printf("tau_gammagamma (Ephcutoff = %g) = %g   in try %i\n",Ephcutoff,taugammagamma,j);
	
	if ( taugammagamma-1.0 < 1e-3 && taugammagamma-1.0 > -1e-3 ) {break;}
	else if ( ( taugammagamma-1.0 > 0.0 && step>0.0 ) || ( taugammagamma-1.0 < 0.0 && step<0.0 ) ) {step = (-step/2.0);}
	else if ( Ephcutoff>1e6 ) {Ephcutoff=1e6; break;} // Artificial cutoff at 1 PeV in the SRF to avoid infinities.
	else if ( Ephcutoff<(511e-6*511e-6)/1e6 ) {Ephcutoff=(511e-6*511e-6)/1e6; break;} // Artificial cut-off to ensure that the tau_gammagamma-integration is always larger than one.
	  
	Ephcutoff = pow(10.0,log10(Ephcutoff)+step);
	}
	
//		printf("tau_gammagamma (Ephcutoff = %g) = %g   in run %i\n",Ephcutoff,taugammagamma,i);
	printf("tau_gammagamma (Ephcutoff = %g) = %g   \n",Ephcutoff,taugammagamma);
	
//		EPHMAX=Ephcutoff; // Only photons below cutoff can escape - and are relevant targets for the escape. The rest cannot be correctly described with our approach.
	
	if( upperintlimit>log10(Ephcutoff) ){upperintlimit=log10(Ephcutoff);}
	else if( upperintlimit<log10(Ephcutoff) && upperintlimit<log10(3e-2*(1.0+redshift)/lorentz) )
	{
		if( Ephcutoff<0.2e-6*(1.0+redshift)/lorentz ){printf("Burst not visible in GBM range!\n"); upperintlimit=log10(3e-2*(1.0+redshift)/lorentz);}
		else if( Ephcutoff<3e-2*(1.0+redshift)/lorentz ){upperintlimit=log10(Ephcutoff);}
		else {upperintlimit=log10(3e-2*(1.0+redshift)/lorentz);}
	}

	
//************	Calculation of the normalization factor Cp - WITHOUT adapting EPMAX. *****************************************************************************************
	
// New approach with integration over logarithmic scale.

	res2= gsl_integration_qag(&fun2,0.0,11.0,0.0,RELERROR,5000,GSL_INTEG_GAUSS41,WSspectra,&intenergyp,&spectraerrorest);
	if(res2!=0)
	{
		printf("Int. error %i occured in normalization of proton spectrum\n",res2);
		intenergyp = 0.0;
	}
		
	gsl_integration_workspace_free(WSspectra);
		
	Cp = 0.1*624.15*luminosity/(lorentz*intenergyp*M_PI*jetradius*jetradius*bulkvelocity);
	printf("Cp = %g\n",Cp);
	
  
// Calculation of all decay and escape rates. Interpolating functions generated at the end.
	i=0;
	double xace[201];
	double yace[201];
	double ybce[201];
	double ycce[201];
	double ydce[201];
	double x,ex,pcool,pesc,ncool,nesc;
	for(x=log10(EPMIN);x<=10.01;x+=0.05)
	{
		xace[i]=x;
		ex=pow(10,x);
		ncoComputeCoolEscRate(NCO_ALL_PION,NCO_PROTON,grb_photons,ex,&pcool,&pesc);
		ncoComputeCoolEscRate(NCO_ALL_PION,NCO_NEUTRON,grb_photons,ex,&ncool,&nesc);
		if(pcool>1e-200) yace[i]=log10(pcool); else yace[i]=-200.0;
		if(pesc>1e-200) ybce[i]=log10(pesc); else ybce[i]=-200.0;
		if(ncool>1e-200) ycce[i]=log10(ncool); else ycce[i]=-200.0;
		if(nesc>1e-200) ydce[i]=log10(nesc); else ydce[i]=-200.0;
		i++;  
	}
	ncoIPAlloc(&s_pcool,xace,yace,i);
	ncoIPAlloc(&s_pesc,xace,ybce,i);
	ncoIPAlloc(&s_ncool,xace,ycce,i);
	ncoIPAlloc(&s_nesc,xace,ydce,i);


// Computation of the secondary particles from photohadronics and direct escape. Additionally with a basic calculation of the pp-interactions.
	
	double FPP(double x, void *params) // Variable x is secondary energy, e.g. ex, divided by the proton energy; x = ex/Ep.
	{
		double inputenergy = *(double *)(params);
		printf("x^-1 = %g, protons = %g, Fpp = %g, inel. x-section = %g\n",1.0/x,grb_protons(inputenergy/x),Fppforpions(x,inputenergy/x),inelcrosssection(inputenergy/x));
		return 1.0/x*Cp*grb_protons(inputenergy/x)*Fppforpions(x,inputenergy/x)*inelcrosssection(inputenergy/x);
	}
	
	double pperrorest;
	int ppres;
	gsl_function funpp;
	funpp.function=&FPP;
	
	gsl_integration_workspace* WSppint=gsl_integration_workspace_alloc(5000);
	
	
	const int PIONIT = NCO_ALL_PION; // Type of interaction - here: all pion interactions.
	double xa[301];
	double ya[301];
	double yb[301];
	double yc[301];
	double yd[301];
	double ypp[301];
	double ypfrompp[301];
	// mioInitOutput("piondata.dat"); // File name has to be set!
	double respp,respm,resk,resn,resp,respfrompp;
	double ppintegration;
	int n=0;
	for(x=0.0;x<=15.01;x+=0.05) // Neutron spectrum extends to higher E. - This computes the secondary spectrum from photons and protons. pi+, pi-, K+ and neutrons are considered as resulting particles. The input spectra are in [GeV^-1 cm^-3], while the resulting spectra of secondary particles are in [GeV^-1 cm^-3 s^-1].
	{
		ex=pow(10,x);
		respp=ncoComputeSecondarySpectrum(PIONIT,NCO_PROTON,NCO_PI_PLUS,grb_protons,grb_photons,ex); // pi+ production
		respm=ncoComputeSecondarySpectrum(PIONIT,NCO_PROTON,NCO_PI_MINUS,grb_protons,grb_photons,ex); // pi- production
		resk=ncoComputeSecondarySpectrum(NCO_K_PLUS_PROD,NCO_PROTON,NCO_K_PLUS,grb_protons,grb_photons,ex); // K+ production
		resn=ncoComputeSecondarySpectrum(PIONIT,NCO_PROTON,NCO_NEUTRON,grb_protons,grb_photons,ex); // neutron production
		resp=grb_protons(ex)*3.0e+5/shellthickness;  // Proton trivial escape
		
		fun.params=(void *)(&ex);
		ppres= gsl_integration_qag(&funpp,ex/EPMAX,1.0,0.0,RELERROR,5000,GSL_INTEG_GAUSS41,WSppint,&ppintegration,&pperrorest);
		if(ppres!=0)
		{
			printf("Int. error %i occured in calculation of pp interaction spectrum.\n",ppres);
			ppintegration = 0.0;
		}
		
  		respfrompp=ncold*3e10*ppintegration;
  		printf("Energy = %g,   Value of pions from pp = %g\n",ex,respfrompp);
		
		xa[n]=x;
		
		if(respp>1e-200) ya[n]=log10(respp); else ya[n]=-200.0;
		if(respm>1e-200) yb[n]=log10(respm); else yb[n]=-200.0;
		if(resk>1e-200) yc[n]=log10(resk); else yc[n]=-200.0;
		if(resn>1e-200) yd[n]=log10(resn); else yd[n]=-200.0;
		if(resp>1e-200) ypp[n]=log10(resp); else ypp[n]=-200.0;
		if(respfrompp>1e-200) ypfrompp[n]=log10(respfrompp); else ypfrompp[n]=-200.0;
		n++;
		// mioAddToOutput5(x,respp*ex*ex,respm*ex*ex,resk*ex*ex,resn*ex); // The spectra multiplied with the energy^2 are in [GeV cm^-3 s^-1].
	}
	// mioCloseOutput();
        printf("Escape rate: %g s^-1\n",3.0e+5/shellthickness);
	
	gsl_integration_workspace_free(WSppint);


// Interpolating functions for the injection spectra (from the secondary spectrum).
	ncoIPAlloc(&s_piplus_in,xa,ya,n); // PI_PLUS
	double w_piplus_in(double energy)
	{
		if (log10(energy) < s_piplus_in.min || log10(energy) > s_piplus_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_piplus_in,log10(energy)));}
	}

	ncoIPAlloc(&s_piminus_in,xa,yb,n); // PI_MINUS
	double w_piminus_in(double energy)
	{
		if (log10(energy) < s_piminus_in.min || log10(energy) > s_piminus_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_piminus_in,log10(energy)));}
	}

	ncoIPAlloc(&s_kplus_in,xa,yc,n); // K_PLUS
	double w_kplus_in(double energy)
	{
		if (log10(energy) < s_kplus_in.min || log10(energy) > s_kplus_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_kplus_in,log10(energy)));}
	}

	ncoIPAlloc(&s_neutron_in,xa,yd,n); // NEUTRON
	double w_neutron_in(double energy)
	{
		if (log10(energy) < s_neutron_in.min || log10(energy) > s_neutron_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_neutron_in,log10(energy)));}
	}

	ncoIPAlloc(&s_proton_in,xa,ypp,n); // PROTONS ESCAPE
	double w_proton_in(double energy)
	{
		if (log10(energy) < s_proton_in.min || log10(energy) > s_proton_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_proton_in,log10(energy)));}
	}

	 
	// The steady state computation for pions, kplus and neutrons. The input spectra for the steady state calculation are in [GeV^-1 cm^-3 s^-1] and the resulting "steady" spectra are in [GeV^-1 cm^-3].
	int nn,nmax=500000;
	double xx[nmax],yy[nmax];

	particlemass=0.140;
	particlecharge=1.0;
	particle=NCO_PI_PLUS;  
	// nn=ncoSteadyState2(w_piplus_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_piplus_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_piplus_steady,xx,yy,nn);

	particlemass=0.140;
	particlecharge=-1.0;
	particle=NCO_PI_MINUS;
	// nn=ncoSteadyState2(w_piminus_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_piminus_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_piminus_steady,xx,yy,nn);

	particlemass=0.494;
	particlecharge=1.0;
	particle=NCO_K_PLUS;
	// nn=ncoSteadyState2(w_kplus_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_kplus_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_kplus_steady,xx,yy,nn);

	particlemass=0.939;
	particlecharge=0.0;
	particle=NCO_NEUTRON;
	// nn=ncoSteadyState2(w_neutron_in,Eloss,Eescape,xx,yy,nmax);
	// ncoIPAlloc(&s_neutron_steady,xx,yy,nn);
	i=0; // Changed calculation to reflect change in nco_model1.
	double ressn;
	for(x=log10(EPMIN);x<=15.01;x+=0.05)
	{
		xa[i]=x;
		ex=pow(10,x);
		// respmr=w_neutron_in(ex)/(ncoComputeDecayEscapeRate(NCO_NEUTRON,ex)+3e5/R); //considers escape of the region and decays only neutrons in region
		ressn=w_neutron_in(ex)/ncoComputeDecayEscapeRate(NCO_NEUTRON,ex); //all neutrons decay do not distinct if in or outside the region
		if(ressn>1e-200) yd[i]=log10(ressn); else yd[i]=-200.0;
		i++;
	}
	ncoIPAlloc(&s_neutron_steady,xa,yd,i);

	 
	double piplus(double energy)
	{
		if ((log10(energy) < s_piplus_steady.min) || (log10(energy) > s_piplus_steady.max)) {return 0;}
			else {return pow(10,ncoIP(&s_piplus_steady,log10(energy)));}
	}

	double piminus(double energy)
	{
		if ((log10(energy) < s_piminus_steady.min) || (log10(energy) > s_piminus_steady.max)) {return 0;}
			else {return pow(10,ncoIP(&s_piminus_steady,log10(energy)));}
	}

	double kplus(double energy)
	{
		if (log10(energy) < s_kplus_steady.min || log10(energy) > s_kplus_steady.max) return 0;
			else return pow(10,ncoIP(&s_kplus_steady,log10(energy)));
	}

	double neutron(double energy)
	{
		if (log10(energy) < s_neutron_steady.min || log10(energy) > s_neutron_steady.max) return 0;
			else return pow(10,ncoIP(&s_neutron_steady,log10(energy)));
	}

	// Decay spectra of neutrons, pions and kaons. The input spectra for the particle decay function need to be in [GeV^-1 cm^-3] while the resulting spectra are in [GeV^-1 cm^-3 s^-1].
	double xb[201];
	double ye[201];
	double yf[201];
	double yg[201];
	double yh[201];
	double nbe1[201];
	double nmu1[201];
	double nbmu2[201];
	double nmu3[201];
	// mioInitOutput("kaondata.dat");
	double resnne,respmr,respml,respn,resmpmr,resmpml,resmpn,reskm;
	int w=0;
	for(x=0;x<=10.01;x+=0.05)
	{
		ex=pow(10,x);
		resnne=ncoComputeParticleDecaySpectrum(neutron,NCO_NEUTRON,NCO_NU_BAR_E,ex); // First set of (anti-)electron neutrinos (from neutron decay) -> nbe1.
		respmr=ncoComputeParticleDecaySpectrum(piplus,NCO_PI_PLUS,NCO_MU_PLUS_R,ex);
		respml=ncoComputeParticleDecaySpectrum(piplus,NCO_PI_PLUS,NCO_MU_PLUS_L,ex);
		respn=ncoComputeParticleDecaySpectrum(piplus,NCO_PI_PLUS,NCO_NU_MU,ex); // First set of muon neutrinos (from piâº decay) -> nmu1.
		resmpmr=ncoComputeParticleDecaySpectrum(piminus,NCO_PI_MINUS,NCO_MU_MINUS_R,ex);
		resmpml=ncoComputeParticleDecaySpectrum(piminus,NCO_PI_MINUS,NCO_MU_MINUS_L,ex);
		resmpn=ncoComputeParticleDecaySpectrum(piminus,NCO_PI_MINUS,NCO_NU_BAR_MU,ex); // Second set of (anti-)muon neutrinos (from piâ» decay ) -> nbmu2.
		reskm=ncoComputeParticleDecaySpectrum(kplus,NCO_K_PLUS,NCO_NU_MU,ex); // Third set of muon neutrinos (from Kâº decay) -> nmu3.
		// mioAddToOutput9(x,resnne*ex*ex,respmr*ex*ex,respml*ex*ex,respn*ex*ex,resmpmr*ex*ex,resmpml*ex*ex,resmpn*ex*ex,reskm*ex*ex); // The spectra multiplied with the energy^2 are in [GeV cm^-3 s^-1].
		xb[w]=x;
		if(respmr>1e-200) ye[w]=log10(respmr); else ye[w]=-200.0;
		if(respml>1e-200) yf[w]=log10(respml); else yf[w]=-200.0;
		if(resmpmr>1e-200) yg[w]=log10(resmpmr); else yg[w]=-200.0;
		if(resmpml>1e-200) yh[w]=log10(resmpml); else yh[w]=-200.0;
		if(resnne>1e-200) nbe1[w]=log10(resnne); else nbe1[w]=-200.0;
		if(respn>1e-200) nmu1[w]=log10(respn); else nmu1[w]=-200.0;
		if(resmpn>1e-200) nbmu2[w]=log10(resmpn); else nbmu2[w]=-200.0;
		if(reskm>1e-200) nmu3[w]=log10(reskm); else nmu3[w]=-200.0;
		w++;
	}
	// mioCloseOutput();


	// Interpolating functions for the injection spectra of the muons.
	ncoIPAlloc(&s_muplusr_in,xb,ye,w); // MU_PLUS_R
	double w_muplusr_in(double energy)
	{
		if (log10(energy) < s_muplusr_in.min || log10(energy) > s_muplusr_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_muplusr_in,log10(energy)));}
	}

	ncoIPAlloc(&s_muplusl_in,xb,yf,w); // MU_PLUS_L
	double w_muplusl_in(double energy)
	{
		if (log10(energy) < s_muplusl_in.min || log10(energy) > s_muplusl_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_muplusl_in,log10(energy)));}
	}

	ncoIPAlloc(&s_muminusr_in,xb,yg,w); // MU_MINUS_R
	double w_muminusr_in(double energy)
	{
		if (log10(energy) < s_muminusr_in.min || log10(energy) > s_muminusr_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_muminusr_in,log10(energy)));}
	}

	ncoIPAlloc(&s_muminusl_in,xb,yh,w); // M_MINUS_L
	double w_muminusl_in(double energy)
	{
		if (log10(energy) < s_muminusl_in.min || log10(energy) > s_muminusl_in.max) {return 0;}
			else {return pow(10,ncoIP(&s_muminusl_in,log10(energy)));}
	}
	 
	 
	//steady state for muons - input in [GeV^-1 cm^-3 s^-1], output in [GeV^-1 cm^-3].
	particlemass=0.105;
	particlecharge=1.0;
	particle=NCO_MU_PLUS_R;  
	// nn=ncoSteadyState2(w_muplusr_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_muplusr_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_muplusr_steady,xx,yy,nn);

	particlemass=0.105;
	particlecharge=1.0;
	particle=NCO_MU_PLUS_L;  
	// nn=ncoSteadyState2(w_muplusl_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_muplusl_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_muplusl_steady,xx,yy,nn);

	particlemass=0.105;
	particlecharge=1.0;
	particle=NCO_MU_MINUS_R;  
	// nn=ncoSteadyState2(w_muminusr_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_muminusr_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_muminusr_steady,xx,yy,nn);

	particlemass=0.105;
	particlecharge=1.0;
	particle=NCO_MU_MINUS_L;  
	// nn=ncoSteadyState2(w_muminusl_in,Eloss,Eescape,xx,yy,nmax);
	nn=ncoSteadyState4(w_muminusl_in,Eloss,Eescape,xx,yy,nmax,EPMAX,particle);
	ncoIPAlloc(&s_muminusl_steady,xx,yy,nn);

	 
	double muplusr(double energy)
	{
		if (log10(energy) < s_muplusr_steady.min || log10(energy) > s_muplusr_steady.max) {return 0;}
			else {return pow(10,ncoIP(&s_muplusr_steady,log10(energy)));}
	}

	double muplusl(double energy)
	{
		if (log10(energy) < s_muplusl_steady.min || log10(energy) > s_muplusl_steady.max) {return 0;}
			else {return pow(10,ncoIP(&s_muplusl_steady,log10(energy)));}
	}

	double muminusr(double energy)
	{
		if (log10(energy) < s_muminusr_steady.min || log10(energy) > s_muminusr_steady.max) {return 0;}
			else {return pow(10,ncoIP(&s_muminusr_steady,log10(energy)));}
	}


	double muminusl(double energy)
	{
		if (log10(energy) < s_muminusl_steady.min || log10(energy) > s_muminusl_steady.max) {return 0;}
			else {return pow(10,ncoIP(&s_muminusl_steady,log10(energy)));}
	}


	// Decay of muons - input in [GeV^-1 cm^-3], output in [GeV^-1 cm^-3 s^-1].
	double xc[201];
	double nbmu4[201];
	double ne2[201];
	double nbmu5[201];
	double ne3[201];
	double nmu6[201];
	double nbe4[201];
	double nmu7[201];
	double nbe5[201];
	// mioInitOutput("muondata.dat");
	double resmlbm,resmle,resmrbm,resmre,resmlm,resmlbe,resmrm,resmrbe,steadymuminusl;
	int v=0;
	for(x=0.00;x<=10.01;x+=0.05)
	{
		ex=pow(10,x);
		resmlbm=ncoComputeParticleDecaySpectrum(muplusl,NCO_MU_PLUS_L,NCO_NU_BAR_MU,ex); // Fourth set of (anti-)muon neutrinos -> nbmu4.
		resmle=ncoComputeParticleDecaySpectrum(muplusl,NCO_MU_PLUS_L,NCO_NU_E,ex); // Second set of electron neutrinos -> ne2.
		resmrbm=ncoComputeParticleDecaySpectrum(muplusr,NCO_MU_PLUS_R,NCO_NU_BAR_MU,ex); // Fifth set of (anti-)muon neutrinos -> nbmu5.
		resmre=ncoComputeParticleDecaySpectrum(muplusr,NCO_MU_PLUS_R,NCO_NU_E,ex); // Third set of electron neutrinos -> ne3.
		resmlm=ncoComputeParticleDecaySpectrum(muminusl,NCO_MU_MINUS_L,NCO_NU_MU,ex); // Sixth set of muon neutrinos -> nmu6.
		resmlbe=ncoComputeParticleDecaySpectrum(muminusl,NCO_MU_MINUS_L,NCO_NU_BAR_E,ex); // Fourth set of (anti-)electron neutrinos -> nbe4.
		resmrm=ncoComputeParticleDecaySpectrum(muminusr,NCO_MU_MINUS_R,NCO_NU_MU,ex); // Seventh set of muon neutrinos -> nmu7.
		resmrbe=ncoComputeParticleDecaySpectrum(muminusr,NCO_MU_MINUS_R,NCO_NU_BAR_E,ex); // Fifth set of (anti-)electron neutrinos -> nbe5.
		steadymuminusl=muminusr(ex);
		// mioAddToOutput10(x,resmlbm*ex*ex,resmle*ex*ex,resmrbm*ex*ex,resmre*ex*ex,resmlm*ex*ex,resmlbe*ex*ex,resmrm*ex*ex,resmrbe*ex*ex,steadymuminusl*ex*ex); // The spectra multiplied by the energy^2 are in [GeV cm^-3 s^-1].
		xc[v]=x;
		if(resmlbm>1e-200) nbmu4[v]=log10(resmlbm); else nbmu4[v]=-200.0;
		if(resmle>1e-200) ne2[v]=log10(resmle); else ne2[v]=-200.0;
		if(resmrbm>1e-200) nbmu5[v]=log10(resmrbm); else nbmu5[v]=-200.0;
		if(resmre>1e-200) ne3[v]=log10(resmre); else ne3[v]=-200.0;
		if(resmlm>1e-200) nmu6[v]=log10(resmlm); else nmu6[v]=-200.0;
		if(resmlbe>1e-200) nbe4[v]=log10(resmlbe); else nbe4[v]=-200.0;
		if(resmrm>1e-200) nmu7[v]=log10(resmrm); else nmu7[v]=-200.0;
		if(resmrbe>1e-200) nbe5[v]=log10(resmrbe); else nbe5[v]=-200.0;
		v++;
	}
	// mioCloseOutput();
	 
	 
	double nbetot[201]; // Summing up the different contributions to a neutrino flavor. All functions in [GeV^-1 cm^-3 s^-1].
	double netot[201];
	double nbmutot[201];
	double nmutot[201];
	// mioInitOutput("prodneutrinos.dat");
	int u=0;
	for(x=0.0;x<=10.01;x+=0.05)
	{
		nbetot[u]=log10(pow(10,nbe1[u])+pow(10,nbe4[u])+pow(10,nbe5[u]));
		netot[u]=log10(pow(10,ne2[u])+pow(10,ne3[u]));
		nbmutot[u]=log10(pow(10,nbmu2[u])+pow(10,nbmu4[u])+pow(10,nbmu5[u]));
		nmutot[u]=log10(pow(10,nmu1[u])+pow(10,nmu3[u])+pow(10,nmu6[u])+pow(10,nmu7[u]));
		// mioAddToOutput5(x,pow(10,nbetot[u])*pow(10.0,2.0*x),pow(10,netot[u])*pow(10.0,2.0*x),pow(10,nbmutot[u])*pow(10.0,2.0*x),pow(10,nmutot[u])*pow(10.0,2.0*x));
		u++;
	}
	// mioCloseOutput();

	nco_ip ipnbe; // Interpolating function for NU_BAR_E.
	ncoIPAlloc(&ipnbe,xc,nbetot,u);
	double neubaretot(double energy)
	{
		if (log10(energy) < ipnbe.min || log10(energy) > ipnbe.max) {return 0;}
			else {return pow(10,ncoIP(&ipnbe,log10(energy)));}
	}

	nco_ip ipne; // Interpolating funtion for NU_E.
	ncoIPAlloc(&ipne,xc,netot,u);
	double neuetot(double energy)
	{
		if (log10(energy) < ipne.min || log10(energy) > ipne.max) {return 0;}
			else {return pow(10,ncoIP(&ipne,log10(energy)));}
	}

	nco_ip ipnbmu; // Interpolating function for NU_BAR_MU.
	ncoIPAlloc(&ipnbmu,xc,nbmutot,u);
	double neubarmutot(double energy)
	{
		if (log10(energy) < ipnbmu.min || log10(energy) > ipnbmu.max) {return 0;}
			else {return pow(10,ncoIP(&ipnbmu,log10(energy)));}
	}

	nco_ip ipnmu; // Interpolating function for NU_MU.
	ncoIPAlloc(&ipnmu,xc,nmutot,u);
	double neumutot(double energy)
	{
		if (log10(energy) < ipnmu.min || log10(energy) > ipnmu.max) {return 0;}
			else {return pow(10,ncoIP(&ipnmu,log10(energy)));}
	}

	// Computation of the effect of the cosmic expansion and the Lorentz boost. Input in [GeV^-1 cm^-3 s^-1] and output in [GeV^-1 cm^-2 s^-1] (ncoComputeGRBCosmicSpectrum() calculates particle flux of one burst).
	double xd[201];
	double boostnbe[201];
	double boostne[201];
	double boostnbmu[201];
	double boostnmu[201];
	double boostneu[201];
	double boostprot[201];
	// mioInitOutput("boostneutrinos.dat");
	double boonbe,boone,boonbmu,boonmu,booneu,booprot;
	int t=0;
	for(x=0;x<=10.01;x+=0.05)
	{
		ex=pow(10,x+log10(doppler(lorentz,viewingangle)/(1.0+redshift)));
		boonbe = ncoComputeGRBCosmicSpectrum(neubaretot,lorentz,viewingangle,redshift,grbradius,shellthickness,ex); // 5e3, (->shellthickness). Depending on function, either ...GRB..., or no "GRB" in name, one has to include the radius and the shell thickness or just the radius.
		boone = ncoComputeGRBCosmicSpectrum(neuetot,lorentz,viewingangle,redshift,grbradius,shellthickness,ex);
		boonbmu = ncoComputeGRBCosmicSpectrum(neubarmutot,lorentz,viewingangle,redshift,grbradius,shellthickness,ex);
		boonmu = ncoComputeGRBCosmicSpectrum(neumutot,lorentz,viewingangle,redshift,grbradius,shellthickness,ex);
		booneu = ncoComputeGRBCosmicSpectrum(w_neutron_in,lorentz,viewingangle,redshift,grbradius,shellthickness,ex);
		// Neutrons -> Protons
		booprot = ncoComputeGRBCosmicSpectrum(w_proton_in,lorentz,viewingangle,redshift,grbradius,shellthickness,ex);
		// Direct Proton escape; compute from injection spectrum!
		// mioAddToOutput5(x+log10(doppler(lorentz,viewingangle)/(1.0+redshift)),boonbe*ex*ex,boone*ex*ex,boonbmu*ex*ex,boonmu*ex*ex); // First column: energy; second: anti-electron neutrinos; thrid: electron neutrinos; fourth: anti-muon neutrinos; fifth: muon neutrinos. Spectra multiplied with the energy^2 are in [GeV cm^-2 s^-1].
		xd[t]=x+log10(doppler(lorentz,viewingangle)/(1.0+redshift));
		if(boonbe>1e-200) boostnbe[t]=log10(boonbe); else boostnbe[t]=-200.0;
		if(boone>1e-200) boostne[t]=log10(boone); else boostne[t]=-200.0;
		if(boonbmu>1e-200) boostnbmu[t]=log10(boonbmu); else boostnbmu[t]=-200.0;
		if(boonmu>1e-200) boostnmu[t]=log10(boonmu); else boostnmu[t]=-200.0;
		if(booneu>1e-200) boostneu[t]=log10(booneu); else boostneu[t]=-200.0;
		if(booprot>1e-200) boostprot[t]=log10(booprot); else boostprot[t]=-200.0;
		t++;
	}
	// mioCloseOutput();

	 
	nco_ip ipbnbe; // Interpolating function for boosted NU_BAR_E.
	ncoIPAlloc(&ipbnbe,xd,boostnbe,u);
	double cosmicneubare(double energy)
	{
		if (log10(energy) < ipbnbe.min || log10(energy) > ipbnbe.max) {return 0;}
			else {return pow(10,ncoIP(&ipbnbe,log10(energy)));}
	}

	nco_ip ipbne; // Interpolating funtion for boosted NU_E.
	ncoIPAlloc(&ipbne,xd,boostne,u);
	double cosmicneue(double energy)
	{
		if (log10(energy) < ipbne.min || log10(energy) > ipbne.max) {return 0;}
			else {return pow(10,ncoIP(&ipbne,log10(energy)));}
	}

	nco_ip ipbnbmu; // Interpolating function for boosted NU_BAR_MU.
	ncoIPAlloc(&ipbnbmu,xd,boostnbmu,u);
	double cosmicneubarmu(double energy)
	{
		if (log10(energy) < ipbnbmu.min || log10(energy) > ipbnbmu.max) {return 0;}
			else {return pow(10,ncoIP(&ipbnbmu,log10(energy)));}
	}

	nco_ip ipbnmu; // Interpolating function for boosted NU_MU.
	ncoIPAlloc(&ipbnmu,xd,boostnmu,u);
	double cosmicneumu(double energy)
	{
		if (log10(energy) < ipbnmu.min || log10(energy) > ipbnmu.max) {return 0;}
			else {return pow(10,ncoIP(&ipbnmu,log10(energy)));}
	}

	nco_ip ipbneu; // Interpolating function for boosted CR neutrons = protons
	ncoIPAlloc(&ipbneu,xd,boostneu,u);
	double cosmicneutr(double energy)
	{
		if (log10(energy) < ipbneu.min || log10(energy) > ipbneu.max) {return 0;}
			else {return pow(10,ncoIP(&ipbneu,log10(energy)));}
	}

		nco_ip ipbprot; // Interpolating function for boosted protons (direct escape)
	ncoIPAlloc(&ipbprot,xd,boostprot,u);
	double cosmicprot(double energy)
	{
		if (log10(energy) < ipbprot.min || log10(energy) > ipbprot.max) {return 0;}
			else {return pow(10,ncoIP(&ipbprot,log10(energy)));}
	}

	// Computation of the oscillation effects after the taking other effects into account. All spectra, input and output, are in [GeV^-1 cm^-2 s^-1].
	double xe[201];
	double oscnbe[201];
	double oscne[201];
	double oscnbmu[201];
	double oscnmu[201];
	double oscnbtau[201];
	double oscntau[201];
	mioInitOutput("Test_oscneutrinos.dat"); // This line sets a filename for the neutrino spectra of one burst. It needs to be changed in each run or the old file will be overwritten.
	double onbe,one,onbmu,onmu,onbtau,ontau;
	 
	double theta12 = 0.599; // Oscillation angles for the three neutrino flavors. (Values from Feb. 2010!)
	double theta13 = 0.0;
	double theta23 = M_PI/4.0;
	 
	int r=0;
	for(x=0;x<=10.01;x+=0.05)
	{
		ex=pow(10,x+log10(doppler(lorentz,viewingangle)/(1.0+redshift)));
		onbe = ncoComputeFlavorSpectrum(NCO_NU_BAR_E,cosmicneubare,cosmicneubarmu,ex,theta12,theta13,theta23,0.0);
		one = ncoComputeFlavorSpectrum(NCO_NU_E,cosmicneue,cosmicneumu,ex,theta12,theta13,theta23,0.0);
		onbmu = ncoComputeFlavorSpectrum(NCO_NU_BAR_MU,cosmicneubare,cosmicneubarmu,ex,theta12,theta13,theta23,0.0);
		onmu = ncoComputeFlavorSpectrum(NCO_NU_MU,cosmicneue,cosmicneumu,ex,theta12,theta13,theta23,0.0);
		onbtau = ncoComputeFlavorSpectrum(NCO_NU_BAR_TAU,cosmicneubare,cosmicneubarmu,ex,theta12,theta13,theta23,0.0);
		ontau = ncoComputeFlavorSpectrum(NCO_NU_TAU,cosmicneue,cosmicneumu,ex,theta12,theta13,theta23,0.0);
		mioAddToOutput9(x+log10(doppler(lorentz,viewingangle)/(1.0+redshift)),Cp*Cgamma*onbe*ex*ex,Cp*Cgamma*one*ex*ex,Cp*Cgamma*onbmu*ex*ex,Cp*Cgamma*onmu*ex*ex,
		Cp*Cgamma*cosmicneutr(ex)*ex*ex,Cp*cosmicprot(ex)*ex*ex,
		Cp*Cgamma*onbtau*ex*ex,Cp*Cgamma*ontau*ex*ex);
		// First column: energy; second: anti-electron neutrinos; thrid: electron neutrinos; fourth: anti-muon neutrinos; fifth: muon neutrinos; sixth: neutron=proton escape; seventh: direct escape (t_esc = delta d') - Spectra multiplied by the energy^2 are in [GeV cm^-2 s^-1].
		xe[r]=x+log10(doppler(lorentz,viewingangle)/(1.0+redshift));
		if(onbe>1e-200) oscnbe[r]=log10(onbe); else oscnbe[r]=-200.0;
		if(one>1e-200) oscne[r]=log10(one); else oscne[r]=-200.0;
		if(onbmu>1e-200) oscnbmu[r]=log10(onbmu); else oscnbmu[r]=-200.0;
		if(onmu>1e-200) oscnmu[r]=log10(onmu); else oscnmu[r]=-200.0;
		if(onbtau>1e-200) oscnbtau[r]=log10(onbtau); else oscnbtau[r]=-200.0;
		if(ontau>1e-200) oscntau[r]=log10(ontau); else oscntau[r]=-200.0;
		r++;
	}
	mioCloseOutput();


	nco_ip iposcnbmu; // Interpolating function for oscillated NU_BAR_MU.
	ncoIPAlloc(&iposcnbmu,xe,oscnbmu,r);
	double oscillatedneubarmu(double energy)
	{
		if (log10(energy) < iposcnbmu.min || log10(energy) > iposcnbmu.max) {return 0;}
			else {return pow(10,ncoIP(&iposcnbmu,log10(energy)));}
	}

	nco_ip iposcnmu; // Interpolating function for oscillated NU_MU.
	ncoIPAlloc(&iposcnmu,xe,oscnmu,r);
	double oscillatedneumu(double energy)
	{
		if (log10(energy) < iposcnmu.min || log10(energy) > iposcnmu.max) {return 0;}
			else {return pow(10,ncoIP(&iposcnmu,log10(energy)));}
	}


	ncoIPFree(&ipnbe);
	ncoIPFree(&ipne);
	ncoIPFree(&ipnbmu);
	ncoIPFree(&ipnmu);
	ncoIPFree(&ipbnbe);
	ncoIPFree(&ipbne);
	ncoIPFree(&ipbneu);
	ncoIPFree(&ipbprot);
	ncoIPFree(&ipbnbmu);
	ncoIPFree(&ipbnmu);
	ncoIPFree(&iposcnbmu);
	ncoIPFree(&iposcnmu);
	// ncoIPFree(&ipAeff);
	// ncoIPFree(&ipsimpAeff);
	ncoDestroySpectra();
}


int main()
{
// Initialising the whole generators
ncoInit();


// Starting here all the properties of the GRB are calculated from the "basic" properties of the GRB. In the last part the normalizations of the photon and the proton spectrum are calculated.

double Viso = 16.0*M_PI*pow(29979245800.0,3.0)*pow(lorentz,5.0)*pow(1.0+redshift,-3.0)*pow(tvar,3.0); // Isotropic volume in the shock rest frame, in [cm^3].
double Eiso = luminosity*624.15*duration/(lorentz*(1.0+redshift)); // Isotropic energy in the shock rest frame, in [GeV].
double nos = duration/tvar; // Number Of Shells, in [1].

B = sqrt(8.0*M_PI*Eiso/(624.15*nos*Viso)); // Magnetic field in [G]. epsilonb/epsilone = 1
grbradius = 2.0*pow(lorentz,2.0)*29979245800.0*tvar/(1.0+redshift)*1e-5; // Collision radius of the shells of the GRB in [km]!!!
shellthickness = lorentz*29979245800.0*tvar/(1.0+redshift)*1e-5; // Shell thickness of the shells of the GRB in [km]!!!

EPHMIN = 1e-6*(1.0+redshift)/lorentz; // Minimal photon energy in [GeV]
EPHMAX = 1e-2*(1.0+redshift)/lorentz; // Maximal photon energy in [GeV]

//maxprotonenergy = 2.00982e11*sqrt(accelerationefficiency/B); // Synchrotron limited case - Magnetic field in [G], accelerationefficiency in [1], and resulting maxprotonenergy in [GeV].
//maxprotonenergy = 1e-9*(lorentz*tvar/(1.0+redshift)*3e8)*accelerationefficiency*3e8*(B*1e-4); // Adiabatic limited case - Magnetic field in [G], accelerationefficiency in [1], Lorentz factor in [1], redshift in [1], variability time scale in [s] and resulting maxprotonenergy in [GeV].
//Structure of formula: (eV->GeV)*(size of region in [m])*(acceleration efficiency)*(magnetic field in [T])

//********* Automatic version to choose synchrotron or adiabatic limited case ********************************************************

	double synchlimit = 2.00982e11*sqrt(accelerationefficiency/B);
	double esclimit = 1e-9*(lorentz*tvar/(1.0+redshift)*3e8)*accelerationefficiency*3e8*(B*1e-4); // (eV->GeV)*(size of region in [m])*(acceleration efficiency)*(magnetic field in [T])

	if(synchlimit<esclimit) maxprotonenergy = synchlimit; 
		else maxprotonenergy = esclimit; 
//************************************************************************************************************************************


double FGAMMA(double y, void* params)
{
return y*grb_photons(y);
}

double FP(double y, void* params)
{
return y*grb_protons(y);
}

const double RELERROR=1e-3;
double spectraerrorest;
int res, res2;
gsl_function fun,fun2;
fun.function=&FGAMMA;
fun.params=NULL;
fun2.function=&FP;
fun2.params=NULL;

gsl_integration_workspace* WSspectra=gsl_integration_workspace_alloc(5000);

res= gsl_integration_qag(&fun,1e-6*(1.0+redshift)/lorentz,1e-2*(1.0+redshift)/lorentz,0.0,RELERROR,5000,GSL_INTEG_GAUSS41,WSspectra,&intenergygamma,&spectraerrorest);
if(res!=0)
{
printf("Int. error %i occured in normalization of photon spectrum\n",res);
intenergygamma = 0.0;
}

res2= gsl_integration_qag(&fun2,1.0,maxprotonenergy,0.0,RELERROR,5000,GSL_INTEG_GAUSS41,WSspectra,&intenergyp,&spectraerrorest);
if(res2!=0)
{
printf("Int. error %i occured in normalization of proton spectrum\n",res2);
intenergyp = 0.0;
}

gsl_integration_workspace_free(WSspectra);

Cgamma = Eiso/(nos*Viso*intenergygamma);
Cp = Eiso/(fe*nos*Viso*intenergyp);

double fluence = duration*luminosity/(4.0*M_PI*luminositydistance(redshift)*luminositydistance(redshift));
printf("Observed fluence: %e\n",fluence);


// In this step finally the spectra are calculated by calling the function ComputeNeutrinoSpectraFromScratch. The additional parameters are needed if a larger number of data sets should be created and one wants to have automatic files for the neutrinos. When wanted this also has to be included in the definition of the function.

ComputeNeutrinoSpectraFromScratch(B,lorentz,viewingangle,redshift,Cp,Cgamma,grbradius,shellthickness); //,oscneutrinos,prodneutrinos,boostneutrinos,piondata,kaondata,muondata);


 
ncoQuit();
exit(0);
}
