#include <math.h>
#include <float.h>
#include <stdlib.h> 
#include <stdio.h>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

int cvode_check_flag (void *flagvalue, const char *funcname, int opt);
double sqrt1px_m1 (const double x);

#define LOG_DBL_MAX 7.0978271289338397e+02

#define CVODE_CHECK(chk,name,val,ret) \
do { \
  if (!cvode_check_flag (chk, name, val)) \
    return ret; \
} while (FALSE)

struct _ExpoConsts
{
  double Rc;
  double lpl;
  double aii;
  double phiii;
  double V0;
  double kappa;
  double ticl;
  double sigma;
  double d;  
  double kpert;
  double H_scale;
  double alpha0;
  double alpha_B;

  
};

typedef struct _ExpoConsts ExpoConsts;

struct _SyState
{
	double alpha;
	double phi;
	double mu;
	double x;
	double t;
	double alpha_B; 
	long double H;
	double V0;
	double zeta_re;
	double dzeta_re;
	double zeta_im;
	double dzeta_im;
	double zeta_ns;
	double pok_ns;
	double tau_ns;
	double w_ns;
	double lab1;
	double lab2;

}; 

typedef struct _SyState SyState;

/*************************************************************************************
 * 
 * Pre-declarations
 * 
 *************************************************************************************/

static double Expo_dmu_p (const double alpha, const double mu_p);
static double Expo_dmu_m (const double alpha, const double mu_p);
static long  double Expo_xqt(const double alpha, const double phi, ExpoConsts*c1);
static long double Expo_dalpha_dt_qt(const double t, const double alpha, const double phi, ExpoConsts*c1);
/*************************************************************************************
 * 
 * Perturbations
 * 
 *************************************************************************************/

/*z''/z == Expo_V*/
static double
Expo_V (const double alpha, const double mu, const long double H)
{
	double x_aux;
	if (H < 0)
	{
		x_aux	= 1.0 - mu;		
	}
	else
	{
		x_aux	= mu - 1.0;;
	}
	
	const double E_a	= exp(alpha);
	const double E2_a	= E_a * E_a ;
	const double H2		= H * H;
	const double x2		= x_aux * x_aux;
	const double x3		= x2 * x_aux;
	const double x4		= x2 * x2;
	
	//~ printf(" E2_a %20.15e \n H2 %20.15e  \n", E2_a, H2);
	
	const double V 		= (E2_a * H2) * (-7.0 + 18.0 * sqrt(2.0) * x_aux- 12.0 * x2 - 18.0 * sqrt(2.0) * x3 + 18.0 * x4);
	
	return V;
}

/*d(z''/z) / d n == Expo_Vl*/
static double 
Expo_Vl(const double alpha, const double mu, const long double H)
{
	double x_aux;
	if(H<0)
	{
		x_aux	= 1.0 - mu;		
	}
	else
	{
		x_aux	= mu - 1.0;;
	}
	
	const double E_a	= exp(alpha);
	const double E2_a	= E_a * E_a ;
	const double H2		= H * H;
	const double x2		= x_aux * x_aux;
	const double x3		= x2 * x_aux;
	const double x4		= x2 * x2;
	const double x5		= x4 * x_aux;
	const double x6		= x4 * x2;
	
	const double Vl		= E2_a * H2 * (40.0 - 54.0 * sqrt(2.0) * x_aux - 126.0 * x2 + 216 * sqrt(2.0)* x3 - 18.0 * x4 - 162.0 * sqrt(2.0) * x5 + 108.0 * x6);
	
	return Vl;
} 

/* |z''/z / k2| == Expo_epsilon*/
static double 
Expo_epsilon (const double alpha, const double mu, const long double H, ExpoConsts*c1 )
{	
	const double V			= Expo_V(alpha, mu, H);
	const double num		= V;
	const double den 		= (c1->kpert * exp (c1->alpha0)) * (c1->kpert * exp (c1->alpha0));
	const double epsilon	= fabs (num / den);
	
	//~ printf(" V %20.15e num %20.15e den %20.15e", V, num, den);
	
	return epsilon;	
}

/* vk_WKB == Expo_vk, pt. Re e Im*/
static double 
Expo_vk (const double alpha, const double mu, const long double H, ExpoConsts*c1 )
{
	const double V		= Expo_V (alpha, mu, H);
	const double k		= c1->kpert * exp (c1->alpha0);
	const double k2		= k * k;
	const double w		= sqrt (k2 - V);
	const double fase	= sqrt (2.0) / 2.0;	
	const double num	= fase;
	const double den	= sqrt (2.0 * w);

	const double vk	= num / den;

	return vk;
}


/* vkl_WKB == Expo_vkl_re, pt. Re*/
static double 
Expo_vkl_re(const double alpha, const double mu, const long double H,  ExpoConsts*c1 )
{
	const double V		= Expo_V (alpha, mu, H);
	const double Vl		= Expo_Vl (alpha, mu, H);
	const double k		= c1->kpert * exp (c1->alpha0);
	const double k2		= k * k;
	const double w		= sqrt (k2 - V);
	const double w2		= k2 - V;
	const double w3		= w2 * w;
	const double w5		= w2 * w3;
	const double sq_w5	= sqrt (w5);
	
	const double vkl_re	= - 0.5 * sqrt(w)  + Vl / (8.0 * sq_w5);
	//~ const double den	= 8.0 * sq_w5;
	
	//~ const double num	= - 4.0 * w3 + Vl;
	//~ const double den	= 8.0 * sq_w5;
	
	//~ printf("w3 %20.15g Vl %20.15g raz %20.15g V %20.15g \n", 4 * w3, Vl, (4 * w3) / Vl, V );
	
	//~ const double vk_l	= num / den;

	return vkl_re;
}
//~ 
static double 
Expo_vkl_im (const double alpha, const double mu, const long double H,  ExpoConsts*c1 )
{
	const double V		= Expo_V (alpha, mu, H);
	const double Vl		= Expo_Vl (alpha, mu, H);
	const double k		= c1->kpert * exp (c1->alpha0);
	const double k2		= k * k;
	const double w		= sqrt (k2 - V);
	const double w2		= k2 - V;
	const double w3		= w2 * w;
	const double w5		= w2 * w3;
	const double sq_w5	= sqrt(w5);
	
	//~ const double num	= 4.0 * w3 + Vl;
	//~ const double den	= 2.0 * sq_w5;
	//~ const double vk_l_im	= num / den;
	
	const double vkl_im	=  0.5 * sqrt(w)  + Vl / (8.0 * sq_w5);
	
	
	return vkl_im;
}

static double
Expo_z (const double alpha, const double mu, const long double H, ExpoConsts*c1)
{
	double x;
	if (H < 0)
	{
		x	= 1.0 - mu;		
	}
	else
	{
		x	= mu - 1.0;
	}

	//~ const double kden	= 1.0 / c1->kappa;
	const double z 			= exp (alpha) * x;

	return z;
}

static double 
Expo_zl_z (const double alpha, const double mu, const long double H)
{
	const double dmu_p = Expo_dmu_p(alpha, mu);
	const double dmu_m = Expo_dmu_m(alpha, mu);
	double zl_z;
	
	if(H<0)
	{
		zl_z	= (exp(alpha) * H * (1.0 - mu - dmu_m)) / (1.0 - mu);		
	}
	else
	{
		zl_z	= (exp(alpha) * H * (mu - 1.0 + dmu_p)) / (mu - 1.0);
	}
	
	return zl_z;
}

static long double 
Expo_zl_z_qt (const double t, const double alpha, const double phi, ExpoConsts*c1)
{
	//~ printf("zl_ze %20.15e %20.15e\n", alpha, phi);
	const long double phi_aux   = (long double)phi;
	const long double alpha_aux = (long double)alpha;
	const long double sigma_aux = (long double)c1->sigma;
		
	const long double sigma2    = sigma_aux * sigma_aux;
    const long double d_aux     = (long double)c1->d;
    const long double arg_tri   = 2.0*alpha_aux*d_aux ;
    const long double arg_hip   = sigma2 * alpha_aux * phi_aux;
    const long double cos_tri   = cosl(arg_tri); 
    const long double cos_tri4  = cosl( 2.0 * arg_tri );
    const long double cos_hip   = coshl(arg_hip);
    const long double sen_tri   = sinl(arg_tri);
    const long double sen_hip   = sinhl(arg_hip);
    
    const long double num1		= -sigma2 * (16.0 * d_aux * d_aux + sigma2) * phi_aux + sigma2 * sigma2 * cos_tri4 * phi - 
	8.0 * d_aux * d_aux * sigma2 * cos_tri * ( sen_hip * alpha_aux + 
    2.0 * cos_hip * phi_aux) + 
 4.0 * d_aux * sen_tri * (sigma2 * sigma2 * cos_hip * alpha_aux * phi_aux + 
    sen_hip * (-4.0 * d_aux * d_aux - sigma2 * sigma2 + sigma2 * sigma2 * phi_aux * phi_aux));
    
    const long double H			= c1->Rc * Expo_dalpha_dt_qt(t, alpha, phi, c1);
    
    const long double den1		= (2.0 * d_aux * sen_hip + sigma2 * sen_tri * phi_aux); 
	const long double den2		= 2.0 * den1 * den1;
	const long double xl		=  num1 / den2;
	
	//~ printf("zl_z2 %20.15Le %20.15Le\n", alpha_aux, phi_aux);
	
	const long double x			= Expo_xqt(alpha_aux, phi_aux, c1);
	const long double num		= exp(alpha_aux) * H * (x + xl);
	const long double den		= x;
	
	const long double zl_z_qt	= num / den;
	
	//~ printf("zl_zs %20.15Le %20.15Le zl_z %20.15Le\n", alpha_aux, phi_aux, zl_z_qt);
	
	//~ printf(">> alpha %20.15Le phi %20.15Le\n", alpha, phi);
	
	return zl_z_qt;	
}

static double
Expo_zeta (const double alpha, const double mu, const long double H, ExpoConsts*c1 )
{	
	const double vk		= Expo_vk (alpha, mu, H, c1 );
	const double z		= Expo_z (alpha, mu, H, c1);
	
	const double zeta	= vk / z;
	
	return zeta;
}

static double
Expo_dzeta_re (const double alpha, const double mu, const long double H,  ExpoConsts*c1 )
{
	const double vk_re			= Expo_vk (alpha, mu, H, c1);
	const double vkl_re			= Expo_vkl_re (alpha, mu, H,c1);
	const double z				= Expo_z (alpha, mu, H, c1);
	const double zl_z			= Expo_zl_z (alpha, mu, H);
	//~ const double alpha_H		= exp(alpha) * H;
	const double zeta_re		= vk_re / z;			
	
	const double dzetal_re	= (vkl_re / z  - zeta_re * zl_z);

	return dzetal_re;
}
//~ 
static double
Expo_dzeta_im (const double alpha, const double mu, const long double H,  ExpoConsts*c1 )
{
	const double vk_im			= Expo_vk (alpha, mu, H, c1);
	const double vkl_im			= Expo_vkl_im (alpha, mu, H,c1);
	const double z				= Expo_z (alpha, mu, H, c1);
	const double zl_z			= Expo_zl_z(alpha, mu, H);
	//~ const double alpha_H		= exp(alpha) * H;
	const double zeta_im		= vk_im / z;			
	
	const double dzeta_im	= ( vkl_im / z  - zeta_im * zl_z);
	
	return dzeta_im;
}

static double 
Expo_dzeta_dalpha (const double alpha, const double mu, const long double H, const double zeta, const double dzeta, ExpoConsts*c1 )
{
	const double alpha_H	= exp (alpha) * H;
	const double zl_z		= Expo_zl_z (alpha, mu, H);
	const double k			= c1->kpert * exp (c1->alpha0);
	const double k2			= k * k;
	
	const double dzeta_dalpha = (- 2.0 * zl_z * dzeta - k2 * zeta) / alpha_H;

	//~ printf("zl_z = %20.15g\n", zl_z);

	return dzeta_dalpha;
}

static double
Expo_dzeta_dt (const double t, const double alpha, const double phi, const double zeta, const double dzeta, ExpoConsts*c1 )
{
	//~ printf("dzdte %20.15e %20.15e %20.15e %20.15e \n", alpha, phi, zeta, dzeta);
	const long double alpha_aux	= (long double)alpha;
	const long double phi_aux	= (long double)phi;
	const long double exp_alpha	= exp (alpha_aux) ;
	const long double zeta_aux	= (long double)zeta;
	const long double dzeta_aux	= (long double)dzeta;
	//~ printf("dzdt2 %20.15e %20.15e\n", alpha_aux, phi_aux);
	const double zl_z_qt	= Expo_zl_z_qt (t, alpha_aux, phi_aux, c1);
	//~ printf("dzdt3 %20.15e %20.15e\n", alpha, phi);
	const double k			= c1->kpert * exp (c1->alpha0);
	const double k2			= k * k;
	//~ const double term1 		= (-2.0 * zl_z_qt * dzeta) / exp_alpha;
	//~ const double term2		= (- k2 * zeta) / exp_alpha;
	//~ const double term0		= 1.0 / exp_alpha;
	
	//~ printf("dzdts %20.15Le %20.15Le %20.15Le %20.15Le\n", alpha_aux, phi_aux, zeta_aux, dzeta_aux);
	
	const double dzeta_dt = (-2.0 * zl_z_qt * dzeta_aux - k2 * zeta_aux) / exp_alpha;
	
	//~ printf(" >>  alpha %20.15e phi %20.15e\n", alpha_aux, phi_aux);
	
	return dzeta_dt;
}

/*************************************************************************************
 * 
 * Background - quantum
 * 
 *************************************************************************************/

static long double
Expo_xqt (const double alpha, const double phi, ExpoConsts*c1)
{   
	//~ printf(" xqte %20.15e %20.15e\n", alpha, phi);
	const long double phi_aux   = (long double)phi;
	const long double alpha_aux = (long double)alpha;
	
	//~ printf("xqt2 %20.15e %20.15e\n", alpha_aux, phi_aux);
	const long double sigma_aux = (long double)c1->sigma;
	const long double sigma2    = sigma_aux * sigma_aux;
    const long double d_aux     = (long double)c1->d;
    const long double arg_tri   = 2.0*alpha_aux*d_aux ;
    const long double arg_hip   = sigma2 * alpha_aux * phi_aux;
    const long double cos_tri   = cosl(arg_tri);
    const long double cos_hip   = coshl(arg_hip);
    const long double sen_tri   = sinl(arg_tri);
    const long double sen_hip   = sinhl(arg_hip);
    const long double num       =  - alpha_aux * sigma2 * sen_tri + 2.0 * d_aux * cos_tri + 2.0 * d_aux * cos_hip;
    const long double den       = phi_aux * sigma2 * sen_tri + 2.0 * d_aux * sen_hip;
    
  
    const long double xqt       = num / den;
  	
	if (isinf (cos_hip) || isinf (sen_hip))
	{
		printf( "# arg_hip_INF = % 20.15Le alpha = % 20.15Le phi = % 20.15Le\n",arg_hip, alpha_aux, phi_aux );
	}
  
	//~ printf("xqts %20.15e %20.15e\n", alpha, phi);
    return xqt;
}

static long double 
Expo_dphi_dt_qt (const double t, const double alpha, const double phi, ExpoConsts*c1)
{  
	const long double phi_aux    = (long double)phi;
	const long double alpha_aux  = (long double)alpha;
	const long double sigma_aux  = (long double)c1->sigma;
	const long double sigma2     = sigma_aux * sigma_aux;
    const long double d_aux      = (long double)c1->d;
    const long double arg_tri    = 2.0 * alpha_aux * d_aux ;
    const long double arg_hip    = sigma2 * alpha_aux * phi_aux;
    const long double cos_tri    = cosl(arg_tri);
    const long double cos_hip    = coshl(arg_hip);
    const long double sen_tri    = sinl(arg_tri);
   	const long double num        =  - alpha_aux * sigma2 * sen_tri + 2.0 * d_aux * cos_tri + 2.0 * d_aux * cos_hip;
	const long double den        = exp(3.0 * alpha_aux) * 2.0 * (cos_tri + cos_hip);
	
	const long double dphi_dt_qt =  num/ den ;

	return dphi_dt_qt;	
}

static long double 
Expo_dalpha_dt_qt (const double t, const double alpha, const double phi, ExpoConsts*c1)
{
	const long double phi_aux      = (long double)phi;
	const long double alpha_aux    = (long double)alpha;
	const long double sigma_aux    = (long double)c1->sigma;
	const long double sigma2       = sigma_aux * sigma_aux;
    const long double d_aux        = (long double)c1->d;
    const long double arg_tri      = 2.0*alpha_aux*d_aux ;
    const long double arg_hip      = sigma2 * alpha_aux * phi_aux;
    const long double cos_tri      = cosl(arg_tri);
    const long double cos_hip      = coshl(arg_hip);
    const long double sen_tri      = sinl(arg_tri);
    const long double sen_hip      = sinhl(arg_hip);
	const long double num          = phi_aux * sigma2 * sen_tri + 2.0 * d_aux * sen_hip;
	const long double den          = exp(3.0 * alpha_aux) * 2.0 * (cos_tri + cos_hip);
	
	const long double dalpha_dt_qt = num / den ;
	
	return dalpha_dt_qt;	
}

/*************************************************************************************
 * 
 * Background - classical
 * 
 *************************************************************************************/


static double
Expo_dmu_m (const double alpha, const double mu_m)
{
	const double sqrt_2_2 = 0.5 * sqrt (2.0);
	const double mu_m2    = mu_m * mu_m; 
	const double dmu_m    = - 3.0 * (2.0 * mu_m - mu_m2) * (- (1.0 - mu_m) + sqrt_2_2);

	return dmu_m;
}

static double
Expo_dmu_p (const double alpha, const double mu_p)
{
	const double sqrt_2_2 = 0.5 * sqrt (2.0);
	const double mu_p2    = mu_p * mu_p;
	const double dmu_p    = + 3.0 * (2.0 * mu_p - mu_p2) * (-(mu_p - 1.0) + sqrt_2_2);

	return dmu_p;
}

static double 
Expo_w (const double alpha, const double x)
{
	const double x2 =  x * x; 
	const double w  = 2.0 * x2 - 1.0;

	return w;	
}

//~ static double 
//~ Expo_a (const double alpha, ExpoConsts*c1)
//~ {
	//~ const double a = c1->lpl * exp (alpha);
//~ 
	//~ return a;
//~ }

static double 
Expo_V0 (const double H, const double phi, const double mu)
{
	const double exp_phi	= exp ( 3.0 * sqrt(2) *  phi );
	const double onemx2		= mu * (2.0 - mu);
	const double V0			= H * H * exp_phi * sqrt (3.0) * onemx2;
	
	return V0;
}

static double 
Expo_mu_approx (const double Hi, const double Hf)
{
	const double num 	= sqrt(2.0) * log(2.0) + log(12.0 - 8.0 * sqrt(2.0))-2.0 * log (Hf / Hi);
	const double den	= sqrt(2.0)-2.0;
	const double mu		= exp(num / den);
	
	return mu;
}

/*************************************************************************************
 * 
 * Background - classical only - EoM
 * 
 *************************************************************************************/

static int 
Expo_f_cl_m_only (realtype alpha, N_Vector y_cl_CI, N_Vector ydot_cl_CI, void *f_data)
{
	const double mu_m		= NV_Ith_S (y_cl_CI, 0);
	const long double H		= NV_Ith_S (y_cl_CI, 1);				
	const double dmu_m		= Expo_dmu_m (alpha, mu_m);
	const double x2			= (1.0 - mu_m) * (1.0 - mu_m);
	
	NV_Ith_S (ydot_cl_CI, 0) = dmu_m;
	NV_Ith_S (ydot_cl_CI, 1) = - 3.0 * x2 * H;
  
	return 0;
}

/*************************************************************************************
 * 
 * Background + Perturbation - classical - contraction - EoM
 * 
 *************************************************************************************/
static int 
Expo_f_cl_m (realtype alpha, N_Vector y_cl, N_Vector ydot_cl, void *f_data)
{
	ExpoConsts *c1					= (ExpoConsts *) f_data;
	const double mu_m				= NV_Ith_S (y_cl, 0);
	const long double H				= NV_Ith_S (y_cl, 1);
	const double zeta_re			= NV_Ith_S (y_cl, 2);		
	const double dzeta_re			= NV_Ith_S (y_cl, 3);
	const double zeta_im			= NV_Ith_S (y_cl, 4);		
	const double dzeta_im			= NV_Ith_S (y_cl, 5);	
	const double dmu_m				= Expo_dmu_m (alpha, mu_m);
	const double dzeta_dalpha_re	= Expo_dzeta_dalpha (alpha, mu_m, H, zeta_re, dzeta_re, c1);
	const double dzeta_dalpha_im	= Expo_dzeta_dalpha (alpha, mu_m, H, zeta_im, dzeta_im, c1);
	const double alpha_H			= exp (alpha) * H;
	const double x2					= (1.0 - mu_m) * (1.0 - mu_m);
	
	NV_Ith_S (ydot_cl, 0) = dmu_m;
	NV_Ith_S (ydot_cl, 1) = - 3.0 * x2 * H;
	NV_Ith_S (ydot_cl, 2) = dzeta_re / alpha_H;
	NV_Ith_S (ydot_cl, 3) = dzeta_dalpha_re;
	NV_Ith_S (ydot_cl, 4) = dzeta_im / alpha_H;
	NV_Ith_S (ydot_cl, 5) = dzeta_dalpha_im;
	
	//~ printf("zeta %20.15e dzeta %20.15e dzeta_alpha2 %20.15e \n", zeta_re, dzeta_re, dzeta_re / alpha_H);
  
	return 0;
}

/*************************************************************************************
 * 
 * Background - classical - expansion - EoM
 * 
 *************************************************************************************/
//~ 
static int 
Expo_f_cl_p_only (realtype alpha, N_Vector y_cl, N_Vector ydot_cl, void *f_data)
{
	const double mu_p		= NV_Ith_S (y_cl, 0);
	const long double H		= NV_Ith_S (y_cl, 1);				
	const double dmu_p		= Expo_dmu_p (alpha, mu_p);
	const double x2			= (mu_p - 1.0)* (mu_p - 1.0);
	
	NV_Ith_S (ydot_cl, 0) = dmu_p;
	NV_Ith_S (ydot_cl, 1) = - 3.0 * x2 * H;
  
	return 0;
}

/*************************************************************************************
 * 
 * Background + Perturbation - classical - expansion
 * 
 *************************************************************************************/

static int 
Expo_f_cl_p (realtype alpha, N_Vector y_cl, N_Vector ydot_cl, void *f_data)
{
	//~ printf("f_cl_p %20.15g %20.15g\n", NV_Ith_S (y_cl, 2) , NV_Ith_S (y_cl, 3));
	ExpoConsts *c1					= (ExpoConsts *) f_data;
	const double mu_p				= NV_Ith_S (y_cl, 0);
	const long double H				= NV_Ith_S (y_cl, 1);
	const double zeta_re			= NV_Ith_S (y_cl, 2);		
	const double dzeta_re			= NV_Ith_S (y_cl, 3);
	const double zeta_im			= NV_Ith_S (y_cl, 4);		
	const double dzeta_im			= NV_Ith_S (y_cl, 5);	
	const double dmu_p				= Expo_dmu_p (alpha, mu_p);
	const double dzeta_dalpha_re	= Expo_dzeta_dalpha(alpha, mu_p, H, zeta_re, dzeta_re, c1);
	const double dzeta_dalpha_im	= Expo_dzeta_dalpha(alpha, mu_p, H, zeta_im, dzeta_im, c1);
	const double alpha_H			= exp(alpha) * H;
	const double x2					= (1.0 - mu_p) * (1.0 - mu_p);
	
	NV_Ith_S (ydot_cl, 0) = dmu_p;
	NV_Ith_S (ydot_cl, 1) = - 3.0 * x2 * H;
	NV_Ith_S (ydot_cl, 2) = dzeta_re / alpha_H ;
	NV_Ith_S (ydot_cl, 3) = dzeta_dalpha_re;
	NV_Ith_S (ydot_cl, 4) = dzeta_im / alpha_H ;
	NV_Ith_S (ydot_cl, 5) = dzeta_dalpha_im ;
  
	return 0;
}

/*************************************************************************************
 * 
 * Background - classical only - H scale root
 * 
 *************************************************************************************/

static int 
Expo_H_scale_de_root (realtype alpha, N_Vector y, realtype *gout, void *user_data)
{
	ExpoConsts *c1			= (ExpoConsts *) user_data;
	const double mu_p		= NV_Ith_S (y, 0);
	const double H			= NV_Ith_S (y, 1);
	const double absH		= fabs (H);
	const double epsilon	= Expo_epsilon(alpha, mu_p, H, c1);

	gout[0] = log (absH / c1->H_scale);
	
	if (H < 0.0)
	{	
		gout[1] = 10.0;
		gout[2]	= epsilon - 1.0e-7;
		//~ printf("Searching for WKB condition eps %20.15e \n", epsilon);
	}
	else
	{
		gout[1] = log (mu_p + 1.0e-7);
		gout[2] = 10.0;
		//~ gout[3] = log(alpha - 19.0) ;
	}

	return 0;
} 

/*************************************************************************************
 * 
 * Background - classical - phi integral
 * 
 *************************************************************************************/

static int
Expo_phi (realtype alpha, N_Vector y_cl, N_Vector yQdot_cl, void *f_data)
{
	NV_Ith_S (yQdot_cl, 0) = 1.0 - NV_Ith_S (y_cl, 0);

	return 0;
}

/*************************************************************************************
 * 
 * Background - quantum
 * 
 *************************************************************************************/

static int 
Expo_f_qt_only (realtype t, N_Vector y_qt, N_Vector ydot_qt, void *f_data)
{
	ExpoConsts *c1		= (ExpoConsts *) f_data;
	const double alpha	= NV_Ith_S (y_qt, 0);
	const double phi	= NV_Ith_S (y_qt, 1); 

	const long double dphi_dt_qt   = Expo_dphi_dt_qt (t, alpha, phi, c1) ;
	const long double dalpha_dt_qt = Expo_dalpha_dt_qt (t, alpha, phi, c1) ;
	
	//~ printf(" CC dalpha %20.15e dphi %20.15e\n", dalpha_dt_qt, dphi_dt_qt);
  
	NV_Ith_S (ydot_qt, 0) = dalpha_dt_qt;
	NV_Ith_S (ydot_qt, 1) = dphi_dt_qt;
	
	return 0;
}

/*************************************************************************************
 * 
 * Background + Perturbation - quantum
 * 
 *************************************************************************************/
	

static int 
Expo_f_qt (realtype t, N_Vector y_qt, N_Vector ydot_qt, void *f_data)
{
	//~ printf("f_qte %20.15e %20.15e %20.15e %20.15e\n", NV_Ith_S (y_qt, 0), NV_Ith_S (y_qt, 1), NV_Ith_S (y_qt, 2), NV_Ith_S (y_qt, 3));
	//~ printf("pre RHS: %20.15e %20.15e %20.15e %20.15e\n", NV_Ith_S (y_qt, 0), NV_Ith_S (y_qt, 1), NV_Ith_S (y_qt, 2), NV_Ith_S (y_qt, 3) );
	ExpoConsts *c1					= (ExpoConsts *) f_data;
	const double alpha				= NV_Ith_S (y_qt, 0);
	const double phi				= NV_Ith_S (y_qt, 1);
	const double zeta_re			= NV_Ith_S (y_qt, 2);		
	const double dzeta_re			= NV_Ith_S (y_qt, 3);
	const double zeta_im			= NV_Ith_S (y_qt, 4);		
	const double dzeta_im			= NV_Ith_S (y_qt, 5);
	const double dphi	 			= Expo_dphi_dt_qt(t, alpha, phi, c1);
	const double da					= Expo_dalpha_dt_qt (t, alpha, phi, c1);
	//~ const double H					= c1->Rc * da;
	
	
	//~ printf("f_qt2 %20.15e %20.15e\n", alpha, phi);

	const double dzeta_dt_re	= Expo_dzeta_dt(t, alpha, phi, zeta_re, dzeta_re, c1);
	const double dzeta_dt_im	= Expo_dzeta_dt(t, alpha, phi, zeta_im, dzeta_im, c1);
	//~ const double teste = Expo_zl_z_qt (t, NV_Ith_S (y_qt, 0), H, NV_Ith_S (y_qt, 1), c1);
	
	//~ printf("f_qts v %20.15e %20.15e %20.15e %20.15e\n", alpha, phi, zeta_re, dzeta_re);
	//~ printf(" RHS4 dzeta_dt %20.15e zl_l %20.15e\n", dzeta_dt, teste);
	//~ printf(" RHS5 dalpha %20.15e dphi %20.15e\n", da, dphi);
	
	NV_Ith_S (ydot_qt, 0) = da;
	NV_Ith_S (ydot_qt, 1) = dphi;
	NV_Ith_S (ydot_qt, 2) = dzeta_re /(c1->Rc * exp(alpha));
	NV_Ith_S (ydot_qt, 3) = dzeta_dt_re/(c1->Rc);
	NV_Ith_S (ydot_qt, 4) = dzeta_im /(c1->Rc * exp(alpha));
	NV_Ith_S (ydot_qt, 5) = dzeta_dt_im/(c1->Rc);
	
	//~ printf("f_qts dv %20.15e %20.15e %20.15e %20.15e\n", c1->Rc * NV_Ith_S (ydot_qt, 0), NV_Ith_S (ydot_qt, 1), NV_Ith_S (ydot_qt, 2), NV_Ith_S (ydot_qt, 3));

	
	return 0;
}
/*************************************************************************************
 * 
 * Background - quantum - bounce condition
 * 
 *************************************************************************************/

static int 
Expo_bounce_H_scale_root (realtype t, N_Vector y, realtype *gout, void *user_data)
{
	ExpoConsts *c1		= (ExpoConsts *) user_data;
	const double alpha	= NV_Ith_S (y, 0);
	const double phi	= NV_Ith_S (y, 1);
	const double H      = c1->Rc * Expo_dalpha_dt_qt(t, alpha, phi, c1);
	
	gout[0] = H;
	gout[1] = log (fabs (H / c1->H_scale));

	return 0;
} 

/*************************************************************************************
 * 
 * Background - classical - effective de w condition
 * 
 *************************************************************************************/

//~ static int 
//~ Exp_eff_de (realtype t, N_Vector y, realtype *gout, void *user_data)
//~ {
	//~ const double mu	= NV_Ith_S (y, 0);
	//~ 
	//~ gout[0] = mu - 1.0;
//~ 
	//~ return 0;
//~ } 

/*************************************************************************************
 * 
 * Background - classical only
 * 
 *************************************************************************************/

SyState
run_evol_cl_only (void *cvode,
					double alpha_i,
					double alpha_f,
					N_Vector y_cl_only, 
					N_Vector yQ_cl_only, 
					ExpoConsts *c1,
					int verbose,
					int wkb,
					int sav, 
					FILE *out[2])
{
	int flag1, flag2;
	double alpha_step,
			w,
			H,
			mu_m,
			x,
			phi,
			tau,
			epsilon,
			zeta_re,
			dzeta_re,
			zeta_im, 
			dzeta_im,
			zetai, 
			dzetai_re,
			dzetai_im,
			alpha_p_d;
	
	int cont1 = 0;
		

	SyState ret0;

	flag1 = CVodeReInit (cvode, alpha_i, y_cl_only);
	CVODE_CHECK (&flag1, "CVodeInit", 1, ret0);

	flag1 = CVodeQuadReInit (cvode, yQ_cl_only);
	CVODE_CHECK (&flag1, "CVodeQuadInit", 1, ret0);

	while (TRUE)
	{
		flag1 = CVode (cvode, alpha_f, y_cl_only, &alpha_step, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag1, "CVode", 1, ret0);
		
		flag2 = CVodeGetQuad (cvode, &alpha_step, yQ_cl_only);
		CVODE_CHECK (&flag2, "CVodeGetQuad", 1,ret0);

		phi			= NV_Ith_S (yQ_cl_only, 0);
		mu_m		= NV_Ith_S (y_cl_only, 0);
		H			= NV_Ith_S (y_cl_only, 1);
		x			= 1.0 - mu_m;
		w			= Expo_w (alpha_step, x);           
		//~ V0			= Expo_V0 (H, phi, mu_m);
		epsilon		= Expo_epsilon (alpha_step, mu_m, H, c1);
		zeta_re		= Expo_zeta (alpha_step, mu_m, H, c1);
		dzeta_re	= Expo_dzeta_re (alpha_step, mu_m, H, c1);
		zeta_im		= Expo_zeta (alpha_step, mu_m, H, c1);
		dzeta_im	= Expo_dzeta_im (alpha_step, mu_m, H, c1);
		tau			= (H / (sqrt(H * H))) * (alpha_step - c1->alpha_B);
		double pok			= c1-> kpert * c1-> kpert * c1-> kpert * (zeta_re*zeta_re + zeta_im * zeta_im);

		if (verbose)
		{
			printf ("# CI % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
			alpha_step,
			mu_m,
			H,
			phi,
			w,
			zeta_re,
			dzeta_re,
			zeta_im,
			dzeta_im
			);
		 }
		 
		 if(sav == 1)
		 {
			 if(out[0]!=NULL)
			 {
				 fprintf(out[0],"% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
							tau,
							x,
							phi,
							H,
							w,
							sqrt(zeta_re * zeta_re + zeta_im * zeta_im), 
							sqrt(dzeta_re * dzeta_re + dzeta_im * dzeta_im),
							sqrt(zeta_re * zeta_re + zeta_im * zeta_im), 
							sqrt(dzeta_re * dzeta_re + dzeta_im * dzeta_im),
							pok, 
							c1->kpert,
							exp(fabs(tau)),
							log(fabs(zeta_re)),
						log(fabs(zeta_im))
							);
				}
		 }
		 
		if (flag1 == CV_ROOT_RETURN)
		{
			int roots[2];
			flag1 = CVodeGetRootInfo (cvode, roots);
			CVODE_CHECK (&flag1, "CVodeGetRootInfo", 1, ret0);
			
			if (verbose)
				printf ("# Some root found!\n");
			if (roots[0])
			{
				if (verbose)
					printf ("# Root found -- H scale!\n");
				if (H < 0.0)
					break;
			}
			if (roots[1])
			{
				if (H > 0.0)
				{
					if (verbose)
						printf ("# Root found -- effective de!\n");
					break;
				}
			}
			
			if (roots[2])
			{
				if(wkb==1)
				{
					if (verbose)
						printf ("# WKB limit found! %20.15g\n", epsilon);
					zetai		= zeta_re;
					dzetai_re	= dzeta_re;
					dzetai_im	= dzeta_im;
					break;
				}
			}
			
			
		
		}
		
			double delta_w = fabs(w -1.0);
			
			
			
			if((H<0)&&(delta_w < 0.5)&&(cont1 == 0))
			{
				alpha_p_d = alpha_step;
				cont1 = cont1 +1;
				
			}
		
	  }
	  	 
	ret0.alpha		= alpha_step;
	ret0.mu			= mu_m;
	ret0.t			= 0.0;
	//~ ret0.alpha_B	= 0.0;
	ret0.H			= H;
	ret0.phi		= phi;
	ret0.zeta_re	= zetai;
	ret0.dzeta_re	= dzetai_re;
	ret0.zeta_im	= zetai;
	ret0.dzeta_im	= dzetai_im;
	ret0.lab1		= alpha_p_d;
	
	return ret0; 
}

/*************************************************************************************
 * 
 * Run_evol - background + perturbation - classical
 * 
 *************************************************************************************/

SyState
run_evol_cl(void *cvode,
			double alpha_i,
			double alpha_f,
			N_Vector y_cl, 
			N_Vector yQ_cl, 
			ExpoConsts *c1,
			int verbose,
			int sav,
			FILE *out[2])
{
	int flag1, flag2;
	double alpha_step,
			w,
			H,
			mu_m,
			x,
			phi,
			zeta_re,
			dzeta_re,
			zeta_im,
			dzeta_im,
			tau,
			//~ alpha_e,
			//~ zeta_e,
			//~ dzeta_e,
			//~ alpha_s,
			//~ zeta_s,
			//~ dzeta_s,
			//~ epsilon, 
			zeta_WKB,
			dzeta_WKB;
				
		
	SyState ret0;		
	flag1 = CVodeReInit (cvode, alpha_i, y_cl);
	CVODE_CHECK (&flag1, "CVodeInit", 1, ret0);

	flag1 = CVodeQuadReInit (cvode, yQ_cl);
	CVODE_CHECK (&flag1, "CVodeQuadInit", 1, ret0);
		
	while (TRUE)
	{
		
		flag1 = CVode (cvode, alpha_f, y_cl, &alpha_step, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag1, "CVode", 1, ret0);
		
		flag2 = CVodeGetQuad (cvode, &alpha_step, yQ_cl);
		CVODE_CHECK (&flag2, "CVodeGetQuad", 1,ret0);

		phi			= NV_Ith_S (yQ_cl, 0);
		mu_m		= NV_Ith_S (y_cl, 0);
		H			= NV_Ith_S (y_cl, 1);
		zeta_re		= NV_Ith_S (y_cl, 2);
		dzeta_re	= NV_Ith_S (y_cl, 3);
		zeta_im		= NV_Ith_S (y_cl, 4);
		dzeta_im	= NV_Ith_S (y_cl, 5);
		
		//~ printf("AA %20.15g %20.15g\n", NV_Ith_S (y_cl, 2), NV_Ith_S (y_cl, 3));
		
		if(H>0)
		{
			x			=  mu_m -1.0;
		}
		else
		{
			x			= 1.0 - mu_m;
		}
		
		w			= Expo_w (alpha_step, x);           
		//~ V0			= Expo_V0 (H, phi, mu_m);
		//~ epsilon		= Expo_epsilon(alpha_step, mu_m, H, c1);
		zeta_WKB	= sqrt(2.0 * Expo_zeta(alpha_step, mu_m, H, c1) * Expo_zeta(alpha_step, mu_m, H, c1));
		dzeta_WKB	= sqrt(Expo_dzeta_re(alpha_step, mu_m, H, c1) * Expo_dzeta_re(alpha_step, mu_m, H, c1) 
		+ Expo_dzeta_im(alpha_step, mu_m, H, c1) * Expo_dzeta_im(alpha_step, mu_m, H, c1));
		tau			= (H / (sqrt(H * H))) * (alpha_step - c1->alpha_B);
		double pok = (c1->kpert) * (c1->kpert) * (c1->kpert) * (zeta_re * zeta_re + zeta_im * zeta_im);
	
		if (verbose)
		{
			printf ("# CP % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
			alpha_step,
			mu_m,
			H,
			phi,
			w,
			zeta_re, 
			dzeta_re,
			zeta_im,
			dzeta_im
			);
		 }
		 
		 if(sav == 1)
		 {
			 if(out[0]!=NULL)
			 {
				 fprintf(out[0],"% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g  % 20.15g % 20.15g % 20.15g  % 20.15g % 20.15g % 20.15g\n",
							tau,
							x,
							phi,
							H,
							w,
							sqrt(zeta_re * zeta_re + zeta_im * zeta_im), 
							sqrt(dzeta_re * dzeta_re + dzeta_im * dzeta_im),
							zeta_WKB,
							dzeta_WKB,
							pok, 
							c1->kpert,
							exp(fabs(tau)),
							log(fabs(zeta_re)),
						log(fabs(zeta_im))
							);
			 }
		}
		 
		if (flag1 == CV_ROOT_RETURN)
		{
			int roots[2];
			flag1 = CVodeGetRootInfo (cvode, roots);
			CVODE_CHECK (&flag1, "CVodeGetRootInfo", 1, ret0);
			
			printf ("# Some root found!\n");
			if (roots[0])
			{
				if (verbose)
					printf ("# Root found -- H scale!\n");
				if (H < 0.0)
					break;
			}
			//~ if (roots[1])
			//~ {
				//~ if (H > 0.0)
				//~ {
					//~ if (verbose)
					//~ {
						//~ printf ("# Root found -- effective de!\n");
					//~ }
					//~ 
					//~ break;
				//~ }
			//~ }
						
			//~ if (roots[3])
			//~ {
				//~ alpha_ns 	= alpha_step;
				//~ zeta_ns  	= sqrt(zeta_re * zeta_re + zeta_im * zeta_im);
				//~ dzeta_ns	= sqrt(dzeta_re * dzeta_re + dzeta_im * zeta_im);
				//~ printf("#Crossing alpha %20.15g", alpha_ns);
			//~ }
								
		}
		
			if((H>0)&&(w < 1.0/3.0))
			{
				printf("new condition, w = %20.15g", w);
				break;
			}
			
			
		
	  }
	  	 
	ret0.alpha		= alpha_step;
	ret0.mu			= mu_m;
	ret0.t			= 0.0;
	//~ ret0.alphab		= 0.0;
	ret0.H			= H;
	ret0.phi		= phi;
	ret0.zeta_re	= zeta_re;
	ret0.dzeta_re	= dzeta_re;

	//~ ret0.zeta_e		= zeta_e;
	//~ ret0.dzeta_e	= dzeta_e;
	//~ ret0.alpha_e	= alpha_e;
	//~ ret0.zeta_s	= zeta_s;
	//~ ret0.dzeta_s	= dzeta_s;
	//~ ret0.alpha_s	= alpha_s;
	
	
	return ret0; 
}

/*************************************************************************************
 * 
 * Run_evol - background - quantum
 * 
 *************************************************************************************/

SyState
run_evol_qt_only (void *cvode,
				double t_i,
				N_Vector y_qt,
				double max_step,
				ExpoConsts *c1,
				int verbose, int cont, double alpha_B, int cont1, int cont2)
{

	double last_t	= t_i;
	double t_f      = 1.0e20;
	double alpha_i = NV_Ith_S (y_qt, 0);
	double phi_i   = NV_Ith_S (y_qt, 1);
	double H_i     = c1->Rc * Expo_dalpha_dt_qt (t_i, alpha_i, phi_i, c1);
	SyState ret;
	double alpha, phi, H, mu_p, x, V0, alpha_H_min, alpha_H_max, dH;
	double t, w = 1.0;
	int flag;
	
	
	flag = CVodeReInit (cvode, t_i, y_qt);
	CVODE_CHECK (&flag, "CVodeInit", 1, ret);

	if (H_i < 0.0)
	{
		flag = CVodeSetMaxStep (cvode, max_step);
		CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, ret);
	}
	else
	{
		flag = CVodeSetMaxStep (cvode, 0.0);
		CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, ret);
	} 
	
	while (TRUE)
	{
		flag = CVode (cvode, t_f, y_qt, &t, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag, "CVode", 1, ret);

		alpha	= NV_Ith_S (y_qt, 0);
		phi		= NV_Ith_S (y_qt, 1);
		H 		= c1->Rc * Expo_dalpha_dt_qt (t, alpha, phi, c1);
		x		= Expo_xqt (alpha, phi, c1);
		mu_p	= 1.0 + x;
		V0		= Expo_V0 (H, phi, mu_p);
		dH		= - 0.5 * Expo_dphi_dt_qt(t, alpha, phi, c1);


		
		//~ printf(">>H %20.15g dH %20.15g alpha_h_max %20.15g\n", H, dH, alpha);
		
		if (verbose)
		{
			 printf ("# QT % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
			 alpha,
			 mu_p,
			 H,
			 phi,
			 w, 
			 V0
			 );
		}
		  		 
		if (flag == CV_ROOT_RETURN)
		{
			int roots[2];
			flag = CVodeGetRootInfo (cvode, roots);
			CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, ret);
			if (roots[0])
			{
				if (verbose)
					{
						printf ("# Bounce found! %20.15g\n", alpha);
					}
				if(cont == 1)
				{
				alpha_B	= alpha;
				cont = cont + 1;					
				}
			}
			if (roots[1])
			{
				if (phi < 0.0 && fabs (phi) > 0.1)
				{
					if (verbose)
						printf ("# Scale found! %20.15g\n", alpha);
					break;
				}
				else
				{
					if (verbose)
						printf ("# Spurious scale found! %20.15g\n", alpha);
				}
			}			
		}
		
		if (fabs ((last_t - t) / t) < 1.0e-3)
		{
			//~ printf("tempo saturou : alpha % 20.15Le phi % 20.15Le t_step % 20.15g \n", alpha, phi, t_step);
			ret.t		= 0.0;
			ret.alpha	= alpha;
			ret.phi		= phi;
			ret.H		= H;
			ret.alpha_B	= alpha_B;
			if(cont==2)
			{
				printf ("##return Bounce found! %20.15g\n", alpha_B);
			return  run_evol_qt_only (cvode, 0.0, y_qt, t * 1.0e-3, c1, verbose, 2, alpha_B, cont1, cont2);
			
			}
			else
			{
			return run_evol_qt_only (cvode, 0.0, y_qt, t * 1.0e-3, c1, verbose, 1, 1, cont1, cont2);
			}
		}
		else
		{
			last_t = t;
		}
		
		if((H<0)&&(dH>0)&&(cont1==0))
		{
			alpha_H_min	= alpha;
			cont1 = cont1 + 1;
			//~ printf("##H %20.15g dH %20.15g alpha_h_min %20.15g\n", H, dH, alpha);
		}

		if((H>0)&&(dH<0)&&(cont2==0))
		{
			alpha_H_max	= alpha;
			cont2 = cont2 + 1;		
			//~ printf("##H %20.15g dH %20.15g alpha_h_max %20.15g\n", H, dH, alpha);
		}

	}

	ret.alpha	= NV_Ith_S(y_qt, 0);
	ret.phi		= NV_Ith_S(y_qt, 1);
	ret.t		= 0.0;
	ret.H		= H;
	ret.mu		= mu_p;
	ret.alpha_B	= alpha_B;
	ret.lab1	= alpha_H_min;
	ret.lab2	= alpha_H_max;
	return ret;
}

/*************************************************************************************
 * 
 * Run_evol - Background + Perturbations - Quantum
 * 
 *************************************************************************************/

SyState
run_evol_qt (void *cvode,
			double t_i,
			N_Vector y_qt,
			double max_step,
			ExpoConsts *c1,
			int verbose,
			int sav,
			FILE *out[2])
{
	double last_t	= t_i;
	double t_f      = 1.0e20;
	double alpha_i 	= NV_Ith_S (y_qt, 0);
	double phi_i   	= NV_Ith_S (y_qt, 1);
	
	
	//~ printf("e_qt %20.15e % 20.15e\n", alpha_i, phi_i);
	
	//~ printf("y_qt %20.15e %20.15e %20.15e %20.15e", NV_Ith_S (y_qt, 0), NV_Ith_S (y_qt, 1), NV_Ith_S (y_qt, 2), NV_Ith_S (y_qt, 3) );
	//~ double zeta_i  	= NV_Ith_S (y_qt, 2);
	//~ double dzeta_i	= NV_Ith_S (y_qt, 3);
	double H_i     = c1->Rc * Expo_dalpha_dt_qt (t_i, alpha_i, phi_i, c1);
	SyState ret;
	double alpha, 
			phi,
			H,
			mu_p,
			x, 
			zeta_re,
			dzeta_re,
			zeta_im,
			dzeta_im,
			tau;
			
	double t, w = 1.0;
	
	int flag;

	flag = CVodeReInit (cvode, t_i, y_qt);
	CVODE_CHECK (&flag, "CVodeInit", 1, ret);

//~ printf("y_qt %20.15e %20.15e %20.15e %20.15e", NV_Ith_S (y_qt, 0), NV_Ith_S (y_qt, 1), NV_Ith_S (y_qt, 2), NV_Ith_S (y_qt, 3) );


	if (H_i < 0.0)
	{
		flag = CVodeSetMaxStep (cvode, max_step);
		CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, ret);
	}
	else
	{
		flag = CVodeSetMaxStep (cvode, 0.0);
		CVODE_CHECK (&flag, "CVodeSetMaxStep", 1, ret);
	} 
	
	while (TRUE)
	{ 
		//~ printf("AA1\n");
		fflush(stdout);
		flag = CVode (cvode, t_f, y_qt, &t, CV_ONE_STEP /* CV_NORMAL */);
		CVODE_CHECK (&flag, "CVode", 1, ret);
		
		
		//~ printf("AA2!");
		alpha	 = NV_Ith_S (y_qt, 0);
		phi		 = NV_Ith_S (y_qt, 1);
		zeta_re	 = NV_Ith_S (y_qt, 2);
		dzeta_re = NV_Ith_S (y_qt, 3);
		zeta_im	 = NV_Ith_S (y_qt, 4);
		dzeta_im = NV_Ith_S (y_qt, 5);
		H 		 = c1->Rc * Expo_dalpha_dt_qt (t, alpha, phi, c1);
 //~ printf("AA3!");		
		x		= Expo_xqt (alpha, phi, c1);
		mu_p	= 1.0 + x;
		tau			= (H / (sqrt(H * H))) * (alpha - c1->alpha_B);
		//~ V0		= Expo_V0 (H, phi, mu_p);
		double zeta_WKB		= 0.0;
		double dzeta_WKB	= 0.0;
		double pok			= c1-> kpert * c1-> kpert * c1-> kpert * (zeta_re*zeta_re + zeta_im * zeta_im);
		
		if (verbose)
		{
			 printf ("# QT % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e % 20.15e\n",
			 alpha,
			 mu_p,
			 H,
			 phi,
			 w, 
			 zeta_re, 
			 dzeta_re,
			 zeta_im, 
			 dzeta_im
			 );
		}
		
		if(sav == 1)
		 {
			 if(out[0]!=NULL)
				{
					fprintf(out[0],"% 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g\n",
						tau,
						x,
						phi,
						H,
						w,
						sqrt(zeta_re*zeta_re + zeta_im * zeta_im), 
						sqrt(dzeta_re*dzeta_re + dzeta_im * dzeta_im),
						zeta_WKB,
						dzeta_WKB,
						pok,
						c1->kpert, 
						exp(fabs(tau)),
						log(fabs(zeta_re)),
						log(fabs(zeta_im))
						);
				}
		}
		 		 
		if (flag == CV_ROOT_RETURN)
		{
			int roots[2];
			flag = CVodeGetRootInfo (cvode, roots);
			CVODE_CHECK (&flag, "CVodeGetRootInfo", 1, ret);
			if (roots[0])
			{
				printf ("# Bounce found! %20.15g\n", alpha);
			}
			if (roots[1])
			{
				if (phi < 0.0 && fabs (phi) > 0.1)
				{
					printf ("# Scale found! %20.15g\n", alpha);
					break;
				}
				else
				{
					printf ("# Spurious scale found! %20.15g\n", alpha);
				}
			}			
		}
		
		if (fabs ((last_t - t) / t) < 1.0e-3)
		{
			//~ printf("tempo saturou : alpha % 20.15Le phi % 20.15Le t_step % 20.15g \n", alpha, phi, t);
			ret.t			= 0.0;
			ret.alpha		= alpha;
			ret.phi			= phi;
			ret.H			= H;
			ret.zeta_re		= zeta_re;
			ret.dzeta_re	= dzeta_re;
			ret.zeta_im		= zeta_im;
			ret.dzeta_im	= dzeta_im;
			
			return  run_evol_qt (cvode, 0.0, y_qt, t * 1.0e-3, c1, verbose, sav, out);
		}
		else
		{
			last_t = t;
		}

	}
	
	//~ printf(">>> %20.15g %20.15g %20.15g\n", alpha_ns, zeta_ns, dzeta_ns);

	ret.alpha	 	= NV_Ith_S(y_qt, 0);
	ret.phi		 	= NV_Ith_S(y_qt, 1);
	ret.zeta_re	 	= NV_Ith_S(y_qt, 2);
	ret.dzeta_re 	= NV_Ith_S(y_qt, 3);
	ret.zeta_im	 	= NV_Ith_S(y_qt, 4);
	ret.dzeta_im	= NV_Ith_S(y_qt, 5);
	ret.t			= 0.0;
	ret.H			= H;
	ret.mu			= mu_p;
	ret.zeta_ns		= log(NV_Ith_S(y_qt, 2) * NV_Ith_S(y_qt, 2) + NV_Ith_S(y_qt, 4) * NV_Ith_S(y_qt, 4));
	ret.pok_ns		= log((c1->kpert * c1->kpert * c1->kpert) * (NV_Ith_S(y_qt, 2) * NV_Ith_S(y_qt, 2) + NV_Ith_S(y_qt, 4) * NV_Ith_S(y_qt, 4)));

	return ret;
}


int 
main (int argc, char *argv[])
{
	void *cvode_cl_only_m	= CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_qt_only		= CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_cl_p_only	= CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_cl_m		= CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_cl_p	 	= CVodeCreate (CV_BDF, CV_NEWTON);
	void *cvode_qt	 		= CVodeCreate (CV_BDF, CV_NEWTON);
	/*void *cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);*/
	unsigned int dim_cl_only	= 2;
	unsigned int dim_qt_only	= 2; 
	unsigned int dim_cl			= 6;
	unsigned int dim_qt			= 6; 
	double alpha_i = 1000.0;
	double t_i = 0.0;
	int verb;
	int i_max;
	int sav = 0;
	//~ double phi_i, mu_i;
	//~ long double H_i;
	double abstol  = 0.0;
	double reltol  = 1.0e-15;
	int flag; 
	SyState  cl_WKB, cl_only, qt_only, cl_p_only, qt_state, cl_state; 
//~ 
	char format1[100];
	char format2[100];
	char format3[100];
	char format4[100];
	FILE *out[2] = {NULL, };
	int nome;
	
	ExpoConsts c1;
	c1.Rc    	= 4.69e62;//1.0e60;2.2e-18;
	c1.lpl    	= 1.616097e-33; /*cm*/
	c1.aii    	= 0.0;//5.3907205e-44  /*s*/;
	c1.kappa  	= sqrt(8.0 * M_PI) /*G=1*/;
	c1.sigma  	= 0.47e0;
	c1.d      	= -7.e-6;
	c1.H_scale  = 1.0e40;
	c1.alpha0   = 0.0;
	//~ c1.kpert	= 1.0e-3;
//~ 
  switch (argc)
  {
	  case 2:
      verb = 1;
      i_max = atof (argv[1]);
      break;
      case 4:
      verb = 0;
      i_max = atof (argv[1]);
      sav	= atof (argv[2]);
      nome	= atof (argv[3]);
      //~ sprintf (format1, "%s-VAR.dat", argv[1]);
      //~ sprintf (format2, "%s-PAR.dat", argv[1]);
      //~ out[0] = fopen (format1, "w");
      //~ out[1] = fopen (format2, "w");
      break;
      default:
      fprintf (stderr, "usage:\n 'main i_max'\n or\n 'main i_max 0' to spectral index\n or\n 'main i_max 1 file_number' to save data\n or\n 'main i_max 2 file_number' to save spectral index data\n");
      exit (-1);
      break;
  }
  
	N_Vector y0_cl_only		= N_VNew_Serial (dim_cl_only);
	N_Vector yQ0_cl_only	= N_VNew_Serial (dim_cl_only);
	N_Vector y0_qt_only		= N_VNew_Serial (dim_qt_only);
	N_Vector y0_cl			= N_VNew_Serial (dim_cl);
	N_Vector yQ0_cl			= N_VNew_Serial (dim_cl);
	N_Vector y0_qt			= N_VNew_Serial (dim_qt);

/*************************************************************************************
 * 
 * Initializing classical background only evolution
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_cl_only_m, &Expo_f_cl_m_only, alpha_i, y0_cl_only);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeQuadInit (cvode_cl_only_m, &Expo_phi, yQ0_cl_only);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, -1);

	flag = CVodeSStolerances (cvode_cl_only_m, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_cl_only_m, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_cl_only_m, dim_cl_only);
	CVODE_CHECK (&flag, "CVDense", 1, -1);

	flag = CVodeRootInit (cvode_cl_only_m, 3, &Expo_H_scale_de_root); 
	CVODE_CHECK (&flag, "CVodeRootInit", 1, -1);
	
	
/*************************************************************************************
 * 
 * Initializing quantum evolution - only
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_qt_only, &Expo_f_qt_only, t_i, y0_qt_only);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeSStolerances (cvode_qt_only, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_qt_only, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_qt_only, dim_qt_only);
	CVODE_CHECK (&flag, "CVDense", 1, -1);

	flag = CVodeRootInit (cvode_qt_only, 2, &Expo_bounce_H_scale_root); 
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, -1);
		
/*************************************************************************************
 * 
 * Initializing classical background only evolution - expanding
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_cl_p_only, &Expo_f_cl_p_only, alpha_i, y0_cl_only);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeQuadInit (cvode_cl_p_only, &Expo_phi, yQ0_cl_only);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, -1);

	flag = CVodeSStolerances (cvode_cl_p_only, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_cl_p_only, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_cl_p_only, dim_cl_only);
	CVODE_CHECK (&flag, "CVDense", 1, -1);

	flag = CVodeRootInit (cvode_cl_p_only, 3, &Expo_H_scale_de_root); 
	CVODE_CHECK (&flag, "CVodeRootInit", 1, -1);

/*************************************************************************************
 * 
 * Initializing classical evolution
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_cl_m, &Expo_f_cl_m, alpha_i, y0_cl);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeQuadInit (cvode_cl_m, &Expo_phi, yQ0_cl);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, -1);

	flag = CVodeSStolerances (cvode_cl_m, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_cl_m, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_cl_m, dim_cl);
	CVODE_CHECK (&flag, "CVDense", 1, -1);
	
	flag = CVodeRootInit (cvode_cl_m, 3, &Expo_H_scale_de_root); 
	CVODE_CHECK (&flag, "CVodeRootInit", 1, -1);

/*************************************************************************************
 * 
 * Initializing quantum evolution
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_qt, &Expo_f_qt, t_i, y0_qt);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeSStolerances (cvode_qt, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_qt, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_qt, dim_qt);
	CVODE_CHECK (&flag, "CVDense", 1, -1);

	flag = CVodeRootInit (cvode_qt, 2, &Expo_bounce_H_scale_root); 
	CVODE_CHECK (&flag, "CVodeSetStopTime", 1, -1);
/*************************************************************************************
 * 
 * Initializing classical evolution - expanding
 * 
 *************************************************************************************/
	flag = CVodeInit (cvode_cl_p, &Expo_f_cl_p, alpha_i, y0_cl);
	CVODE_CHECK (&flag, "CVodeInit", 1, -1);

	flag = CVodeQuadInit (cvode_cl_p, &Expo_phi, yQ0_cl);
	CVODE_CHECK (&flag, "CVodeQuadInit", 1, -1);

	flag = CVodeSStolerances (cvode_cl_p, reltol, abstol);
	CVODE_CHECK (&flag, "CVodeSStolerances", 1, -1);

	flag = CVodeSetUserData (cvode_cl_p, &c1);
	CVODE_CHECK (&flag, "CVodeSetUserData", 1, -1);

	flag = CVDense (cvode_cl_p, dim_cl);
	CVODE_CHECK (&flag, "CVDense", 1, -1);

	flag = CVodeRootInit (cvode_cl_p, 3, &Expo_H_scale_de_root); 
	CVODE_CHECK (&flag, "CVodeRootInit", 1, -1);

/*************************************************************************************
 * 
 * Grid I
 * 
 *************************************************************************************/	
int i;

if(sav == 2)
	{
		sprintf (format3, "%d-ns-VAR.dat", nome);
		sprintf (format4, "%d-ns-PAR.dat", nome);
		out[0] = fopen (format3, "w");
		out[1] = fopen (format4, "w");
		
		
	}
	
	const double cp_dist  	= 1.0e-12;
	const double sqrt_2_2	= 0.5 * sqrt(2.0);
	double phi_i    		= 0.0;
	double mu_m_i			= 1.0 - (sqrt_2_2 + cp_dist);
	double H_i				= -1.0e-16;
	double del_phi_i		= 25.0;

/******lab parameters********/
long double cond1, cond2, cond3, cond4;
double tau_p_d, tau_m_1, H_1, delta_H_1, delta_x_1, 
	 tau_m_2, H_2,  tau_0;
double Calpha_i, Cphi_i, Cmu_m_i, CH_i;

/******lab parameters********/

for (i=0; i<= i_max; i++)
{
	double lnk_i;
	
	if(i_max!=0)
	{
		lnk_i		= -3.0 * M_LN10 + 6.0 * M_LN10 * i/(i_max);
	}
	else
	{
		lnk_i		= -3.0 * M_LN10;
	}
	
	
   	//~ double alpha_p, phi_p, zeta_p, dzeta_p;
	c1.kpert	= exp(lnk_i);
	printf("# Entrei no for, k = % 20.15g verb = % 20.15d i_max = % 20.15d \n",  c1.kpert, verb, i_max);
	
	if(sav == 1)
	{
		sprintf (format1, "%d-VAR.dat", nome);
		sprintf (format2, "%d-PAR.dat", nome);
		out[0] = fopen (format1, "w");
		out[1] = fopen (format2, "w");
	}
	
	
/*************************************************************************************
 * 
 * Calibrating contraction phase and perturbations initial conditions
 * 
 *************************************************************************************/
	{		
		//~ double zeta_i_re		= 0.0;
		//~ double dzeta_i_re		= 0.0;
		
		alpha_i = 1000.0;
    
		NV_Ith_S (y0_cl_only, 0)	= mu_m_i;
		NV_Ith_S (y0_cl_only, 1)	= H_i;
		NV_Ith_S (yQ0_cl_only, 0)	= phi_i;

		cl_only = run_evol_cl_only (cvode_cl_only_m, alpha_i, 0.0, y0_cl_only, yQ0_cl_only, &c1, 0, 0, 0, out);
		{
			double H			= cl_only.H;
			double x            = 1.0 - cl_only.mu;
			double alpha_q_ini	= (1.0 / 3.0) * log (( c1.d * c1.Rc) / (x * H));
			double phi_q_ini   	= del_phi_i / (c1.sigma * c1.sigma * alpha_q_ini);
			
			double delta_alpha  = cl_only.alpha - alpha_i;
			double delta_phi    = cl_only.phi   - phi_i;
			
			alpha_i		= alpha_q_ini - delta_alpha;
			phi_i  		= phi_q_ini - delta_phi;
			
			Calpha_i = alpha_i;
			Cmu_m_i  = mu_m_i;
			CH_i     = H_i;
			Cphi_i   = phi_i;
			
			NV_Ith_S (y0_cl_only, 0)	= mu_m_i;
			NV_Ith_S (y0_cl_only, 1)	= H_i;
			NV_Ith_S (yQ0_cl_only, 0)	= phi_i;
		
			c1.V0 = Expo_V0 (cl_only.H, phi_q_ini, cl_only.mu);
			
			
			{
				const long double arg_tri   = 2.0 * alpha_q_ini * c1.d;
				const long double arg_hip   = c1.sigma * c1.sigma * alpha_q_ini * phi_q_ini;
				const long double cos_tri   = cosl (arg_tri);
				const long double cos_hip   = coshl (arg_hip);
				const long double sen_tri   = sinl (arg_tri);
				const long double sen_hip   = sinhl (arg_hip);
				cond1     = (2.0 * c1.d * sen_hip) / (phi_q_ini * c1.sigma * c1.sigma * sen_tri);
				cond2     = (2.0 * c1.d * cos_hip) / (- alpha_q_ini * c1.sigma * c1.sigma * sen_tri + 2.0 * c1.d * cos_tri);
				cond3     = cos_hip / cos_tri ;
				cond4		= tanhl (arg_hip);
				printf ("# cond1 % 20.15Le cond2 % 20.15Le cond3 % 20.15Le cond4 % 20.15Le \n# alpha_match = % 20.15e phi_match = % 20.15e \n", 
					cond1,
					cond2,
					cond3,
					cond4,
					alpha_q_ini,
					phi_q_ini);
			}
		}
/*************************************************************************************
 * 
 * Contracting phase + Exansion phase -- ONLY
 * 
 *************************************************************************************/
		{
			printf("d %20.15g\n sigma %20.15g\n H_i %20.15g\n delta_ phi %20.15g\n", c1.d, c1.sigma,H_i,
					del_phi_i);
			
			printf("cond1 %20.15Lg\n cond2 %20.15Lg\n cond3 %20.15Lg\n cond4 %20.15Lg\n alpha_i %20.15g\n phi_i %20.15g\n", cond1,cond2, cond3,
					cond4, Calpha_i, Cphi_i);
			
			cl_WKB = run_evol_cl_only (cvode_cl_only_m, alpha_i, 0.0, y0_cl_only, yQ0_cl_only, &c1, 0, 1, 0, out);
			
			alpha_i						= cl_WKB.alpha;
			NV_Ith_S (y0_qt_only, 0)	= cl_WKB.alpha;
			NV_Ith_S (y0_qt_only, 1)	= cl_WKB.phi;
			
			H_1 = cl_WKB.H;
			
			
			delta_H_1 = fabs(c1.Rc *Expo_dalpha_dt_qt(t_i, cl_WKB.alpha, cl_WKB.phi, &c1) - cl_WKB.H) / fabs(cl_WKB.H);
			delta_x_1 = fabs(Expo_xqt(cl_WKB.alpha, cl_WKB.phi, &c1) + (cl_WKB.mu -1)) / fabs(cl_WKB.mu -1);
									
			qt_only	= run_evol_qt_only (cvode_qt_only, t_i, y0_qt_only, 0.0, &c1, 0  /*verbose*/, 1 /* cont*/, 0 /* alpha_B*/, 0 /* cont1*/, 0  /* cont2*/);

			NV_Ith_S (y0_cl_only, 0)	= Expo_mu_approx (qt_only.H, 1.0);
			NV_Ith_S (y0_cl_only, 1)	= qt_only.H;
			NV_Ith_S (yQ0_cl_only, 0)	= qt_only.phi;
			alpha_i 	= qt_only.alpha;
			c1.alpha_B	= qt_only.alpha_B;
			
			printf("qt_only_alpha_B  %20.15g\n", qt_only.alpha_B);
			
			tau_p_d 	= - (cl_WKB.lab1 - qt_only.alpha_B);
			printf("tau_p_d  %20.15g\n",tau_p_d );
			tau_m_1 	= - (cl_WKB.alpha - qt_only.alpha_B);
			printf("tau_m_1  %20.15g\n",tau_m_1 );
			printf("H_1 %20.15g\n delta_H_1 %20.15g\n delta_x_1 %20.15g\n", H_1, delta_H_1, delta_x_1);
			//~ tau_H_min	= - (qt_only.lab1 -  qt_only.alpha_B);
			//~ printf("tau_H_min  %20.15g\n",tau_H_min );
			printf("alpha_B  %20.15g\n",c1.alpha_B );
			//~ tau_H_max	= (qt_only.lab2 -  qt_only.alpha_B);
			//~ printf("tau_H_max  %20.15g\n",tau_H_max );
			tau_m_2		= (qt_only.alpha - qt_only.alpha_B);
			printf("tau_m_2 %20.15g\n",tau_m_2 );
			
			H_2			= qt_only.H;
			printf("H_2 %20.15g\n", H_2);
			
			cl_p_only = run_evol_cl_only (cvode_cl_p_only, alpha_i, 1.0e3, y0_cl_only, yQ0_cl_only, &c1, 0, 0, 0, out);	

			c1.alpha0	= cl_p_only.alpha;		
			tau_0		= cl_p_only.alpha - qt_only.alpha_B;
			printf("alpha0 %20.15g\n tau0 %20.15g\n", c1.alpha0, tau_0);
		
		}
	}
/*************************************************************************************
 * 
 * Contracting phase : Background + Perturbation
 * 
 *************************************************************************************/
	{
		alpha_i = Calpha_i;
		NV_Ith_S (y0_cl_only, 0)	= Cmu_m_i;
		NV_Ith_S (y0_cl_only, 1)	= CH_i;
		NV_Ith_S (yQ0_cl_only, 0)	= Cphi_i;
		
		//~ printf ("# PreWKB %20.15e %20.15e %20.15e %20.15e\n",
				//~ alpha_i,
				//~ NV_Ith_S (y0_cl_only, 0),
				//~ NV_Ith_S (y0_cl_only, 1),
				//~ NV_Ith_S (yQ0_cl_only, 0));
		
		cl_WKB = run_evol_cl_only (cvode_cl_only_m, alpha_i, 0.0, y0_cl_only, yQ0_cl_only, &c1, verb, 1, sav, out);
//~ 
		//~ printf ("# PosWKB %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ cl_WKB.alpha,
				//~ NV_Ith_S (y0_cl_only, 0),
				//~ NV_Ith_S (y0_cl_only, 1),
				//~ NV_Ith_S (yQ0_cl_only, 0),
				//~ cl_WKB.zeta_re,
				//~ cl_WKB.dzeta_re,
				//~ cl_WKB.zeta_im,
				//~ cl_WKB.dzeta_im);

		
		NV_Ith_S (y0_cl, 0) = cl_WKB.mu;
		NV_Ith_S (y0_cl, 1) = cl_WKB.H;
		NV_Ith_S (y0_cl, 2) = cl_WKB.zeta_re;
		NV_Ith_S (y0_cl, 3) = cl_WKB.dzeta_re;
		NV_Ith_S (y0_cl, 4) = cl_WKB.zeta_im;
		NV_Ith_S (y0_cl, 5) = cl_WKB.dzeta_im;
		NV_Ith_S (yQ0_cl, 0) = cl_WKB.phi;
		alpha_i	= cl_WKB.alpha;
		
		//~ printf ("# PreCLC %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ alpha_i,
				//~ NV_Ith_S (y0_cl, 0),
				//~ NV_Ith_S (y0_cl, 1),
				//~ NV_Ith_S (yQ0_cl, 0),
				//~ NV_Ith_S (y0_cl, 2), 
				//~ NV_Ith_S (y0_cl, 3),
				//~ NV_Ith_S (y0_cl, 4), 
				//~ NV_Ith_S (y0_cl, 5));				
				
		cl_state = run_evol_cl (cvode_cl_m, alpha_i, 0.0, y0_cl, yQ0_cl, &c1, verb, sav, out);
			
		NV_Ith_S (y0_qt, 0)	= cl_state.alpha;
		NV_Ith_S (y0_qt, 1)	= NV_Ith_S (yQ0_cl, 0);
		NV_Ith_S (y0_qt, 2)	= NV_Ith_S (y0_cl, 2);
		NV_Ith_S (y0_qt, 3)	= NV_Ith_S (y0_cl, 3);
		NV_Ith_S (y0_qt, 4)	= NV_Ith_S (y0_cl, 4);
		NV_Ith_S (y0_qt, 5)	= NV_Ith_S (y0_cl, 5);
			
		//~ printf ("# PosCLC %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ cl_state.alpha,
				//~ NV_Ith_S (y0_cl, 0),
				//~ NV_Ith_S (y0_cl, 1),
				//~ NV_Ith_S (yQ0_cl, 0),
				//~ NV_Ith_S (y0_cl, 2), 
				//~ NV_Ith_S (y0_cl, 3),
				//~ NV_Ith_S (y0_cl, 4), 
				//~ NV_Ith_S (y0_cl, 5));		
				
				//~ printf ("# PosCLC2 %20.15e %20.15e %20.15Le %20.15e %20.15e %20.15e\n",
				//~ cl_state.alpha,
				//~ cl_state.mu,
				//~ cl_state.H,
				//~ cl_state.phi,
				//~ cl_state.zeta_re, 
				//~ cl_state.dzeta_re);		

	
/*************************************************************************************
 * 
 * Quantum phase : Background + Perturbation
 * 
 *************************************************************************************/
		//~ printf ("# PreQTC %20.15e %20.15e %20.15Le %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ NV_Ith_S (y0_qt,0),
				//~ cl_state.mu,
				//~ cl_state.H,
				//~ NV_Ith_S (y0_qt,1),
				//~ NV_Ith_S (y0_qt,2), 
				//~ NV_Ith_S (y0_qt,3),
				//~ NV_Ith_S (y0_qt,4), 
				//~ NV_Ith_S (y0_qt,5));				

		qt_state = run_evol_qt (cvode_qt, 0.0, y0_qt, 0.0, &c1, verb, sav, out);
		
		//~ printf ("# PosQTC %20.15e %20.15e %20.15Le %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ NV_Ith_S (y0_qt,0),
				//~ Expo_mu_approx(qt_state.H, 1.0),
				//~ qt_state.H,
				//~ NV_Ith_S (y0_qt,1),
				//~ NV_Ith_S (y0_qt,2), 
				//~ NV_Ith_S (y0_qt,3),
				//~ NV_Ith_S (y0_qt,4), 
				//~ NV_Ith_S (y0_qt,5));
				
		//~ 
			//~ printf(">>> %20.15g %20.15g\n", qt_state.zeta_ns, qt_state.pok_ns);
		
		
		if(sav==2)
		{
			if(out[0]!=NULL)
			{
				fprintf(out[0],"%20.15g %20.15g %20.15g\n", log(c1.kpert),
				qt_state.zeta_ns,
				qt_state.pok_ns);
			}
		}
	
/*************************************************************************************
 * 
 * Classical Expanding phase 
 * 
 *************************************************************************************/
	
		NV_Ith_S (y0_cl, 0)		= Expo_mu_approx(qt_state.H, 1.0);
		NV_Ith_S (y0_cl, 1)		= qt_state.H;
		NV_Ith_S (yQ0_cl, 0)	= qt_state.phi;
		NV_Ith_S (y0_cl, 2)		= NV_Ith_S (y0_qt,2);
		NV_Ith_S (y0_cl, 3)		= NV_Ith_S (y0_qt,3);
		NV_Ith_S (y0_cl, 4)		= NV_Ith_S (y0_qt,4);
		NV_Ith_S (y0_cl, 5)		= NV_Ith_S (y0_qt,5);
		alpha_i 				= qt_state.alpha;
		
		
		//~ printf ("# PreCLE %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n",
				//~ alpha_i,
				//~ NV_Ith_S (y0_cl, 0),
				//~ NV_Ith_S (y0_cl, 1),
				//~ NV_Ith_S (y0_cl, 2), 
				//~ NV_Ith_S (y0_cl, 3),
				//~ NV_Ith_S (y0_cl, 4), 
				//~ NV_Ith_S (y0_cl, 5));
		
		cl_state = run_evol_cl (cvode_cl_p, alpha_i, 1.0e3, y0_cl, yQ0_cl, &c1, verb, sav, out);
		
		if(sav==2)
		{
			if(out[0]!=NULL)
			{
				fprintf(out[0],"%20.15g %20.15g %20.15g\n", log(c1.kpert),
				cl_state.zeta_ns,
				cl_state.pok_ns);
			}
			
			if(out[0]!=NULL)
			{
				fprintf(out[1],"tau_ns % 20.15g\n w_ns % 20.15g\n ",
					cl_state.tau_ns,
					cl_state.w_ns);
			}
		}
		
	}

}
//~fim do if
{ 
	if(out[0]!=NULL)
			{
				fprintf(out[1],"% 20.15g % 20.15g % 20.15g % 20.15g % 20.15Lg % 20.15Lg % 20.15Lg % 20.15Lg % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g % 20.15g",
					c1.d,
					c1.sigma,
					H_i,
					del_phi_i,
					cond1, 
					cond2,
					cond3,
					cond4,
					Calpha_i,
					Cphi_i,
					tau_p_d,
					tau_m_1,
					H_1, 
					delta_H_1,
					delta_x_1,
					c1.alpha_B,
					tau_m_2,
					H_2,			 
					c1.alpha0, 
					tau_0
					);
				}
}

  return (0);
}

int 
cvode_check_flag (void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  switch (opt)
  {
    case 0:
    {
      if (flagvalue == NULL)
      {
        fprintf (stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE;
      }
      break;
    }
    case 1:
    {
      errflag = (int *) flagvalue;
      if (*errflag < 0)
      {
        fprintf (stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
        return FALSE;
      }
      break;
    }
    case 2:
    {
      if (flagvalue == NULL)
      {
        fprintf (stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return FALSE;
      }
      break;
    }
    default:
    {
      fprintf (stderr, "Error: should not be reached!\n");
      exit (-1);
    }
  }
  return TRUE;
}

double 
sqrt1px_m1 (const double x)
{
  double binfact = 1.0;
  double res = 0;
  double xn = 1;
  int n = 0;

  if (fabs(x) > 1e-1)
    return sqrt (1.0 + x) - 1.0;

  while (TRUE)
  {
    binfact *= 3.0 / (2.0 * (1.0 + n++)) - 1.0;
    xn *= x;
    res += binfact * xn;

    if (fabs (binfact * xn / res) < DBL_EPSILON)
      break;
  }

  return res;
}
