x 'cd \\cstatnas.cstat.msu.edu\redirects\zhangz19\Desktop\BUS';
LIBNAME RH '\\cstatnas.cstat.msu.edu\redirects\zhangz19\Desktop\BUS\';
ODS HTML CLOSE;         
FILENAME Results 'BUS_Beta_Mixed_Regression_2013-06-24.htm';
ODS HTML BODY=Results;
ODS PDF FILE='BUS_Beta_Mixed_Regression_2013-06-24.pdf'; 
TITLE1 "Beta Mixed Regression for Ford brand&nameplates data";

/*===========================================================================*/
/* DEFINE MACROS */

/* Beta Regression Residual Plots */ 

%Macro Beta_Regression_Residuals (mu, phi, depvar);

data mu_results;
    set &mu;
	resid = &depvar - pred;
    if &depvar = . then mu_hat = .;
       else mu_hat = pred;
    record=_n_;
    keep record &depvar mu_hat resid;
run;

proc sgplot data=mu_results;
	scatter x=mu_hat y=resid;
run;

data phi_results;
    set &phi;
    if &depvar = . then phi_hat = .;
       else phi_hat = pred;
    record=_n_;
    keep record &depvar phi_hat;
run;

data residual;
    merge mu_results phi_results;
    by record;
    pearson = (&depvar - mu_hat) / sqrt( (mu_hat * (1-mu_hat)) / (phi_hat + 1) );
    mustar = digamma(mu_hat * phi_hat) - digamma( (1-mu_hat) * phi_hat);
    astar = trigamma(mu_hat * phi_hat) + trigamma( (1-mu_hat) * phi_hat);
    logit_y = log(&depvar / (1 - &depvar) );
    scoreres = ( logit_y - mustar) / sqrt(astar);
run;

proc rank data=residual out=qqplot normal=blom;
    var pearson scoreres;
    ranks pearsonrank scorerank;
run;

proc sql;
    Create table Range As
    select ceil(abs(max(pearson))) as Max1,
           ceil(abs(min(pearson))) as Min1,
           ceil(abs(max(scoreres))) as Max2,
           ceil(abs(min(scoreres))) as Min2
    from qqplot ;

Data Range;
    set Range;
    Range1=Max(Max1,Min1);
    Range2=Max(Max2,Min2);
    call symput("Range1",Range1);
    call symput("Range2",Range2);
run;

ODS GRAPHICS / HEIGHT = 7IN WIDTH = 7IN;

PROC SGPLOT DATA = qqplot;
    XAXIS LABEL = "Inverse Normal Distribution" MIN = -&Range1 MAX = &Range1;
	YAXIS LABEL = "Pearson Residual" MIN = -&Range1 MAX = &Range1;
	SCATTER X = pearsonrank Y = pearson;
	SERIES X = pearsonrank Y = pearsonrank / 
      LINEATTRS= GraphPrediction (PATTERN = 2 COLOR = RED);
RUN;

PROC SGPLOT DATA = qqplot;
    XAXIS LABEL = "Inverse Normal Distribution" MIN = -&Range2 MAX = &Range2;
	YAXIS LABEL = "Standardized Residual" MIN = -&Range2 MAX = &Range2;
	SCATTER X = scorerank Y = scoreres;
	SERIES X = scorerank  Y = scorerank  / 
      LINEATTRS= GraphPrediction (PATTERN = 2 COLOR = RED);
RUN;

ODS GRAPHICS / RESET = all;

%mend; /* End of Macro Beta_Regression_Residuals*/

TITLE2 "Beta regression, logit link, random intercept, location & dispersion submodels"; 
TITLE3 "Run the mixed beta regression"; 
/* The input data set must be clustered according to the SUBJECT= variable*/
proc sort data=dat; by PI_NAMEPLATE;
run;

/*===========================================================================*/
/* MODEL 5: MIXED BETA REGRESSION MODEL W/LOGIT LINK, RANDOM INTERCEPT, & 
   SIMPLE DISPERSION MODEL. */
PROC NLMIXED DATA = dat METHOD=GAUSS QPOINTS = 20;
PARMS /* Fixed Effects for Location Submodel (estimates from model 3 GLM output) */
      b0 = 0, b1 = 0, b2 = 0, b3 = 0,
	  g0 = 0, g1 = 0, g2 = 0, g3 = 0,
      su0 = 1;      
xb = (b0+u0) + b1*PI1_NAMEPLATE + b2*PI2_EXPERIENCE + b3*PI2_SERVICE;  
mu = EXP(xb)/(1 + EXP(xb));          
wgamma = g0 + g1*PI1_NAMEPLATE + g2*PI2_EXPERIENCE  + g3*PI2_SERVICE;  
phi = exp(wgamma);
var_y = (mu*(1 - mu))/(1 + phi);
p = mu*phi;
q = phi - mu*phi;
ll = lgamma(p+q) - lgamma(p) - lgamma(q) + (p-1)*log(PI1_BRANDtrans) + (q-1)*log(1 - PI1_BRANDtrans); 
MODEL PI1_BRANDtrans ~ GENERAL(ll);      
RANDOM u0 ~ NORMAL(0, su0) SUBJECT = PI_NAMEPLATE out = solr;      
PREDICT mu OUT = m5pred_mu;   
PREDICT phi OUT = m5pred_phi; 
PREDICT var_y OUT = m5pred_vary; 
RUN;

/* Residual Plots*/ 
%Beta_Regression_Residuals(m5pred_mu, m5pred_phi, PI1_BRANDtrans);

PROC CORR DATA = m5pred_mu;
VAR PI1_BRANDtrans pred;
RUN;

/*
PROC MEANS DATA=m5pred_mu;
  VAR Voltage MA PMA TMA LogitTMA BMI CentBMI pred;
RUN;

PROC EXPORT DATA= WORK.m5pred_mu 
            OUTFILE= "RCPMa_Model5_Predicted_Values_2013-04-24.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
*/

TITLE2 "Random Effects Ranking";
proc sort data=solr; by Estimate;
data solr2; set solr; 
     obs  = _n_;
run;
proc sgplot data=solr2;
   scatter x=obs y=estimate /
            markerchar  = PI_NAMEPLATE
            yerrorupper = upper
            yerrorlower = lower;
      xaxis grid label='PI_NAMEPLATE Rank' values=(1 4 7 10 13 16);
      yaxis grid label='Predicted PI_NAMEPLATE Effect';
run;
ODS  _ALL_  CLOSE; /* Close all output destinations */
