/*http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_glimmix_sect023.htm*/
x 'cd \\cstatnas.cstat.msu.edu\redirects\zhangz19\Desktop\BUS';
LIBNAME RH '\\cstatnas.cstat.msu.edu\redirects\zhangz19\Desktop\BUS\';

ODS HTML CLOSE;         
FILENAME Results 'BUS_proportional_odds_model_TO.htm';
ODS HTML BODY=Results;
ODS PDF FILE='BUS_proportional_odds_model_TO.pdf'; 
TITLE1 "Proportional Odds Model with random effects for brand&nameplates data";

TITLE2 "Exploratory Analysis";
PROC FREQ DATA=dat;
  TABLES PI1_BRAND PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE;
  /*TABLES BAC11 NA1 NA2 NA11 NAMEPLATE; */
RUN;

TITLE2 "Proportional Odds Model with Random Effects";
proc glimmix data=dat method=LAPLACE; /*ASYCOV*/
      class PI_NAMEPLATE;
      model PI1_BRAND = PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE/ s cl link=CUMLOGIT  dist=multinomial; /*iB5_OWN*/
      random int / sub=PI_NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
	  covtest 0/ cl;
run;

/*proc glimmix data=dat method=LAPLACE ASYCOV;
      class NAMEPLATE;
      model BAC11 = NA1 NA2 NA11/ s cl link=CUMLOGIT  dist=multinomial; iB5_OWN
      random int / sub=NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
	  covtest 0/ cl;
run;*/

/*proc glimmix data=dat method=LAPLACE ASYCOV;
      class NAMEPLATE B5_OWN;
      model BAC11 = B5_OWN NA1 NA2 NA11 B5_OWN*NA11/ s cl link=CUMLOGIT  dist=multinomial;
      random int / sub=NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
	  covtest 0/ cl;
run;*/

/*proc glimmix data=dat method=LAPLACE;
      class PI_NAMEPLATE;
      model PI1_BRAND = PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE/ s cl link=CUMCLOGLOG  dist=multinomial;
      random int / sub=PI_NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
	  covtest 0/ cl;
run;*/

/*proc glimmix data=dat method=LAPLACE;
      class PI_NAMEPLATE;
      model PI1_BRAND = PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE/ s cl link=CUMLOGIT  dist=multinomial;
      random int PI2_EXPERIENCE/ sub=PI_NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
run;*/

/*proc glimmix data=dat method=quad;
      class PI_NAMEPLATE;
      model PI1_BRAND = PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE/ s cl link=cumprobit dist=multinomial;
      random int / sub=PI_NAMEPLATE s cl;
      ods output Solutionr=solr ParameterEstimates=solf;
	  covtest 0/ cl;
run;

TITLE2 "Random Effects Assessment";
ods select FitStatistics CovParms Covtests;
proc glimmix data=dat method=quad;
      class PI_NAMEPLATE;
      model PI1_BRAND = PI1_NAMEPLATE PI2_EXPERIENCE PI2_SERVICE/ link=cumprobit dist=multinomial;
      random int / sub=PI_NAMEPLATE;
	  covtest GLM;
run;
*/

/*proc sort data=solr; by Estimate;*/

/*data solr2; set solr; 
     length NAMEPLATE $16;
     obs  = _n_;
     NAMEPLATE = left(substr(Subject,11,16));
run;*/
/*
TITLE2 "Random Effects Ranking";
proc sgplot data=solr2;
   scatter x=obs y=estimate /
            markerchar  = PI_NAMEPLATE
            yerrorupper = upper
            yerrorlower = lower;
      xaxis grid label='PI_NAMEPLATE Rank' values=(1 9); 10 13 16  1 4 7
      yaxis grid label='Predicted PI_NAMEPLATE Effect';
run;*/
PROC EXPORT DATA= solr
            OUTFILE= "re_propodds.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
PROC EXPORT DATA= solf
            OUTFILE= "fe_propodds.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

ODS  _ALL_  CLOSE; /* Close all output destinations */
