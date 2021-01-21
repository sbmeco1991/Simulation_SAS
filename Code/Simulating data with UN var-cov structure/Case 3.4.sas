%let homefolder=~/EST142/data;
libname bios6624 "&homefolder";

****************************************************************************
****************************************************************************
/* Simulate datasets from linear mixed models with a CS Var-Cov Structure */
@ Soumak Basumallik
****************************************************************************
****************************************************************************

**************************************************
/* Case 2.4, N = 200, Obs = 10, AR(1) Structure */
**************************************************;

%let NumIterations = 1000;  /* Number of datasets/iterations */
%let NumSubjects = 200;   /* Number of subjects (i) */
%let NumObs = 10;        /* Number of observations (j) */

/* Under the Alternate Hypothesis */

proc iml;
 cor={1.258948 1.216972 1.263182 1.284801 1.318455 1.346830 1.395617 1.472833 1.495551 1.511026,
      1.216972 1.569884 1.576514 1.620383 1.643137 1.657910 1.744566 1.798470 1.831529 1.821206,
      1.263182 1.576514 1.871055 1.850499 1.876640 1.904133 1.989870 2.021830 2.101136 2.127854,
      1.284801 1.620383 1.850499 2.184748 2.108270 2.137798 2.208366 2.282520 2.354493 2.381030,
      1.318455 1.643137 1.876640 2.108270 2.344584 2.309922 2.391588 2.454625 2.568323 2.638423,
      1.346830 1.657910 1.904133 2.137798 2.309922 2.593798 2.559027 2.611183 2.695155 2.777488,
      1.395617 1.744566 1.989870 2.208366 2.391588 2.559027 2.820401 2.834706 2.938639 3.025025,
      1.472833 1.798470 2.021830 2.282520 2.454625 2.611183 2.834706 3.220984 3.212444 3.307028,
      1.495551 1.831529 2.101136 2.354493 2.568323 2.695155 2.938639 3.212444 3.441534 3.462288,
      1.511026 1.821206 2.127854 2.381030 2.638423 2.777488 3.025025 3.307028 3.462288 3.694063};
 mean={0 0 0 0 0 0 0 0 0 0};
 N=&NumIterations*&NumSubjects;      /* NumIterations*NumSubjects */
  call randseed(61345);
  z=RandNormal(n,mean,cor);
  create USErrors from z;
  append from z;
quit;

data simulation_US;
 call streaminit(689127);
  do SampleID=1 to &NumIterations;
   do i=1 to &NumSubjects; 
     array e{10} col1-col10;
     set USErrors;
     if i<=(&NumSubjects/2) then Trt=0;
     else if i>(&NumSubjects/2) then Trt=1;
    do j=1 to &NumObs;  
      Time=j-1;
      Trt_Time=Trt*Time;
      Y = 0.5 + 1*Trt + 0.5*Time + 0.5*Trt_Time + e{j};   
     output; 
    end; 
  end;
end;

proc print data=simulation_US;
 where SampleID=1000;
run;

/* Under the Null Hypothesis */

proc iml;
 cor={1.258948 1.216972 1.263182 1.284801 1.318455 1.346830 1.395617 1.472833 1.495551 1.511026,
      1.216972 1.569884 1.576514 1.620383 1.643137 1.657910 1.744566 1.798470 1.831529 1.821206,
      1.263182 1.576514 1.871055 1.850499 1.876640 1.904133 1.989870 2.021830 2.101136 2.127854,
      1.284801 1.620383 1.850499 2.184748 2.108270 2.137798 2.208366 2.282520 2.354493 2.381030,
      1.318455 1.643137 1.876640 2.108270 2.344584 2.309922 2.391588 2.454625 2.568323 2.638423,
      1.346830 1.657910 1.904133 2.137798 2.309922 2.593798 2.559027 2.611183 2.695155 2.777488,
      1.395617 1.744566 1.989870 2.208366 2.391588 2.559027 2.820401 2.834706 2.938639 3.025025,
      1.472833 1.798470 2.021830 2.282520 2.454625 2.611183 2.834706 3.220984 3.212444 3.307028,
      1.495551 1.831529 2.101136 2.354493 2.568323 2.695155 2.938639 3.212444 3.441534 3.462288,
      1.511026 1.821206 2.127854 2.381030 2.638423 2.777488 3.025025 3.307028 3.462288 3.694063};
 mean={0 0 0 0 0 0 0 0 0 0};
 N=&NumIterations*&NumSubjects;      /* NumIterations*NumSubjects */
  call randseed(61345);
  z=RandNormal(n,mean,cor);
  create USErrors_Null from z;
  append from z;
quit;

data simulation_US_Null;
 call streaminit(689127);
  do SampleID=1 to &NumIterations;
   do i=1 to &NumSubjects; 
     array e{10} col1-col10;
     set USErrors_Null;
     if i<=(&NumSubjects/2) then Trt=0;
     else if i>(&NumSubjects/2) then Trt=1;
    do j=1 to &NumObs;  
      Time=j-1;
      Trt_Time=Trt*Time;
      Y = 0.5 + 1*Trt + 0.5*Time + 0*Trt_Time + e{j};   
     output; 
    end; 
  end;
end;

proc print data=simulation_US_Null;
 where SampleID=1000;
run;


*****************
/* MODEL FIT 1 */
*****************

************************************************************************
************************************************************************
/* Fitting a CS model based on the simulation obtained from Case 3.4 */
************************************************************************
************************************************************************;

***********************************************************************
/* Fixed Effects: Parameter Estimates, Standard Error, P-Value, Bias */
***********************************************************************;

ods graphics off;
ods exclude all; 

proc mixed data=simulation_US;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 repeated /type=CS subject=i;
 ods output SolutionF=FE (keep=Effect Estimate StdErr Probt);
run;

data FE_v2;
 set FE;
 if Effect="Intercept" then Bias=(0.5-Estimate);
 else if Effect="Trt" then Bias=(1-Estimate);
 else if Effect="Time" then Bias=(0.5-Estimate);
 else if Effect="Trt_Time" then Bias=(0.5-Estimate);
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_v2;
 by Effect;
run;

proc means data=FE_v2 mean;
 var Estimate StdErr Probt Bias;
 by Effect;
 ods output Summary=Results (keep=Effect Estimate_Mean StdErr_Mean Probt_Mean Bias_Mean);
run;

ods output close;
ods exclude none;

/* Renaming Variables in the Result dataset */

data Results_v2;
 set Results;
 rename Estimate_Mean = Average_Estimate;
 rename StdErr_Mean = Average_SE;
 rename Probt_Mean = Average_Pvalue;
 rename Bias_Mean = Average_Bias;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Results of Fixed Effects from the Simulations obtained under Case 3.4';
title3 'Fitting the model with a CS var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_CSFit_FE.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


*******************
/* Type II Error */
*******************;

/* In Type II error we are accepting a false Null hypothesis. Basically here we need the 
   proportion of p-values which are more than 0.5. */
  
proc freq data=FE_v2;
 tables P_less;
 by Effect;
 ods output OneWayFreqs=Results (keep=Effect P_less Frequency Percent);
run;
  
/* Renaming Variables in the Type II dataset */

data Results_v2;
 set Results;
 if P_less ne 1 then Type_II_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type II Error';
title3 'Fitting the model with a CS var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_CSFit_Type2.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


******************
/* Type I Error */
******************;

/* In Type I error we are rejecting a true Null hypothesis. Basically here we need the 
   proportion of p-values which are less than 0.5 under the Null */
  
/* Fitting the CS Model under the Null */

ods graphics off;
ods exclude all; 

proc mixed data=simulation_US_Null;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 repeated /type=CS subject=i;
 ods output SolutionF=FE_Null (keep=Effect Estimate StdErr Probt);
run;

data FE_Null_v2;
 set FE_Null;
 if Effect ne "Trt_Time" then delete;
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_Null_v2;
 by Effect;
run;

proc freq data=FE_Null_v2;
 tables P_less;
 by Effect;
 ods output OneWayFreqs=Results (keep=Effect P_less Frequency Percent);
run;
  
ods exclude none;
  
/* Renaming Variables in the Type I dataset */

data Results_v2;
 set Results;
 if P_less eq 1 then Type_I_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type I Error';
title3 'Fitting the model with a CS var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_CSFit_Type1.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


*****************
/* MODEL FIT 2 */
*****************

***************************************************************************
***************************************************************************
/* Fitting a Linear model based on the simulation obtained from Case 2.4 */
***************************************************************************
***************************************************************************;

***********************************************************************
/* Fixed Effects: Parameter Estimates, Standard Error, P-Value, Bias */
***********************************************************************;

ods graphics off;
ods exclude all; 

proc glm data=simulation_US;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 ods output ParameterEstimates=FE (keep=Parameter Estimate StdErr Probt);
run;

data FE_v2;
 set FE;
 if Parameter="Intercept" then Bias=(0.5-Estimate);
 else if Parameter="Trt" then Bias=(1-Estimate);
 else if Parameter="Time" then Bias=(0.5-Estimate);
 else if Parameter="Trt_Time" then Bias=(0.5-Estimate);
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_v2;
 by Parameter;
run;

proc means data=FE_v2 mean;
 var Estimate StdErr Probt Bias;
 by Parameter;
 ods output Summary=Results (keep=Parameter Estimate_Mean StdErr_Mean Probt_Mean Bias_Mean);
run;

ods output close;
ods exclude none;

/* Renaming Variables in the Result dataset */

data Results_v2;
 set Results;
 rename Estimate_Mean = Average_Estimate;
 rename StdErr_Mean = Average_SE;
 rename Probt_Mean = Average_Pvalue;
 rename Bias_Mean = Average_Bias;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Results of Fixed Effects from the Simulations obtained under Case 3.4';
title3 'Fitting the model with an Uncorrelated var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_UNCORRFit_FE.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


*******************
/* Type II Error */
*******************;

/* In Type II error we are accepting a false Null hypothesis. Basically here we need the 
   proportion of p-values which are more than 0.5. */
  
proc freq data=FE_v2;
 tables P_less;
 by Parameter;
 ods output OneWayFreqs=Results (keep=Parameter P_less Frequency Percent);
run;
  
/* Renaming Variables in the Type II dataset */

data Results_v2;
 set Results;
 if P_less ne 1 then Type_II_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type II Error';
title3 'Fitting the model with an Uncorrelated var-cov Structure';

ods rtf file='~/EST142/data/Tab_3.4_UNCORRFit_Type2.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


******************
/* Type I Error */
******************;

/* In Type I error we are rejecting a true Null hypothesis. Basically here we need the 
   proportion of p-values which are less than 0.5 under the Null */
  
/* Fitting the Uncorrelated Model under the Null */

ods graphics off;
ods exclude all; 

proc glm data=simulation_US_Null;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 ods output ParameterEstimates=FE_Null (keep=Parameter Estimate StdErr Probt);
run;

data FE_Null_v2;
 set FE_Null;
 if Parameter ne "Trt_Time" then delete;
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_Null_v2;
 by Parameter;
run;

proc freq data=FE_Null_v2;
 tables P_less;
 by Parameter;
 ods output OneWayFreqs=Results (keep=Parameter P_less Frequency Percent);
run;
  
ods exclude none;
  
/* Renaming Variables in the Type I dataset */

data Results_v2;
 set Results;
 if P_less eq 1 then Type_I_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type I Error';
title3 'Fitting the model with a Uncorrelated var-cov Structure';

ods rtf file='~/EST142/data/Tab_3.4_UNCORRFit_Type1.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;



*****************
/* MODEL FIT 3 */
*****************

************************************************************************
************************************************************************
/* Fitting a US model based on the simulation obtained from Case 2.4 */
************************************************************************
************************************************************************;

***********************************************************************
/* Fixed Effects: Parameter Estimates, Standard Error, P-Value, Bias */
***********************************************************************;

ods graphics off;
ods exclude all; 

proc mixed data=simulation_US;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 repeated /type=UN subject=i;
 ods output SolutionF=FE (keep=Effect Estimate StdErr Probt);
run;

data FE_v2;
 set FE;
 if Effect="Intercept" then Bias=(0.5-Estimate);
 else if Effect="Trt" then Bias=(1-Estimate);
 else if Effect="Time" then Bias=(0.5-Estimate);
 else if Effect="Trt_Time" then Bias=(0.5-Estimate);
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_v2;
 by Effect;
run;

proc means data=FE_v2 mean;
 var Estimate StdErr Probt Bias;
 by Effect;
 ods output Summary=Results (keep=Effect Estimate_Mean StdErr_Mean Probt_Mean Bias_Mean);
run;

ods output close;
ods exclude none;

/* Renaming Variables in the Result dataset */

data Results_v2;
 set Results;
 rename Estimate_Mean = Average_Estimate;
 rename StdErr_Mean = Average_SE;
 rename Probt_Mean = Average_Pvalue;
 rename Bias_Mean = Average_Bias;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Results of Fixed Effects from the Simulations obtained under Case 3.4';
title3 'Fitting the model with an US var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_USFit_FE.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


*******************
/* Type II Error */
*******************;

/* In Type II error we are accepting a false Null hypothesis. Basically here we need the 
   proportion of p-values which are more than 0.5. */
  
proc freq data=FE_v2;
 tables P_less;
 by Effect;
 ods output OneWayFreqs=Results (keep=Effect P_less Frequency Percent);
run;
  
/* Renaming Variables in the Type II dataset */

data Results_v2;
 set Results;
 if P_less ne 1 then Type_II_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type II Error';
title3 'Fitting the model with an US var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_USFit_Type2.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;


******************
/* Type I Error */
******************;

/* In Type I error we are rejecting a true Null hypothesis. Basically here we need the 
   proportion of p-values which are less than 0.5 under the Null */
  
/* Fitting the US Model under the Null */

ods graphics off;
ods exclude all; 

proc mixed data=simulation_US_Null;
 by SampleID;
 class i;
 model y=Trt Time Trt_Time/ solution;
 repeated /type=UN subject=i;
 ods output SolutionF=FE_Null (keep=Effect Estimate StdErr Probt);
run;

data FE_Null_v2;
 set FE_Null;
 if Effect ne "Trt_Time" then delete;
 if Probt <= 0.05 then P_less = 1;
 else if Probt > 0.05 then P_less = 0;
run;
 
proc sort data=FE_Null_v2;
 by Effect;
run;

proc freq data=FE_Null_v2;
 tables P_less;
 by Effect;
 ods output OneWayFreqs=Results (keep=Effect P_less Frequency Percent);
run;
  
ods exclude none;
  
/* Renaming Variables in the Type I dataset */

data Results_v2;
 set Results;
 if P_less eq 1 then Type_I_Error=Percent;
run;

options nocenter nodate nonumber nolabel;
title1 'Table 3.4';
title2 'Type I Error';
title3 'Fitting the model with an US var-cov structure';

ods rtf file='~/EST142/data/Tab_3.4_USFit_Type1.rtf' bodytitle_aux style=Journal3;

proc print data=Results_v2;
run;

ods rtf close;








