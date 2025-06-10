/*Set up the library where you store the datasets*/
libname sas "C:\Users\henry\OneDrive\Desktop\Research\Programming\SAS\PS method workshop_2025\sample_output";


/***************************************************************************************************/
/*Step 1: Calculating the Propensity Score with imputed data*/
proc import datafile = "C:\Users\henry\Downloads\PS method\Data_dictionary_sample.xlsx"
	out = sas.cov dbms = xlsx replace;
	getnames = yes;
	sheet = "DEMO";
run;


proc import datafile = "C:\Users\henry\Downloads\PS method\sample_dataset.xlsx"
	out = sas.sample dbms = xlsx replace;
	getnames = yes;
	sheet = "sample";
run;


proc sql noprint; select distinct (variable_name) into :classcov separated by ' ' from sas.cov where Variable_type in ("binary"); quit;
%put &classcov.;


proc sql noprint; select distinct (variable_name) into :categocov separated by ' ' from sas.cov where Variable_type in ("binary") and Variable_name not in ("EXPOSE","MACE_OUTCOME"); quit;
%put &categocov.;


proc sql noprint; select distinct (variable_name) into :continucov separated by ' ' from sas.cov where Variable_type in ("continuous"); quit;
%put &continucov.;



/*Calculate Propensity score using logistic regression*/
proc logistic data=sas.Sample descending;
	class &classcov.;
	model EXPOSE= &categocov. &continucov.;
  output out = sas.Sample_ps (drop = _LEVEL_) p = denom;
run;


/***************************************************************************************************/
/*Step 2: check the distribution of PS before Matching/Weighting*/
data ps_before_weighting; 
	set sas.Sample_ps;
  exp_ps = .;
  non_ps= .;
  if EXPOSE = 1 then exp_ps = denom;
  if EXPOSE = 0 then non_ps = denom;
  keep exp_ps non_ps;
run;


proc sgplot data=ps_before_weighting;
   histogram exp_ps / transparency=0.75 fillattrs=(color="#0096A0") binwidth=0.05 legendlabel='Treatment A';
   density exp_ps / lineattrs=(color="#007B82" thickness=2) legendlabel='Treatment A';

   histogram non_ps / transparency=0.75 fillattrs=(color="#D85F58") binwidth=0.05 legendlabel='Treatment B';
   density non_ps / lineattrs=(color="#B5483D" thickness=2) legendlabel='Treatment B';

   keylegend / location=outside position=bottom;
   xaxis label="Propensity Score Distribution (before applying PS method)" values=(0 to 1.5 by 0.1);
run;


/***************************************************************************************************/
/**************************************/
/* Step 3: Apply different PS methods */
/**************************************/
/***************************************************************************************************/
/*Method 1:Covariate adjustment using PS*/
proc logistic data=sas.Sample_ps;
   class EXPOSE; 
   model MACE_OUTCOME = EXPOSE denom;
run;
/*Covariate adjustment using PS: OR 0.824(0.709-0.956)*/



/***************************************************************************************************/
/*Method 2:PS matching-Greedy (Nearest Neighbor) Matching*/
proc freq data=sas.Sample_ps;
	tables EXPOSE;
run;

/*For each treated unit, find the closest control (by PS distance)
that hasn’t yet been matched (and within caliper if specified)*/


proc psmatch data=sas.Sample region=cs;
   class &classcov.;
   psmodel EXPOSE= &categocov. &continucov.;
   match method=greedy(k=1) distance=lps caliper=0.2   
   /*a caliper width of 0.2 × SD of the logit(PS) is often recommended*/
   /*“logit PS” means take your propensity score p, form the odds p/(1-p),
   and then take the natural log of that odds.*/
         weight=none;
   assess lps allcov;
   output out(obs=match)=Final_dataset_ps_greedy_match lps=_Lps matchid=_MatchID;
run;


proc sort data=Final_dataset_ps_greedy_match out=Final_dataset_ps_greedy_match_m;
	by _MatchID;
run;


proc logistic data=Final_dataset_ps_greedy_match_m;
   class EXPOSE; 
   model MACE_OUTCOME = EXPOSE;
   strata _MatchID;
run;
/*Greedy matching: OR 0.853(0.723-1.006)*/


/***************************************************************************************************/
/*Method 3:Inverse probability of treatment weighting (IPTW)*/
/*Generate inverse propensity score weights */
proc logistic data=sas.Sample_ps descending;
	model EXPOSE=;
	output out=stable_ps p=num;
run;


data Final_dataset_ps_iptw;
	set  stable_ps;
	if EXPOSE=1 then uw=1/denom; else if EXPOSE=0 then uw=1/(1-denom);
	if EXPOSE=1 then sw=num/denom; else if EXPOSE=0 then sw=(1-num)/(1-denom);
run;


proc means data = Final_dataset_ps_iptw min p5 median mean p95 std max;
	var uw;
run;


proc means data = Final_dataset_ps_iptw min p5 median mean p95 std max;
	var sw;
run;


*** Fit logistic regression model without truncated weights ***;	
proc logistic data = Final_dataset_ps_iptw;
   	class EXPOSE(ref="0") MACE_OUTCOME(ref="0"); 
   	model MACE_OUTCOME = EXPOSE;
	weight uw;
run;
/*Unstabilized IP weighting without truncation: OR 0.789(0.721-0.862)*/


*** Truncate stabilized weights at the 1th and 99th percentile ***;
*Extract 1th and 99th percentile of stabilized weights;
proc univariate data = Final_dataset_ps_iptw noprint;
	var sw;
	output out=pctl pctlpts=1 99 pctlpre=p;
run;

* Save 5th, 95th percentile cutoff;
data temp2;
	set pctl;
	call symput ('cutoff_99', p99);
	call symput ('cutoff_1', p1);
run;


data Final_dataset_ps_sw_1_99;
	set Final_dataset_ps_iptw;
	sw_p1_99 = sw;
	if sw_p1_99 > %sysevalf(&cutoff_99) then do;
		sw_p1_99 = %sysevalf(&cutoff_99);
	end;
	else if sw_p1_99 < %sysevalf(&cutoff_1) then do;
		sw_p1_99 = %sysevalf(&cutoff_1);
	end;
run;


proc means data = Final_dataset_ps_sw_1_99 min p5 median mean p95 std max;
	var sw_p1_99;
run;


/*Check the distribution of PS after stabilized IP Weighting*/
data ps_after_sw; 
	set Final_dataset_ps_sw_1_99;
  exp_ps = .;
  non_ps= .;
  if EXPOSE = 1 then exp_ps = denom;
  if EXPOSE = 0 then non_ps = denom;
  keep exp_ps non_ps sw_p1_99;
run;


proc sgplot data=ps_after_sw;
   histogram exp_ps / weight= sw_p1_99 transparency=0.75 fillattrs=(color="#0096A0") binwidth=0.05 legendlabel='Treatment A';
   density exp_ps / weight= sw_p1_99 lineattrs=(color="#007B82" thickness=2) legendlabel='Treatment A';

   histogram non_ps / weight= sw_p1_99 transparency=0.75 fillattrs=(color="#D85F58") binwidth=0.05 legendlabel='Treatment B';
   density non_ps / weight= sw_p1_99 lineattrs=(color="#B5483D" thickness=2) legendlabel='Treatment B';

   keylegend / location=outside position=bottom;
   xaxis label="Propensity Score Distribution (after stabilized IP Weighting)" values=(0 to 1.5 by 0.1);
run;


*** Fit logistic regression model with truncated weights ***;	
proc logistic data = Final_dataset_ps_sw_1_99;
   	class EXPOSE(ref="0") MACE_OUTCOME(ref="0"); 
   	model MACE_OUTCOME= EXPOSE;
	weight sw_p1_99;
run;
/*Stabilized IP weighting with truncation: OR 0.793(0.695-0.904)*/



/***************************************************************************************************/
/*Method 4-1: Standardized mortality ratio weighting (SMRW)*/

data Final_dataset_ps_smrw;
	set sas.Sample_ps;
	if EXPOSE=1 then smrw=1; else if EXPOSE=0 then smrw=denom/(1-denom);
run;


proc means data = Final_dataset_ps_smrw min p5 median mean p95 std max;
	var smrw;
run;


*** Truncate weights at the 1th and 99th percentile ***;
*Extract 5th and 95th percentile of SMR weights;
proc univariate data = Final_dataset_ps_smrw noprint;
	var smrw;
	output out=pctl pctlpts=1 99 pctlpre=p;
run;

* Save 5th, 95th percentile cutoff;
data temp3;
	set pctl;
	call symput ('cutoff_smrw_99', p99);
	call symput ('cutoff_smrw_1', p1);
run;


*Truncate weights at the 1th and 99th percentile;
data Final_dataset_ps_smrw_1_99;
	set Final_dataset_ps_smrw;
	smrw_p1_99 = smrw;
	if smrw_p1_99 > %sysevalf(&cutoff_smrw_99) then do;
		smrw_p1_99 = %sysevalf(&cutoff_smrw_99);
	end;
	else if smrw_p1_99 < %sysevalf(&cutoff_smrw_1) then do;
		smrw_p1_99 = %sysevalf(&cutoff_smrw_1);
	end;
run;


proc means data = Final_dataset_ps_smrw_1_99 min p5 median mean p95 std max;
	var smrw_p1_99;
run;


*** Fit logistic regression model with truncated weights ***;	
proc logistic data = Final_dataset_ps_smrw_1_99;
   	class EXPOSE(ref="0") MACE_OUTCOME(ref="0");
   	model MACE_OUTCOME = EXPOSE;
	weight smrw_p1_99;
run;
/*Standardized mortality ratio weighting: OR 0.830(0.720-0.958)*/


 
