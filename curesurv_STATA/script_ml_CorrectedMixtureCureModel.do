*---------------------------------Maximum likelihood estimation using individual and grouped data-----------------------
/*
 To estimate the parameter using Maximum likelihood on individual data run the lines 10-106
 To estimate the parameter using Maximum likelihood on grouped data run the lines 23-67, 114-196 and 198-200
 To verify that both methods work run all the code
*/

*DEFINITION OF THE PROGRAMS

program define extcuremod_individual
version 14.0
args lnf lambda gamma pcure alpha delta
tempvar usurv uhaz cured
quietly replace cured= 1/(1+exp(-`pcure'))
quietly replace usurv=(exp(-`lambda'*($ML_y1)^`gamma'))^exp(-`delta')
quietly replace uhaz=((1-cured)*exp(-`delta')*`lambda'*`gamma'*(($ML_y1)^(`gamma'-1))*usurv) / (cured + (1-cured)*usurv)
quietly replace `lnf'= $ML_y2*ln($ML_y3 * `alpha'+ uhaz)+ln(cured + (1-cured)*usurv) +ln($ML_y4)*`alpha'
end


****   User specific path  ****
* CHANGE IT ACCORDING TO YOUR SETTING, individual and lifetable files have to be saved in the directory listed in cd   
adopath ++ "C:\ado"
cd "C:\Users"


pause on

use "individual.dta"
/* 
Variables included in individual
i= progressive indicator to identify each patient
agediag=age in years
tfup=time to event 
lifestat= failure event 
exphazard= expected hazard of death at attained age
expsurv= expected survival at attained age
*/

*---------------------------------Maximum likelihood estimation using individual data-----------------------

generate meth="Individual"
generate censored_time=15
gen agestand=(agediag-60)/15	
gen double uhaz=0
gen double usurv=0
gen double cured=0
generate likesave=.
generate alfa_mle=.
generate cured60_mle=.
generate scale_mle= .
generate shape_mle= .
generate beta_mle=.
generate delta_mle= .
generate alfase_mle=.
generate cured60se_mle= .
generate scalese_mle= .
generate shapese_mle= .
generate betase_mle= .
generate deltase_mle= .
generate cured_mle=.
generate curedse_mle= .	

* Stset the individual data

stset tfup, failure (lifestat==1) id(i) exit(time censored_time)

**** saving dataset ****

save "temp.dta" , replace

****************** MAXIMUM LIKELIHOOD ESTIMATION_Individual data  ****

ml model lf extcuremod_individual (lambda: _t _d exphaz expsurv  = ) (gamma: ) (pcure: agestand) (alpha: )  (delta:agestand, nocons) 
ml init 0.8  1.  0.3 0.0 1.0 0.0 ,copy
ml maximize 

pause type q to continue

****   Postestimation values  ****

replace likesave =e(ll)
replace alfa_mle=_b[alpha:_cons] 
replace cured_mle=_b[pcure:_cons] 
replace cured60_mle=(1/(1+exp(-cured_mle))) 
replace scale_mle=_b[lambda:_cons]  
replace shape_mle=_b[gamma:_cons] 
replace beta_mle=_b[pcure:agestand] 
replace delta_mle=_b[delta:agestand] 

replace alfase_mle=_se[alpha:_cons] 
replace curedse_mle=_se[pcure:_cons] 
replace cured60se_mle= curedse_mle*(exp(cured_mle))/(1+exp(cured_mle))^2 
replace scalese_mle=_se[lambda:_cons] 
replace shapese_mle=_se[gamma:_cons] 
replace  betase_mle=_se[pcure:agestand] 
replace deltase_mle=_se[delta:agestand] 

****  Exporting estimates   ****
	
	keep meth likesave *mle*
	duplicates drop *, force
	save "extcuremod_estimates.dta" , replace
	
pause off	

***************************************************GROUPED DATA 
****  re-importing the dataset  *****

use "temp.dta" , clear

*DEFINITION OF THE PROGRAM 
program define extcuremod_grouped
version 14.0
args lnf lambda gamma pcure alpha delta
tempvar usurv1 usurv0  cured surv1 surv0
quietly replace  cured= 1/(1+exp(-`pcure'))
quietly replace  usurv1=(exp(-`lambda'*($ML_y1)^`gamma'))^exp(-`delta')
quietly replace  usurv0=(exp(-`lambda'*($ML_y1 - 1)^`gamma'))^exp(-`delta')
quietly replace  surv1=cured+ (1-cured)*usurv1
quietly replace  surv0=cured+ (1-cured)*usurv0
quietly replace `lnf'= $ML_y2*ln(1-($ML_y4 ^ `alpha')*surv1/surv0)+ ($ML_y3-$ML_y2-0.5*$ML_y5)*(ln(surv1) - ln(surv0) + ln($ML_y4)*`alpha') 
end

pause on
*Generate Year since diagnosis to correspond to _year in your lifetable. In this example year_dg=1 is used

generate year_dg=1

*Generate the relative survival lifetable using the strs comand.
strs using "lifetable",	br(0(1)15) mergeby( _year _age) diagage(agediag) diagyear(year_dg) maxage(100) by(agestand) list(n d w r cp cp_e2 cr_e2 se_cr_e2 d_star) 	savgroup(grouped_sample, replace)

*use the grouped_sample file
use "grouped_sample.dta", clear

gen double usurv0=0
gen double usurv1=0
gen double surv0=0
gen double surv1=0
gen double cured=0
generate fup=start+1	
generate meth="Grouped"
generate likesave=.
generate alfa_mle=.
generate cured60_mle=.
generate scale_mle= .
generate shape_mle= .
generate beta_mle=.
generate delta_mle= .
generate alfase_mle=.
generate cured60se_mle= .
generate scalese_mle= .
generate shapese_mle= .
generate betase_mle= .
generate deltase_mle= .
generate cured_mle=.
generate curedse_mle= .	

****  MAXIMUM LIKELIHOOD ESTIMATION_Grouped data  ****
/*DESCRIPTION OF THE VARIBLES USED
d= observed number of deaths during the interval
p_star= intreval specific expected survival
fup=end of the observed interval 
n=number of cases alive at the begining of the interval
w= withdrawals during the interval
*/
ml model lf extcuremod_grouped (lambda: fup d n p_star w = ) (gamma: ) (pcure: agestand) (alpha: )  (delta:agestand, nocons) 
* set the initial values
ml init 0.8  1.  -0.0 0.0 1.0 -0.0 ,copy
ml maximize

pause type q to continue

****   Postestimation values  ****

replace likesave =e(ll)
replace alfa_mle=_b[alpha:_cons] 
replace cured_mle=_b[pcure:_cons] 
replace cured60_mle=(1/(1+exp(-cured_mle))) 
replace scale_mle=_b[lambda:_cons]  
replace shape_mle=_b[gamma:_cons] 
replace beta_mle=_b[pcure:agestand] 
replace delta_mle=_b[delta:agestand] 

replace alfase_mle=_se[alpha:_cons] 
replace curedse_mle=_se[pcure:_cons] 
replace cured60se_mle= curedse_mle*(exp(cured_mle))/(1+exp(cured_mle))^2 
replace scalese_mle=_se[lambda:_cons] 
replace shapese_mle=_se[gamma:_cons] 
replace  betase_mle=_se[pcure:agestand] 
replace deltase_mle=_se[delta:agestand] 

****  Exporting estimates   ****

	keep meth likesave *mle* 
	duplicates drop *, force
	append using  "extcuremod_estimates.dta" 
	save  "extcuremod_estimates.dta"  , replace
pause off
exit