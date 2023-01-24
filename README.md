# Curesurv
Codes and examples for the application of the corrected mixture cure model in STATA and R
In this repository  you can find two folders. 
One named “curesurv_R” which includes:
o	A pdf file entitled "How to estimate a new mixture cure model for increased risk of non-cancer death" that illustrates the data analysis using the conventional and new models presented in the article and guides the user through the installation part of the R package.
o	A folder called “curesurv” includes the R package curesurv in tar.gz version. This package allows to fit a variety of cure models using the excess-risk modelling methodology, including the new mixed cure model.

The second, called 'curesurv_Stata', comprises:
o  The do file named  “script_ml_CorrectedMixtureCureModel “ that includes the programmes to  define the Maximum likelihood  to estimate model parameters using individual and grouped data. An example using the following two .dta files is also included in this script.
o A .dta file named “individual” that includes simulated individual data 
o A .dta file named “'lifetable” that includes the simulated lifetable of a general population 
o A .dta file named "extcuremod_estimates" that includes the results of the example for final comparison

Reference 
Botta, L., J. Goungounga, R. Capocaccia, G. Romain, M. Colonna, G. Gatta, O. Boussari, and V. Jooste
“A new cure model that corrects for increased risk of non-cancer death. Analysis of reliability and robustness, and application to real-life data.”  
