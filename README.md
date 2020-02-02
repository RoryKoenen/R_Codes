# Slow_Tight
Topic:
-Biochemistry
-Enzymology
-Thrombosis and hemostasis
This repository contains an R script that allows the correction of enzymatic substrate conversion in time for substrate consumption. Subsequently, the curves are fitted for slow-tight binding. The script returns the corrected curves, the fitted curves, the parameters for substrate correction, the kinetic parameters V0, Vs and kobs, and the percentage of enzyme activity versus time. All files are saved in time-labelled .csv formats in folders labelled as date and can be readily imported in spreadsheet software (provided that the locale is US or UK, ie. periods as decimal separators). 
The script was designed to correct the inactivation of factor Xa by TFPI, monitored by chromogenic substrate conversion (S2765, I1165), but can be used for correction of substrate conversion of any enzyme and for fitting slow-tight binding of any enzyme-inhibitor pair that obeys these kinetics.
