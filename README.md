# Slow_Tight
Topic:
-Biochemistry
-Enzymology
-Thrombosis and hemostasis.
This repository contains an R script that allows the correction of enzymatic substrate conversion for substrate consumption, recorded in time. 
The curves are then fitted for slow-tight binding. The script returns the corrected curves, the fitted curves, the parameters for substrate correction, the kinetic parameters V0, Vs and kobs, and the percentage of enzyme activity versus time. 
All files are saved in time-labelled .csv formats in folders labelled as date and can be readily imported in spreadsheet software (provided that the locale is US or UK, ie. periods as decimal separators). 
The script was designed to correct the inactivation of factor Xa by TFPI, monitored by chromogenic substrate conversion (S2765, CS-1165), but can be used for correction of substrate conversion of any enzyme and for fitting slow-tight binding of any enzyme-inhibitor pair that obeys these kinetics.

The data can be imported from the clipboard. Importante considerations are:
- There is a script for Mac OSX and for Windows. This is because the clipboard is addressed differently in these environments.
- The first column of the data always contains the substrate conversion curve on which all corrections are based. This generally means that the enzymatic substrate conversion is recorded in the absence of an inhibitor but else under the same conditions as all other recorded curves.
- Data should be "rectangular", meaning only the named columns with the respective readings. No notes, dates, filenames.
- The titles or annotations are always in the first row of the data-array. These will be converted into syntactically correct titles for R.
- The script corrects for empty cells, but be aware that all readings will be truncated to the length of the shortest measurement. This is because the used R data.frame class can only contain vectors with the same length.


