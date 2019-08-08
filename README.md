# PhD_ClaireVeillard
Contain matlab files from my PhD
This repository contains 6 functions and 3 classes:

Functions

1. matthewDW.m: calibration between T, d18Odol and d18Ow (Matthew et al., 1977)
2. RomanekcalcCO2.m: calibration between T, d13Ccalc and d13Cco2 (Romanek et al., 1992)
3. plot_isolines_d18Odol.m: Plot curves of constant values of d18Odol on a T=f(d18Ow) plot with the calibration of Matthew et al., 1977.
4. plot_isolines_d18Ow.m: Plot curves of constant values of d18Ow on a T=f(d18Odol) plot with the calibration of Matthew et al., 1977). 
5. daviesD2T.m: convert D47 to T with the calibration of Davies et al., 2018aa.
6. daviesT2D.m: convert T to D47 with the calibration of Davies et al., 2018aa.
    
Classes

1. cst.m: model constants
2. ic.m: user chooses the initial conditions
3. vectorcolor.m: a vector which contains hundreds of values between 0 and 1

Scripts
1. Recrystallization.m: Recrystallization model via dissolution reprecipitation
2. Recrystallization_IncTemp.m: Recrystallization model via dissolution reprecipitation at increasing temperature
3. FluidMixing.m: Dolomitization model from a mixture of fluid A and B
4. Cementation.m: Model of dolomite cementation
