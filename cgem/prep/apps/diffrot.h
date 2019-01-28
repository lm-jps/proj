// Differential rotation coefficient

#define CARR_RATE   (2.86532908457)

// /cvs/JSOC/proj/lev1.5_hmi/libs/lev15/rotcoef_file.txt
//double diffrot[3] = {2.71390, -0.4050, -0.4220};

// Snodgrass '82/'84 coeffs from Rick Bogart's mtrack.c
double diffrot[3] = {-0.02893 + CARR_RATE, -0.3441, -0.5037};
