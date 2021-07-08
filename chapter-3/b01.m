function trans01 = b01(xi, epsil, bmax, kappa)
% A function to compute beta01 transmission rate
trans01 = bmax*(1 - epsil)*xi/(kappa + xi);