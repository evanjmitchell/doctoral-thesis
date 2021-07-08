function trans10 = b10(xi, epsil, bmax, kappa)
% A function to compute beta10 transmission rate
trans10 = bmax*(1 - epsil)*xi/(kappa + xi);