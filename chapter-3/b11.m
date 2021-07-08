function trans11 = b11(xi, epsil, bmax, kappa)
% A function to compute beta11 transmission rate
trans11 = bmax*(1 - epsil)^2*xi/(kappa + xi);