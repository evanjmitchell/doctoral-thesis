# doctoral-thesis
This respository contains all of the Matlab code used in the thesis of Evan Mitchell, entitled "Coevolution of hosts and pathogens in the presence of multiple types of hosts", completed as part of the requirements for the degree Doctor of Philosophy at the University of Western Ontario. All of the code contained in this repository was written by EM.

### Files in chapter-2
- checkESSAlpha_VT.m --> function to check that the equilibrium values of alpha satisfy the evolutionary stability condition
- checkESSGamma_VT.m --> function to check that the equilibrium values of gamma satisfy the evolutionary stability condition
- findCSS_VT.m --> function to approximate the convergence stable evolutionary equilibrium
- findEE_VT.m --> function to approximate the endemic equilibrium of the resident system
- fitnessAlpha_VT.m --> function to compute the pathogen fitness
- fitnessGamma_VT.m --> function to compute the host fitness
- resident_VT.m --> function to define the system of differential equations governing the resident population
- sgradAlphaF_VT.m --> function to define the selection gradient for alpha in female hosts
- sgradAlphaM_VT.m --> function to define the selection gradient for alpha in male hosts
- sgradGammaF_VT.m --> function to define the selection gradient for gamma in female hosts
- sgradGammaM_VT.m --> function to define the selection gradient for gamma in male hosts
    
### Files in chapter-3
- b00.m --> function to compute the transmission rate beta_00
- b01.m --> function to compute the transmission rate beta_01
- b10.m --> function to compute the transmission rate beta_10
- b11.m --> function to compute the transmission rate beta_11
- dFitness.m --> function to compute the selection gradient
- findCSS.m --> function to approximate the CSS level of pathogen exploitation
- g.m --> function to compute the recovery rate gamma
- resident.m --> function to define the system of differential equations governing the resident population
