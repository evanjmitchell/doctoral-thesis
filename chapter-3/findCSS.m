function css = findCSS(xires, tol, step, maxit, k, ...
    epsil, bmax, cost, kappa)
% A function to find a candidate CSS exploitation level; 
% input xires should be a vector of the lower and upper
% initial resident exploitation levels; step is the step 
% size used for choosing a new exploitation level and 
% should be bigger than tol, the tolerance level for 
% ending the program; maxit is the maximum number of 
% iterations to compute before ending the program; other 
% parameters are as detailed in the model in 
% the main paper.
r0 = bmax/(1 + sqrt(kappa))^2;
usir = 1/r0;
vsir = (1/(1 + sqrt(kappa)))*(1 - 1/r0);
for i = 1:2
    % find equilibrium using Euler's method
    x0 = [usir, vsir, 0.01, 0.02]; % initialize
    dxidt = dFitness(xires(i), x0(1), x0(3), k, ...
        epsil, bmax, cost, kappa);
    itn = 1;
    if abs(dxidt) < tol
        flagxi = 0;
    else
        flagxi = 1;
    end
    while flagxi == 1
        if itn <= maxit % check to make sure maxit hasn't 
        % been exceeded
            flagx0 = 1; % flag to indicate resident system is 
            % away from equilibrium
            while flagx0 == 1
                dx0dt = resident(x0, xires(i), k, ...
                        epsil, bmax, cost, kappa);
                if max(abs(dx0dt)) < 1e-03
                    flagx0 = 0;
                end
                x0 = x0 + (1e-05)*dx0dt;
                for j = 1:4
                    if x0(j) < 0
                        x0(j) = 0;
                    end
                end
                if x0(1) + x0(2) > 1
                    x0(1) = 0.5*(1 + x0(1) - x0(2));
                    x0(2) = 1 - x0(1);
                    x0(1) = 0.9*x0(1);
                    x0(2) = 0.9*x0(2);
                end
                for j = 3:4
                    if x0(j) > 1
                        x0(j) = 1 - 1e-05;
                    end
                end
            end
            dxidt = dFitness(xires(i), x0(1), x0(3), ...
                k, epsil, bmax, cost, kappa);
            if abs(dxidt) < tol
                flagxi = 0;
            end
            xires(i) = xires(i) + step*dxidt;
            itn = itn + 1; % increase iteration counter
        else
            % print message if maximum number of iterations 
            % has been exceeded
            out = sprintf('%f ', xires);
            fprintf(['Maximum number of iterations ', ...
                'exceeded; increase maxit and try ', ...
                'again. Most recent xi values: %s.\n'], ...
                out); % print message to console
            fprintf('Number of iterations: %d.\n', itn);
            return
        end
    end
end
% return average of both trajectories and take steps on either side 
% to approximate second derivative
css = [mean(xires) + step, mean(xires), mean(xires) - step];
% compute equilibrium values of u, v, x, and y using Euler's method
x0 = [[usir, vsir, 0.01, 0.02]; [usir, vsir, 0.01, 0.02]; ...
    [usir, vsir, 0.01, 0.02]]; % initialize
for i = 1:3
    flagx0 = 1; % flag to indicate resident system is 
    % away from equilibrium
    while flagx0 == 1
        dx0dt = resident(x0(i, :), css(i), k, epsil, ...
            bmax, cost, kappa);
        if max(abs(dx0dt)) < 1e-03
            flagx0 = 0;
        end
        x0(i, :) = x0(i, :) + (1e-05)*dx0dt;
        for j = 1:4
            if x0(i, j) < 0
                x0(i, j) = 0;
            end
        end
        if x0(i, 1) + x0(i, 2) > 1
            x0(i, 1) = 0.5*(1 + x0(i, 1) - x0(i, 2));
            x0(i, 2) = 1 - x0(i, 1);
            x0(i, 1) = x0(i, 1) - 1e-05;
            x0(i, 2) = x0(i, 2) - 1e-05;
        end
        for j = 3:4
            if x0(i, j) > 1
                x0(i, j) = 1 - 1e-05;
            end
        end
    end
end
% evaluate payoff function at these points
dWcssplus = dFitness(css(1), x0(1, 1), x0(1, 3), k, epsil, ...
    bmax, cost, kappa);
dWcssminus = dFitness(css(3), x0(3, 1), x0(3, 3), k, epsil, ...
    bmax, cost, kappa);
% use these values to evaluate the ESS condition
essCondition = (dWcssplus - dWcssminus)/(2*step);
% save the equilibrium values and CSS value, along with an 
% indicator of whether the ESS condition is satisfied
if essCondition < (step^2)*10
    css = [x0(2, :), css(2), 1, essCondition];
else
    css = [x0(2, :), css(2), 0, essCondition];
end
out = sprintf('%f ', css);
fprintf(['Equilibrium values of u, v, x, y at the ', ...
    'CSS, the CSS level of xi, and whether the ESS condition ', ...
    'is satisfied (with value of second derivative): %s.\n'], ...
    out);