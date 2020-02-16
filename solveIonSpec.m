% FUNCTION SOLVER FOR IONIC SPECIATION pH CALCULATION AND IONIC STRENGTH CORRECTION
% INPUTS
% Ic_ini: Initial Ic value to be use for convergence (if>0) and methods selector (if<0)
% pH_ini: Initial pH value to be used for pH loop convergence.
% StM:  Vector of total concentrations of components (typically from state variables)
% KdM:  Vector of equlibrium dehydration contants.
% KaM:  Vector of equlibrium acid-base contants.
% chrM: Matrix of species charge information.

% OUTPUTS
% pH:   Solved pH value based on proton activity
% Ic:   Solved ionic strength.
% Sh:   Solved concentration of protons in solution.
% spaM: Matrix of activities of all ionic species.
% spcM: Matrix of concentrations of all ionic species.

function [Ic, pH, Sh, spaM, spcM] = solveIonSpec(Ic_ini, pH_ini, StM, KdM, KaM, chrM, T, flag_solver) 
ah_ini = 10^(-pH_ini);                          % Activity of protons.
gammaM = ones(size(chrM));      gammaH = 1;     % initialisation of activity coefficients at 1 by default.
tolIc = 1e-10;  tolSh = 1e-15;                  % Convergence criteria paratemers.
niterSh = 0;  niterIc = 0;                      % Iterations counters
pH_limit = mod(pH_ini + 7, 14);            % This is the last value for convergence of pH.
errorIc = 1e4;    errorSh = 1e4;                % Initialising loop controllers.
maxIterSh = 25; maxIterIc = 25;                 % Iteration control parameters.
a_h2o = 1;                                      % Activity of water is one considering dilluted solution.
j=0;                                            % Counter of number of maxiter reached.


% This function can be called in three different forms. "flag_solver" name is used
% instead of "flag" because the latter is used when calling main SFunction.
switch flag_solver
    % Both pH and Ic are calcualted for the given system.
    case 0
    % Ionic strength Ic is constant at provided value Ic_ini.
    case 1
        % High Tol prevents Ic convergence and keeps provided value.
        tolIc = 1e3;
    % pH is constant at provided value pH_ini.
    case 2
        % High Tol prevents pH convergence and keeps provided value.
        tolSh = 1e3;
end

Ic_new = Ic_ini;
ah_new = ah_ini;

% Calculation of A and epsilon. They are parameters function of temperature for Gamma calculations.
    epsilon = 87.74 - 0.4008 * (T-273) + 9.398e-4 * (T-273)^2 - 1.410e-6 * (T-273)^3;
    A = 1.82e6 * (epsilon * T)^(-3/2);
     
% Convergence loop for ionic strength Ic.
while abs(errorIc)> tolIc
    niterIc = niterIc + 1;  % Iterations counter for Ic.
    Ic = Ic_new;
    %validpH has to be restarted as 0 to force the new Gammas to enter
    %again in inside the pH loop. The only exception that will not happen
    %is in the case we fix the pH (flag_solver = 2)
    if flag_solver ~= 2, validpH = 0; end 
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to be commented if gammas are to be used at 1.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of activity coefficients based on current Ic value.
    gammaM = 10.^(-A * chrM.^2 * (( sqrt(Ic)/(1+sqrt(Ic))   ) - 0.3*Ic ));
    gammaH = 10.^(-A * 1.^2    * (( sqrt(Ic)/(1+sqrt(Ic))   ) - 0.3*Ic ));
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iterative pH and ionic speciation solution for a given set of activity coefficients.
while abs(errorSh) > tolSh || validpH == 0
    niterSh = niterSh + 1;  % Iterations counter.
    % Activities Matrix initialization
    spaM = zeros(size(chrM));
    ah = ah_new;
    pH = -log10(ah);
    % Generalised common denominator for up to three deprotonations.
    DN = (KdM./(gammaM(:,1).*a_h2o))                      .* ah^3 ...
       + (1                                ./gammaM(:,2)) .* ah^3 ...
       + (KaM(:,1)                         ./gammaM(:,3)) .* ah^2 ...
       + (KaM(:,1) .* KaM(:,2)             ./gammaM(:,4)) .* ah   ...
       + (KaM(:,1) .* KaM(:,2) .* KaM(:,3) ./gammaM(:,5))         ;

    % Calculation of species activities function of the proton activity.
    spaM(:,1) = (StM .* ah^3 .* KdM)                             ./(DN * a_h2o);
    spaM(:,2) = (StM .* ah^3)                                    ./ DN;
    spaM(:,3) = (StM .* ah^2 .* KaM(:,1))                        ./ DN;
    spaM(:,4) = (StM .* ah   .* KaM(:,1) .* KaM(:,2))            ./ DN;
    spaM(:,5) = (StM         .* KaM(:,1) .* KaM(:,2) .* KaM(:,3))./ DN;

    % Calculation of species concentrations from activities.
    spcM = spaM ./ gammaM;
    Sh   = ah    / gammaH;

    % Evaluation of the charge balance error function for the current Sh value, F(Sh)
    F = Sh + sum(sum(spcM.*chrM)); 

    % Computing the derivative of the error function for Newton-Raphson algorithm.
    % Sizing of the species concentration derivatives matrix.
    dspcM = zeros(size(spcM));

    % Analytic derivative of the common denominator of the spcM expressions DN.
    dDN = (KdM ./ (gammaM(:,1).* a_h2o))       .* 3 * gammaH * ah^2 ...
       + (1                   ./  gammaM(:,2)) .* 3 * gammaH * ah^2 ...
       + KaM(:,1)             ./  gammaM(:,3)  .* 2 * gammaH * ah   ...
       + KaM(:,1) .* KaM(:,2) ./  gammaM(:,4)       * gammaH        ;
   
    % Generalised analytical derivative computation as per rule of (u/v)' = du/v - u*dv/v^2
    dspcM(:,1) = (1./(gammaM(:,1)*a_h2o)).*((3*StM .* gammaH .* ah^2 .* KdM)               ./ DN - (StM .* ah^3 .* KdM                            .* dDN) ./ (DN.^2));
    dspcM(:,2) = (1./gammaM(:,2))     .* ((3 * StM .* gammaH .* ah^2)                      ./ DN - (StM .* ah^3                                   .* dDN) ./ (DN.^2));
    dspcM(:,3) = (1./gammaM(:,3))     .* ((2 * StM .* gammaH .* ah .* KaM(:,1))            ./ DN - (StM .* ah^2 .* KaM(:,1)                       .* dDN) ./ (DN.^2));
    dspcM(:,4) = (1./gammaM(:,4))     .* ((    StM .* gammaH .*       KaM(:,1) .* KaM(:,2))./ DN - (StM .* ah   .* KaM(:,1) .* KaM(:,2)           .* dDN) ./ (DN.^2));
    dspcM(:,5) = (1./gammaM(:,5))     .* (                                                      - (StM         .* KaM(:,1) .* KaM(:,2) .* KaM(:,3).* dDN) ./ (DN.^2));

    % Derivative of the charge balance error respect to Sh dF(Sh)
    dF = 1 + sum(sum(dspcM.*chrM));
    % New Sh as per Newton-Raphson algorithm
    Sh_new = Sh - F/dF;
    ah_new = Sh_new * gammaH;    % New proton activity calculated from new proton concentration
    errorSh = F/dF;              % Difference between last to Sh values.

    % Checking if a valid pH is obtained
    if (ah > 1e-14) && (ah < 1),    
        validpH = 1;
    else
        validpH = 0;
    end

    % If number of iterations is exceeded then an new initial pH is used.
    if niterSh>=maxIterSh,
        % Counter of number of new pH_ini required.
        j = j + 1;
        pH_ini = mod(pH_ini + (-1)^j * j * 0.25, 14);
        %Restart of the errorSh and number of iterations for the new pH to test
        errorSh = 1;            niterSh = 0;        ah_new = 10^(-pH_ini);
        
        % If after checking all the pH values within the 0-14 range
        % we stop as no solution was found.
        if (pH_ini == pH_limit), 
            fprintf('\n> NO pH SOLUTION WAS REACHED. NO CONVERGENCE...\n');
            pH = -1
            break           
        end
    end
end


% Take into account the charge correction in the case that pH is fixed (flag_solver = 2 ).
% The correction is assumed for a monovalent species. If charge balance (F)
% is positive, the remaining charge is assigned to anions. If it is negative, the remaining
% charge is assigned to cations instead.
if flag_solver == 2, 
    Scat_id = find(strcmp('Scat',St.StNames));   San_id = find(strcmp('San', St.StNames)); % Consider moving from here
    if F < 0,
        spcM(Scat_id, 2) = spcM(Scat_id, 2) - F;
        spaM(Scat_id, 2) = spcM(Scat_id, 2) * gammaH;  %Assuming cations are monoprotic
    elseif F > 0,
        spcM(San_id, 2) = spcM(San_id, 2) + F;
        spaM(San_id, 2) = spcM(San_id, 2) * gammaH;  %Assuming anions are monoprotic
    end
end

% New value for ionic strength based on spcM for the given pH. It adds also
% the proton concentration. In case pH is fixed, the correction for the
% charge must be added
Ic_new = 0.5 * sum(sum(spcM.*chrM.^2)) + 0.5 * Sh_new * 1^2;
errorIc = Ic_new - Ic;  

%Checking the number of iterations of Ionic Strength
    if niterIc > maxIterIc, 
        fprintf('\n> NO Ic SOLUTION WAS REACHED. NO CONVERGENCE...\n');
        Ic = -1
        break        
    end
end


