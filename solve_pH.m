% Typically charge balance solution for pH calculation loop
% Function dealing with the numerical solution of the charge balance by means of the Newton method.
function [fSh, fspcM] = solve_pH(pHini, StV, KdM, KaM, chrM)
% Initialising variables.
pH_last = 100;  
% Counter of convergences & Break indicator.
j=1;    brk = 0;		
% Main loop to be repeated until a valid pH value is obtained.
while brk==0,
    pH=pHini;   F = 1;  i=1;    relTol = 1.e-5;     maxIter = 30;
% Setting convergence and break criteria.
while (abs((pH-pH_last)/pH_last)>relTol) && (i<maxIter) && (isnan(abs((pH-pH_last)/pH_last))==0),
    % First the F value for x+dx for numeric derivative computation.
    pH_orig = pH;       eps = 0.0001;    pH = pH*(1+eps);
    % Concentration for the Sh plus delta Sh for numerical derivative calculation.
    Sh = 10^(-pH);
    spcM = zeros(size(chrM));
    DNM = ( KaM(:,1).*Sh^2  +   KaM(:,1).*KaM(:,2).*Sh  + KaM(:,1).*KaM(:,2).*KaM(:,3)  +  Sh^3 .* (1 + KdM/1));
    spcM(:,1) = (StV .* Sh^3 .* KdM)                             ./ DNM;
    spcM(:,2) = (StV .* Sh^3)                                    ./ DNM;
    spcM(:,3) = (StV .* Sh^2 .* KaM(:,1))                        ./ DNM;
    spcM(:,4) = (StV .* Sh   .* KaM(:,1) .* KaM(:,2))            ./ DNM;
    spcM(:,5) = (StV         .* KaM(:,1) .* KaM(:,2) .* KaM(:,3))./ DNM;
    % Evaluation of the charge balance for the current Sh value, F(Sh).
    Fplus = Sh + sum(sum(spcM.*chrM));

    % Restoration of the original Sh value
    pH = pH_orig;
    % Calculation of all the concentration for the actual Sh value.
    Sh = 10^(-pH);
    spcM = zeros(size(chrM));
    DNM = ( KaM(:,1).*Sh^2  +   KaM(:,1).*KaM(:,2).*Sh  + KaM(:,1).*KaM(:,2).*KaM(:,3)  +  Sh^3 .* (1+KdM/1));
    spcM(:,1) = (StV .* Sh^3 .* KdM)                             ./ DNM;
    spcM(:,2) = (StV .* Sh^3)                                    ./ DNM;
    spcM(:,3) = (StV .* Sh^2 .* KaM(:,1))                        ./ DNM;
    spcM(:,4) = (StV .* Sh   .* KaM(:,1) .* KaM(:,2))            ./ DNM;
    spcM(:,5) = (StV         .* KaM(:,1) .* KaM(:,2) .* KaM(:,3))./ DNM;
    % Evaluation of the charge balance for the current Sh value, F(Sh).
    F = Sh + sum(sum(spcM.*chrM));

    % Numerical evaluation of the derivative of the charge balance for the current value of Sh, F'(Sh).
    dF_dpH = (Fplus-F)/(eps*pH);
    % Newton Raphson algorithm.
    Fp = dF_dpH;
    % Storing the last value of Sh achieved.
    pH_last = pH;
    % Newton-Raphson algorithm.
    pH = pH_last - F/Fp;
    i=i+1;
end

% In case of non convergence we keep going with different starting point.
if i>=maxIter || (isnan(abs((pH-pH_last)/pH_last))==1),
    % Initial pH point changes alternatively to positive and to negative directions.
    pHini = pHini + (-1)^j * j * 0.25;
    pH_last = 100
    j=j+1;
    if pHini>14, brk = 1;   end
else
% Checking if a valid pH was obtained.
pH = -log10(Sh);
if (Sh>1e-14) && (Sh<1),
    brk = 1;
    % Returning the species matrix.
    Sh = 10^(-pH);
    spcM = zeros(size(chrM));
    DNM = ( KaM(:,1).*Sh^2  +   KaM(:,1).*KaM(:,2).*Sh  + KaM(:,1).*KaM(:,2).*KaM(:,3)  +  Sh^3 .* (1+KdM/1));
    spcM(:,1) = (StV .* Sh^3 .* KdM)                             ./ DNM;
    spcM(:,2) = (StV .* Sh^3)                                    ./ DNM;
    spcM(:,3) = (StV .* Sh^2 .* KaM(:,1))                        ./ DNM;
    spcM(:,4) = (StV .* Sh   .* KaM(:,1) .* KaM(:,2))            ./ DNM;
    spcM(:,5) = (StV         .* KaM(:,1) .* KaM(:,2) .* KaM(:,3))./ DNM;
end
end
end
% Returning the result obtained.
fSh = Sh;
fspcM = spcM(1:size(spcM,1)-1,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%