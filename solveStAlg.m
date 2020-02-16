function [F] = solveStAlg(StAlg, R, idR)

%Flag for not writing names inside the function
R(idR).upd = 0;
if StAlg < 0, StAlg = StAlg^2; end

% Calculation of the dilution rate and influent vector for the mass balance equation.
D = R(idR).Fd.Qinf / R(idR).pOp.Vliq; 
InfV = R(idR).Fd.FdV(3:length(R(idR).Fd.FdV));
    
% Finding variable(s) position to be calculated algebraically
idsp = find(R(idR).St.Alg == 1);

R(idR).StM(idsp) = StAlg;
R(idR).St.(char(R(idR).St.StNames(idsp))) = StAlg;
R(idR).AlgSt.spaM(idsp,2) = StAlg;

% First calculation of thermo, kinetics and stoichiometry rates to give the
% initial point
    %Recalculate Gibbs energy values
    [R(idR).AlgSt.DGV] = my_thermo(idR, R);
    
    %Recalculate kinetics
    [R(idR).rRc, R(idR).rTr, R(idR).rRcM, R(idR).rTrM] = my_kinetics(idR, R);

    %Recalculate rates
    [R(idR).Rt, R(idR).RtM] = my_stoichiometry(idR, R);

%Loop to calculate algebraically the state via Newton-Raphson method. Error is

%Calculate species mass balance in function of the obtained rates
    F = D * (InfV(idsp) - R(idR).StM(idsp)) + R(idR).RtM(idsp);
