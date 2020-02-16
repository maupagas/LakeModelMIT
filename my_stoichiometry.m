% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae

function [Rt, RtM] = my_stoichiometry(idR, R)
% Reading all parameters and variables for this reactor.
stoM = R(idR).stoM;             trpM  = R(idR).trpM;
rRcM = R(idR).rRcM;             rTrM  = R(idR).rTrM;
St   = R(idR).St;               AlgSt = R(idR).AlgSt;
numGrd = R(idR).Dim.numGrd;     numStG = R(idR).Dim.numStG;

%% Concatenate Stoichiometric Matrix (no water nor H+ included) with TranspMatrix
vStoM = [stoM trpM];
rM    = [rRcM; rTrM];
% Calculation of the net generation layer by layer.
RtM = vStoM * rM;

% Structure with the net formation of each state species.
if R(idR).upd ~= 0,
    for i=1:(numStG),   Rt.(char(St.StNames(i))) = RtM(i,:);    end
else
    Rt = [];
end