% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae

%Last version released: 10/10/16

%%% This script calculates the Gibbs energy for each one of the reactions
%%% included in the model. For energy calculations, activiies instead of
%%% concentrations are used as required.
function [DGV] = my_thermo(idR, R, numGrd)

% Reading all parameters required for this function
St    = R(idR).St;          StM    = R(idR).StM;       
AlgSt = R(idR).AlgSt;       pTh    = R(idR).pTh;           
pOp   = R(idR).pOp;         stoMg  = R(idR).stoMg; 
pPhys = R(idR).pPhys;

%Defining pH, species activity Matrix and Temperature for calculations later on
pH   = AlgSt.pH(numGrd);
spaM = AlgSt.spaM(:,:,numGrd);
T = pPhys.T(numGrd);

% Construction of the Gibbs (corrected with temperature) matrix for each species.
% A line of zeros is added in the last line in order to account for the
% Gibbs activity of the proton, which is zero.
G0ftM = [pTh.G0ftM(:,:,numGrd); zeros(1, length(pTh.G0ftM(1,:,1)))] ;

% Construction of the species activity matrix for all of the components
% considered: The number of the states is the size of the States Matrix plus
% 2 (corresponds to water and proton activities). For the particulate
% components, their activities for Gibbs calculations is 1 as they are in a solid non-reactive state.
% Gases components will be included with their own concentrations.
% The species activity matrix need to be deconstructed to reposition proton
% and water at the very bottom of the calculation.

% Builds a matrix of zeros with same size as G0tfM
spaM_DG = zeros(size(G0ftM));

% Selecting the activities for the dissolved components except water.
% Identification of liquid states in all compartments.
idL = (strcmp(St.Phase,'L'));    
spaM_DG(idL,:) = spaM(idL,:);

% For the non-dissolved components, the activity of them for THERMODYNAMICS CALCULATIONS
% is assumed as 1. They are positioned in column 2 (fully protonated and hydrated)
spaM_DG(idL==0,2) = 1;

%Identification of gaseous states in all compartments
idG = (strcmp(St.Phase,'G'));      

% Overwrites the 1s assigned to non-liquid components by the pressure of the gas considered.
% Pressure of gases are included to calculate energetics of reactions where gaseous phase
% components are involved (i.e. Henry Constant calculation)
spaM_DG(idG,2) = StM(idG);

% Appends the activity of the water and the OH-
spaM_DG(length(spaM_DG(:,1))-1,:) = spaM(length(spaM(:,1)),:);
% Appends the activity of the proton in the second column
spaM_DG(length(spaM_DG(:,1)),2) = 10^-pH;

% Creates a vector with the Gibbs formation values (corrected by temperature) for each one of 
% the species defined in ThermoFlag (Number 2 indicates that it sums the columns of the matrix)
% This vector is used later on to calculate the chemical potentials of each component. 
DG0V  = sum(G0ftM .* pTh.spcSelM, 2);

% Creates a vector for the activities of the species used in the metabolic reactions. 
% The species are those defined in myModel Excel file by the ThermoFlag
% label in ReactMatrix tab.
% (Number 2 indicates that it sums the columns of each file of the matrix)
spaV_DG = sum(spaM_DG  .* pTh.spcSelM, 2);

% Checking if some component concentration is equal to zero to avoid ln(0)= -Inf
% In those cases 0 is replaced by eps :
for i = 1:length(spaV_DG),
    if spaV_DG(i) == 0, spaV_DG(i) = eps; end
end
%%% BETTER ALTERNATIVE??? spaV_DG = (spaV_DG==0)*eps + spaV_DG;


% Calculates chemical potentials mu for each species involved in the reaction
my_a = DG0V + pOp.Rth * T * (log(spaV_DG));

% Multiplies every column of stoichiometric matrix by the Gibbs of every
% individual component and then sums all the elements to calculate the
% overall value of Gibbs reaction
DGV = sum(bsxfun(@times, my_a, stoMg));

DGV = DGV';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    