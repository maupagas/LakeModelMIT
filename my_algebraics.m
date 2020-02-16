% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae

% my_algebraics.m - CALCULATES THE ALGEBRAIC VARIABLES FROM THE STATE VARIABLES.
% This function is almost fully customizable by the user on an easy manner 
% since the names from the structre variables can be used to write the equations 
% and the algebraic variables are created inside a structre called AlgSt.
function [AlgSt, AlgStM, St, StM] = my_algebraics(idR, R) 
% %%%%%%%%%%%%%%%%%% COMMON CODE SECTION, GENERALLY DO NOT EDIT THIS SECTION %%%%%%%%%%%%%%% %

%Shortening the names for the variables
St     = R(idR).St;                    StM   = R(idR).StM;
AlgSt  = R(idR).AlgSt;                 stoMg = R(idR).stoMg; 
pTh    = R(idR).pTh;                   pOp   = R(idR).pOp; 
numGrd = R(idR).Dim.numGrd;            pPhys  = R(idR).pPhys;
  
% Equilibrium constants are calculated for each layer (if applicable) for different temperatures and put in a 3D matrix.
% KdM  = pTh.KdaM(:,1);
% KaM  = pTh.KdaM(:,2:length(pTh.KdaM(1,:)));
% chrM = pTh.chrM;
%         
% Identification of liquid states in all compartments.
idL = (strcmp(St.Phase,'L'));     % Row indexes of dissolved states.
%Selection only the state variables that are on liquid phase
StM_l = StM(idL==1,:,:);
% Adding water concentration as 1 to the States Matrix to compute the pH.
StM_lw = [StM_l; ones(1,numGrd)];

% Passing the KaM, KdM and chrM for each layer 
%(POSSIBILITY OF RECALCULATING Kas AS FUNCTION OF TEMPERATURE HERE FROM DATA FROM pPhys)
% for i=1:numGrd,
%     KaMw(:,:,i)  = KaM;
%     KdMw(:,:,i)  = KdM;
%     chrMw(:,:,i) = chrM;     
% end 

% Definition of the IonSolver Mode: (0) Solve all; (1) Contant Ic, solve the rest; (2) pH given solve the rest.
IonSolvMode = 0;    % IF MODE 1, AS IT SI UNCOMMON JUST MODIFY THIS LINE OF DEVELOP SPECIFIC CODE FOR HANDLING.

%Time identifier to take the feeding from the flow matrix. Identifies the position of the first
%value that is greater than t in the time matrix. 
% if Prf.flag == 1,
%     t_id = find(R(idR).t >= Prf.Time,1,'last');
%     AlgSt.pH = Prf.pH(t_id);
%     IonSolvMode = 2;
% end


% Calling the ionic speciation solver for each layer (if applicable).
% Last value is a flag for pH, 0 being calculated pH and Ic, 1, fixed
% Ic and 2 fixed pH. (look at solveIonSpec for more details). 

%%% NEED TO INDEX THE OUTPUTS TO CALCULATE FOR EVERY LAYER AND PREDEFINE
%%% ITS SIZE %%%

Ic = zeros(1,numGrd); pH = Ic; Sh = Ic;
spaM = zeros(length(StM_lw(:,1)),length(pTh.KaM(1,:,1))+1,numGrd); spcM = spaM;

for i=1:numGrd,
    [Ic(i), pH(i), Sh(i), spaM(:,:,i), spcM(:,:,i)] = solveIonSpec(AlgSt.Ic(i), AlgSt.pH(i), StM_lw(:,i), pTh.KdM(:,:,i), pTh.KaM(:,:,i), pTh.chrM, pPhys.T(i), IonSolvMode);
    AlgSt.Ic(i) = Ic(i);
    AlgSt.pH(i) = pH(i);
    AlgSt.Sh(i) = Sh(i);
    AlgSt.spaM(:,:,i) = spaM(:,:,i);
    AlgSt.spcM(:,:,i) = spcM(:,:,i);
end

% %%%%%%%%%%%%%%%%%% END OF COMMON CODE, GENERALLY DO NOT EDIT ABOVE THIS LINE %%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% %%%%%%%%%%% CUSTOMISABLE AREA TO DEFINE THE ALGEBRAIC VARIABLES %%%%%

    
    % %%%%%%%% IONIC SPECIES AS FUNCTION OF pH %%%%%%%%
    % Naming algebraic states in the anode from the species matrix.
    
    %Permute matrix to obtain every algebraic state as a vector
    spcM_perm = permute(spcM,[1,3,2]);

                                     AlgSt.Sach   = spcM_perm(2,:,2);       AlgSt.Sac_   = spcM_perm(2,:,3);
    AlgSt.Sco2 = spcM_perm(4,:,1);   AlgSt.Sh2co3 = spcM_perm(4,:,2);       AlgSt.Shco3_ = spcM_perm(4,:,3);        AlgSt.Sco3_ = spcM_perm(4,:,4);
                                     AlgSt.Snh4   = spcM_perm(6,:,2);       AlgSt.Snh3   = spcM_perm(6,:,3);
                                     AlgSt.Shno2  = spcM_perm(7,:,2);       AlgSt.Sno2_  = spcM_perm(7,:,3);
                                     AlgSt.Ssh2   = spcM_perm(10,:,2);      AlgSt.Ssh_   = spcM_perm(10,:,3);
                                                                            AlgSt.Soh    = spcM_perm(17,:,3);
                                % Returning the speciation vlues to the global variable.
%     AlgSt.spcM = spcM;        
        
    % %%%%%%%%%%% GAS VOLUME PER LAYER FUNCTION OF TOTAL GAS MOLES AND T AND P %%%%%%%%%%%%%%%%%%%
    idG = find(strcmp(St.Phase,'G'));
    GgasT = zeros(1,numGrd);
    for i=1:length(idG),
        GgasT = GgasT + StM(idG(i),:);  % molG/Lliq
    end
    AlgSt.Vgas = GgasT * pOp.Rg .* pPhys.T ./ pPhys.P;   % Lg/Lliq
    
    
    % %%%%%%% END OF CUSTOMISABLE AREA %%%%%%%%%%%
R(idR).AlgSt.spcM = spcM;
R(idR).AlgSt.spaM = spaM;
R(idR).AlgSt.pH   = pH;    
% %%%%%%%SOLVER ALGEBRAICALLY ONLY WITH FZERO (DOUBLE CHECK FOR LAKE MODEL
% // IGNORED FOR NOW #MPG: 258/3/18)
if any(R(idR).St.Alg == 1) && R(idR).ini == 0,
  
    idsp = find(R(idR).St.Alg == 1);
    StAlg_ini = StM(idsp);
    
    %Define the options for the fzero solver (show graphics, results of evaluations or not)
    %options = optimset('Display', 'iter', 'FunValCheck', 'on', 'PlotFcns', @optimplotfval, 'TolX', eps);
    options = 0;
    
    %Passed into the fzero function the variable R as local and the index of reactor
    % [StAlgCalc, fval, exitflag, output] = fzero(@solveStAlg_fzero, [0, 1], options, R, idR);
    [StAlgCalc] = fzero(@solveStAlg, StAlg_ini, options, R, idR);
    
    %Update into the species activity matrix (in both global and local variable) the calculation of the
    %algebraic component (Watch out if component is ionized!)
    R(idR).AlgSt.spaM(idsp,2) = StAlgCalc;
    R(idR).AlgSt.spcM(idsp,2) = StAlgCalc;
    AlgSt.spaM(idsp,2) = StAlgCalc;
    AlgSt.spcM(idsp,2) = StAlgCalc;

    %Update also into the States values by name and inside the matrix
    St.(char(R(idR).St.StNames(idsp))) = StAlgCalc;
    StM(idsp) = StAlgCalc;
end

% Calls the Gibbs calculation function (applied to each layer
for i=1:numGrd,
    [AlgSt.DGM(:,i)] = my_thermo(idR, R, i);    
end    

% Passing the values of the variables to the Algebraic States Matrix (NEED TO MAKE THE MATRIX BIGGER????) 
AlgStNames_length = length(AlgSt.AlgStNames);
AlgStM = zeros(AlgStNames_length, numGrd);
DG_id = AlgStNames_length - length(AlgSt.DGM(:,1));  %Returns the numer of Algebraic states that are NOT Gibbs calculations
% DGM_perm = permute(AlgSt.DGM,[1,3,2]);

for i=1:AlgStNames_length,      
    if i<=DG_id,
        AlgStM(i,:) = AlgSt.(char(AlgSt.AlgStNames(i,:)));  
    else
        AlgStM(i,:) = AlgSt.DGM(i-DG_id,:);  
    end    
end


%%%%%%%%%%%ALTERNATIVE
% for i=1:length(AlgSt.AlgStNames),      
%     AlgStMv2(i,:) = AlgSt.(char(AlgSt.AlgStNames(i,1)))
% end

% Output to the global R.
R(idR).AlgSt = AlgSt;
R(idR).AlgStM = AlgStM;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%