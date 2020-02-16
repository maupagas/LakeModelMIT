% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae
function [sys,x0,str,ts] = SFunction_Lake(t, x, u, flag, idR, Conn)
global R;

% Depending on the flag value the simulator conducts a different task. 
switch flag,
    case 0,
        % Initialisation function.
        [sys, x0, str, ts] = mdlInitializeSizes(idR, Conn);
    case 1,
        [sys] = mdlDerivatives(t, x, u, idR, Conn);
    case 3,
        sys = mdlOutputs(t, x, u, idR, Conn);
    case {2, 4, 9}; 	% Unused flags.
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S-FUNCTION TASK OF INITIALISATION
function [sys,x0,str,ts] = mdlInitializeSizes (idR, Conn)
global R;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If there is no update in mymodelstructure.xls file, the .mat file is loaded much faster if available.
clc;
FileModelXlsx = dir('my*.xlsx');   %%%% FIX THIS READING THE NAME OF THE EXCEL FILE CALLED myXXX.xlsx
R.FileModelXlsx = FileModelXlsx.name;

fileinfoMat = dir('model_loaded.mat');
    if isempty([fileinfoMat.date])==0 && (datenum(fileinfoMat.date)>=datenum(FileModelXlsx.date)),
fprintf('\n>>>> Reading the Model Structure and Parameters from PREVIOUS RUN...\n');
load('model_loaded.mat');
    else
fprintf('\n>>>> Reading the Model Structure and Parameters from EXCEL FILE...\n');
loadModelXls(idR);
    end
fprintf('\n... Model Loaded, Ready to RUN>>>>>\n');        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags of initialisation.
R(idR).ini = 1;
R(idR).t = 0;

% Dimensions of model. Number of Grid points and layers.
R(idR).Dim.numGrd = size(R(idR).StM,2);
R(idR).Dim.numStG = size(R(idR).StM,1);
R(idR).Dim.numLyr = R(idR).Dim.numGrd - 2;

% Function executed here to generate the names and size of the algebraic variables.
[R(idR).AlgSt, R(idR).AlgStM, R(idR).St, R(idR).StM] = my_algebraics(idR, R);    

% Number of states is grid points times grid states.
R(idR).Dim.numSt    = numel(R(idR).StM);
R(idR).Dim.numAlgSt = numel(R(idR).AlgStM);
R(idR).Dim.numAlgStG = size(R(idR).AlgStM,1);
R(idR).Dim.numRcG = size(R(idR).stoM,2);
R(idR).Dim.numTrG = size(R(idR).trpM,2);
% R(idR).Dim.numDGG  = size(R(idR).stoM,2);        %Number of reactions DG to be calculated for each timestep

%Number of states differentiated by its phase
R(idR).Dim.numStG_S = sum(1*strcmp(R(idR).St.Phase,'L'));     % Number of dissolved states.
R(idR).Dim.numStG_B = sum(1*strcmp(R(idR).St.Phase,'B'));     % Number of active biomass states.
R(idR).Dim.numStG_X = sum(1*strcmp(R(idR).St.Phase,'S'));     % Number of solid states.
R(idR).Dim.numStG_G = sum(1*strcmp(R(idR).St.Phase,'G'));     % Number of gas states.

% Number of outputs of the model. States plus thermodynamic states AlgSt plus reactions and transports.
R(idR).Dim.NOutG = R(idR).Dim.numStG + R(idR).Dim.numAlgStG+  R(idR).Dim.numRcG + ...
                   R(idR).Dim.numTrG;% + R(idR).Dim.numDGG + R(idR).Dim.numRcG + R(idR).Dim.numTrG;
R(idR).Dim.NOut  = R(idR).Dim.numSt  + R(idR).Dim.numAlgSt + R(idR).Dim.numRcG*R(idR).Dim.numGrd + R(idR).Dim.numTrG*R(idR).Dim.numGrd;% + R(idR).Dim.numDGG*R(idR).Dim.numGrd;

% Defining the variable dimensions of the S-Function.
sizes = simsizes;
sizes.NumContStates  = R(idR).Dim.numSt;     % States.
sizes.NumDiscStates  = 0;
sizes.NumInputs      = -1;              % Automatic (-1).
sizes.NumOutputs     = R(idR).Dim.NOut;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;               % At least one needed.
sys = simsizes(sizes);


% Initial states vector. 
x0 = [];    str = [];   ts  = [0 0];
for i=1:(R(idR).Dim.numGrd),    
    x0 = [x0; R(idR).StM(:,i)];     
end

% Update of the state variables by name.
for idSt=1:(R(idR).Dim.numStG),
    R(idR).St.(char(R(idR).St.StNames(idSt))) = R(idR).StM(idSt,:);
end
% END OF INITIALISATION TASK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S-FUNCTION TASK OF DERIVATIVES CALCULATION
function [sys] = mdlDerivatives(t, x, u, idR, Conn)
global R;

% Setting initialisation flag to 0 and enabling time to be read for inflow function.
R(idR).ini = 0;        R(idR).t = t;  R(idR).upd = 1;

% Prints time in screen only every 24 units.
if mod(t,24)==0 && t>0,
    clc;
    t
%     pause(1)
%     plotLake3D(tout, allOut,{'pH'});      %I DO NOT KNOW HOW TO USE THIS
%     YET (Mauricio 20/3/18)
end

% %%%%%%%%%% CUSTOM AREA %%%%%%%%%%%%%%
% % If connected for inputs we read them here.
% switch Conn,
% case 1,
% case 2,
% end
% %%%%%%% END OF CUSTOM AREA %%%%%%%%%%%

% Zero or negative values are corrected to avoid negative concentrations and log(0) events.
x = (x>=0).*x + (x<0)*eps;

% Update of the state variables from x system state variables.
%First term keeps the calculated value in x for the variables labeled as
%dynamic and cancels its value for those labeled as not dynamic. Second
%term of the sum adds the value calculated for the non-dynamic variables to
%update the States Matrix. 

for idGrd=1:R(idR).Dim.numGrd,          
    R(idR).StM(:,idGrd) = x(((idGrd-1)*R(idR).Dim.numStG+1)   :  idGrd*R(idR).Dim.numStG) .* R(idR).St.Dyn ...
                          + R(idR).StM(:,idGrd) .* (R(idR).St.Dyn == 0); 
end
                      
%Passing up the values from States matrix to States Names.                    
for idSt =1:R(idR).Dim.numStG,      
    R(idR).St.(char(R(idR).St.StNames(idSt))) = R(idR).StM(idSt,:);     
end

%% Update of the feeding stream.  (PREVIOUS VERSION OF LOADING FEED)
% % Depending on whether it is connected to other unit or from file.
% [Fd, FdV] = my_feeding(t, u, idR, Conn);
% Xin = FdV;              % kmol/m3


%% Loading functions to be used for the derivatives
%%Update of the feeding stream.
% Depending on whether it is connected to other unit or from file.
% [R(idR).Fd] = my_inflow(t, R, idR, Conn);

% All algebraic calculations necessary for kinetics or stoichiometry, are carried out.
[R(idR).AlgSt, R(idR).AlgStM, R(idR).St, R(idR).StM] = my_algebraics(idR, R);              

%Reaction and transport rates.  
[R(idR).rRc, R(idR).rTr, R(idR).rRcM, R(idR).rTrM] = my_kinetics(idR, R);

% Net generation rates by application of the stoichiometry matrix to the
% vector of reaction rates previously calculated.
[R(idR).Rt, R(idR).RtM] = my_stoichiometry(idR, R);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the reactor mass balances.
[dStM_dt] = my_reactor(idR, R);

% Derivatives assignation for all the state variables.
sys = zeros(R(idR).Dim.numSt,1);
sys(1:R(idR).Dim.numSt) = dStM_dt;
for idGrd =1:R(idR).Dim.numGrd,
    sys(((idGrd-1)*R(idR).Dim.numStG+1):idGrd*R(idR).Dim.numStG) = dStM_dt(:,idGrd);
end
% END OF DERIVATIVES CALCULATION TASK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONLY THE FIRST TIME AND ONCE TO WRITE THE NAMES OF THE OUTPUT VARIABLES
% IN EXCEL (FOR WHAT IS THIS? MAURICIO (20/3/18)
% Definition of the vector of names vector of the output variables and rates.
% clear hd;
% if R(idR).ini==1,
%     hd(1:((R(idR).Dim.NOut/R(idR).Dim.numGrd)+1)) = {''};
%     % After reading the headers are writen in the excel file for output later.
%     hd(1) = {'Time'};
% 
%     % Names of the state variables.
%     for i=2:R(idR).Dim.numStG+1,    hd(i) = {char(R(idR).St.StNames(i-1))};    end
% 
%     % Number of the algebraic state variables.
%     for j=(i+1):(i+R(idR).Dim.numAlgStG), hd(j) = {char(R(idR).AlgSt.AlgStNames(j-i))};    end
% 
% %     % Names of the reaction rates.
% %     for i=(j+1):(j+R(idR).Dim.numRcG),  hd(i) = {char(R(idR).rRc.RcNames(i-1))};    end
% % 
% %     % Names of the transport rates.
% %     for i=(j+1):(j+R(idR).Dim.numTrG),  hd(i) = {char(R(idR).rTr.RcNames(i-1))};    end
% % 
%     % Write in the Excel output file the variables names.
%     for idGrd=1:R(idR).Dim.numGrd,
% %         xlswrite('simResults.xls', hd, strcat('OutputL', num2str(idGrd)));
% %         xlswrite('simResults.xls', {strcat('L', num2str(idGrd))}, strcat('OutputBF', num2str(idGrd)), 'M1:M1');
%     end
%     R(idR).ini = 0;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF DERIVATIVES CALCULATION TASK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S-FUNCTION TASK OF OUTPUTS ASSIGNATION
function sys = mdlOutputs(t, x, u, idR, Conn)
global R;

% Update of the state variables
% for idGrd=1:R(idR).Dim.numGrd,
%     R(idR).StM(:,idGrd) = x(((idGrd-1)*R(idR).Dim.numStG+1)   :  idGrd*R(idR).Dim.numStG);
% end
% for idSt=1:(R(idR).Dim.numStG),
%     R(idR).St.(char(R(idR).St.StNames(idSt))) = R(idR).StM(idSt,:);
% end

% If it is the first time at t=0, the output is only the state and the rest of variables are set to zero.
if t==0,
    sys = zeros(R(idR).Dim.NOut,1);
    sys(1:R(idR).Dim.numSt) = x(1:R(idR).Dim.numSt);
% The output contains already the whole information available.
else

%Prints the output for the STATE VARIABLES and for the ALGEBRAIC STATES. If
%other outputs are desired to be printed, add at the endd of the matrix
%developed inside the for loop. Here rates can be added, conc. of
%individual species and other things as well.
    ROut = zeros(R(idR).Dim.NOut,1);   %Prepares the output for every layer in a matrix
        for idGrd=1:R(idR).Dim.numGrd,
            ROut(((idGrd-1)*R(idR).Dim.NOutG+1):(idGrd*R(idR).Dim.NOutG)) = [R(idR).StM(:,idGrd); R(idR).AlgStM(:,idGrd); R(idR).rRcM(:,idGrd); R(idR).rTrM(:,idGrd)];% R(idR).AlgSt.DGV(:,idGrd)
           % Printing out rates along with Gibbs calculations
            if t>0,
                ROut(((idGrd-1)*R(idR).Dim.NOutG+1):(idGrd*R(idR).Dim.NOutG)) = [R(idR).StM(:,idGrd); R(idR).AlgStM(:,idGrd); R(idR).rRcM(:,idGrd); R(idR).rTrM(:,idGrd)];% R(idR).AlgSt.DGV(:,idGrd)
            % R(idR).Dim.NOut  = R(idR).Dim.numSt  + R(idR).Dim.numAlgSt +  R(idR).Dim.numRcG*R(idR).Dim.numGrd + R(idR).Dim.numTrG*R(idR).Dim.numGrd + R(idR).Dim.numDG*R(idR).Dim.numGrd;
            end
        end     
        
    sys = ROut;
end
% end
% END OF OUTPUTS ASSIGNATION TASK
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**********************************************************************************************
% Update of the state variables
% for idGrd=1:R(idR).Dim.numGrd,
%     R(idR).StM(:,idGrd) = x(((idGrd-1)*R(idR).Dim.numStG+1)   :  idGrd*R(idR).Dim.numStG);
% end
% for idSt=1:(R(idR).Dim.numStG),
%     R(idR).St.(char(R(idR).St.StNames(idSt))) = R(idR).StM(idSt,:);
% end

% If it is the first time at t=0, the output is only the state and the rest of variables are set to zero.
% if t==0,
%     sys = zeros(R(idR).Dim.NOut,1);
%     sys(1:R(idR).Dim.numSt) = x(1:R(idR).Dim.numSt);
% % The output contains already the whole information available.
% else
%**************************************************************************************************
