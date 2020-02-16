% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae

function [] = loadModelXls(idR)
global R;
warning off;

% If only one reactor is used is set to '' and no extra for the names of the xls sheets is needed.
if idR==1,  id = '';    else    id = num2str(idR);   end
FileModelXlsx = R.FileModelXlsx;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL MODEL PARAMETERS from the Excel file.
fprintf('\n> Loading and creating GENERAL MODEL PARAMETERS...');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL STATES structure and initial values.
[States,      StNamesP]   = xlsread(FileModelXlsx, strcat('States', id), '', 'basic');
aux = '';
StNames = StNamesP(:,1);
for i=1:size(States,1),
    % Building the string command for the structure, last without ', '.
    if i<size(States,1),   
        aux = strcat(aux, char(39), StNames(i), char(39), ', 0, ');
    else
        aux = strcat(aux, char(39), StNames(i), char(39), ', 0');
    end 
end

aux = char(aux);
%Creates an structure with every variable name with its own value.
eval(strcat('St = struct(', aux, ');'));
St.StNames = StNames;

%Identify by its label if species are to be calculated dinamically (d) or not.
%If it is calculated algebraically, it will be read under St.Alg
St.Dyn = strcmp(StNamesP(:,size(StNamesP,2)),'d');
St.Alg = strcmp(StNamesP(:,size(StNamesP,2)),'a');
%Reading the phase of the species (L: Liquid, S: Solid, X: Biomass, G: Gaseous
St.Phase = StNamesP(:,size(StNamesP,2)-1);
StM = States;

% Giving initial values vector (element for each layer) of each state variable. 
for i=1:length(StNames),
    eval(strcat('St.', char(StNames(i)), ' = [', char(num2str(StM(i,:))), '];'));
end
clear States StNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALGEBRAIC STATES structure and initial values.
[AlgStates,   AlgStNames]= xlsread(FileModelXlsx, strcat('AlgStates',  id), '', 'basic');
aux = '';
for i=1:size(AlgStates,1),
    % Building the string command for the structure, last without ', '.
    if i<size(AlgStates,1),   
        aux = strcat(aux, char(39), AlgStNames(i,1), char(39), ', 0, ');
    else
        aux = strcat(aux, char(39), AlgStNames(i,1), char(39), ', 0');
    end
end
aux = char(aux);
eval(strcat('AlgSt = struct(', aux, ');'));
AlgSt.AlgStNames = AlgStNames(:,1);
AlgStNames = AlgSt.AlgStNames;
AlgStM = zeros(size(AlgStates));
% Giving initial values vector (element for each layer) of each state variable. 
for i=1:length(AlgStNames),
    eval(strcat('AlgSt.', char(AlgStNames(i)), ' = [', char(num2str(AlgStates(i,:))), '];'));
end
clear AlgStates AlgStNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % REACTION MATRIX structure (OLD VERSION LAKE MODEL)
% [ReactMatrixDG,text] = xlsread(FileModelXlsx, strcat('ReactMatrix', id));
% ReactMatrix = ReactMatrixDG(1:size(ReactMatrixDG,1)-1,:);
% 
% % targeting the right excel cells for the names.
% if size(ReactMatrix,2)<=25,
%     rng = strcat('B1:', char(65+size(ReactMatrix,2)), '1');
%     [N, RcNames] = xlsread(FileModelXlsx, strcat('ReactMatrix', id), rng);
% else
%     a = char(floor((size(ReactMatrix,2)+1)/26) +64);
%     b = char(mod  ((size(ReactMatrix,2)+1),26) +64);
%     [N, RcNames] = xlsread(FileModelXlsx, strcat('ReactMatrix', id), strcat('B1:', a, b ,'1'));
% end
% RcNames = RcNames';
% stoM = ReactMatrix;
% rRc.RcNames = RcNames;
% clear N ReactMatrix RcNames;

%% REACTION MATRIX structure.
[ReactMatrix ,ReactText] = xlsread(FileModelXlsx, strcat('ReactMatrix', id), '', 'basic');

%Selects the names for the transport rates. Ignores the last column as it
%is assigned to ThermoFlag
RcNames = ReactText(1,2:length(ReactText(1,:)) -1);
RcNames = RcNames';

% The stoichiometry matrix used for the mass balance of the model is
% read. It does not take into account the water and proton stoichiometry
stoM = ReactMatrix(1:(length(ReactMatrix(:,1))-2), 1:(length(ReactMatrix(1,:))-1));

%Stoichiometry matrix with water and proton included. Will be
%used when loading Thermodynamic Parameters and  for Gibbs 
% reaction calculations in my_algebraics
stoMg = ReactMatrix(:,1:(length(ReactMatrix(1,:))-1));

% Flag value for Gibbs energy columns to use. This is read from the stoichiometry matrix. It is used to identify the
% components involved in the Gibbs reaction calculation. Will be used when 
% loading Thermodynamic parameters to construct a matrix for species
% selection. 2 is added to provide a valid index for the matrix that will
% be created.
spcThFlag = ReactMatrix(:,(length(ReactMatrix(1,:)))) +2;     

rRc.RcNames = RcNames;
clear N ReactMatrix RcNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% TRANSPORT MATRIX structure. (OLD VERSION LAKE MODEL)
% [TranspMatrix,TrNames] = xlsread(FileModelXlsx, strcat('TranspMatrix', id));
% if size(TranspMatrix,2)<=25,
%     rng = strcat('B1:', char(65+size(TranspMatrix,2)), '1');
%     [N, TrNames] = xlsread(FileModelXlsx, strcat('TranspMatrix', id), rng);
% else
%     a = char(floor((size(TranspMatrix,2)+1)/26) +64);
%     b = char(mod  ((size(TranspMatrix,2)+1),26) +64);
%     [N, TrNames] = xlsread(FileModelXlsx, strcat('TranspMatrix', id), strcat('B1:', a, b ,'1'));
% end
% TrNames = TrNames';
% trpM = TranspMatrix;
% rTr.TrNames = TrNames;
% clear N TranspMatrix TrNames;

%% TRANSPORT MATRIX structure.
[TranspMatrix, TrText] = xlsread(FileModelXlsx, strcat('TranspMatrix', id), '', 'basic');

%Selects the names for the transport rates. Ignores the last column as it
%is assigned to ThermoFlag
TrNames = TrText(1,2:length(TrText(1,:)) -1);
TrNames = TrNames';

%Transport Matrix for the mass balance of the model. Analogous to stoM,
%doesn't take into account H2O and H+ transport
trpM =  TranspMatrix(1:(length(TranspMatrix(:,1))-2),1:(length(TranspMatrix(1,:))-1));
%Transport Matrix taking into consideration water and proton, mainly for gas transfer or precipitation
%processes where water or proton can be involved
trpMg = TranspMatrix(:,1:(length(TranspMatrix(1,:))-1));

rTr.TrNames = TrNames;
clear N TranspMatrix TrNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRAULIC FLOW TRANSPORT concentrations structure.
% [Flow,      FlwNamesP]   = xlsread(FileModelXlsx, strcat('FlowProgram',       id));
% aux = '';
% FlwNames = strcat(FlwNamesP(:,1),'_inf');
% for i=1:size(Flow,1),
%     % Building the string command for the structure, last without ', '.
%     if i<size(Flow,1),   
%         aux = strcat(aux, char(39), FlwNames(i), char(39), ', 0, ');
%     else
%         aux = strcat(aux, char(39), FlwNames(i), char(39), ', 0');
%     end 
% end
% aux = char(aux);
% eval(strcat('Flw = struct(', aux, ');'));
% Flw.FlwNames = FlwNames;
% FlwM = Flow;
% Giving initial values vector (element for each layer) of each state variable. 
% for i=1:length(FlwNames),
%     eval(strcat('Flw.', char(FlwNames(i)), ' = [', char(num2str(FlwM(i,:))), '];'));
% end
% FlwM = Flow(2:size(Flow,1),:);
% clear Flow FlwNames;

%% HYDRAULIC FLOW TRANSPORT concentrations structure. 
% CONSTANT INFLOW INPUT
% Modified from Lake model in order to obtain Matrix of Time with influent concentrations.
[Flow,      FlwNamesP]   = xlsread(FileModelXlsx, strcat('Inflow', id), '', 'basic');
aux = '';
FlwNames = strcat(FlwNamesP(:,1),'_inf');
for i=1:size(Flow,1),
    % Building the string command for the structure, last without ', '.
    if i<size(Flow,1),   
        aux = strcat(aux, char(39), FlwNames(i), char(39), ', 0, ');
    else
        aux = strcat(aux, char(39), FlwNames(i), char(39), ', 0');
    end 
end
aux = char(aux);
eval(strcat('Flw = struct(', aux, ');'));
Flw.FlwNames = FlwNames;
FlwM = Flow;
% Giving initial values vector (element for each layer) of each state variable. 
for i=1:length(FlwNames),
    eval(strcat('Flw.', char(FlwNames(i)), ' = [', char(num2str(FlwM(i,:))), '];'));
end
FlwM = Flow(2:size(Flow,1),:);
clear Flow FlwNames;


%ALTERNATIVE FOR MATRIX WITH INPUT FLOW
% [Inflow, FlwNames] = xlsread(FileModelXlsx, strcat('Inflow', id), '', 'basic');
% aux = '';
% FlwNames = FlwNames(2:size(FlwNames,2));
% for i=1:length(FlwNames),
%     if i<length(FlwNames), 
%         aux = strcat(aux, char(39), FlwNames(i), char(39), ', Inflow(:,', num2str(i+1), '), '); 
%     else
%         aux = strcat(aux, char(39), FlwNames(i), char(39), ', Inflow(:,', num2str(i+1), ')'); 
%     end   
% end
% aux = char(aux);
% eval(strcat('Flw = struct(', aux, ');'));
% Flw.FlwM = Inflow(:,2:size(Inflow,2)); 
% Flw.Time = Inflow(:,1);
% Flw.FlwNames = FlwNames;
% clear Flow FlwNames;


% Output to the global variable R.
R(idR).St    = St;      clear St;
R(idR).StM   = StM;     clear StM;
R(idR).AlgSt = AlgSt;   clear AlgSt;
R(idR).AlgStM= AlgStM;  clear AlgStV;
R(idR).trpM  = trpM;    clear trpM;
R(idR).trpMg = trpMg;   clear trpMg;
R(idR).stoM  = stoM;    clear stoM;
R(idR).stoMg = stoMg;   clear stoMg;
R(idR).rRc   = rRc;     clear rRc;
R(idR).rTr   = rTr;     clear rTr;
R(idR).Flw   = Flw;     clear Flw;
R(idR).FlwM  = FlwM;    clear FlwM;
R(idR).ini = 1;         % Indicates first evaluation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Structure with the OPERATIONAL PARAMETERS.
[OperatParam, pOpNames]   = xlsread(FileModelXlsx, strcat('OperatParam', id), '', 'basic');
aux = '';
for i=1:length(OperatParam),
    % Building the string command for the structure, last without ', '.
    if i<length(OperatParam),   
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)),', '); 
    else
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i))); 
    end
end
aux = char(aux);
eval(strcat('pOp = struct(', aux, ');'));
pOp.pOpNames = pOpNames(:,1);
clear OperatParam pOpNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Structure with the PHYSICAL PARAMETERS. (CHECK FOR BASIC PROCEDURE)
[PhysParam,    pPhysNames]  = xlsread(FileModelXlsx, strcat('PhysParam', id), '', 'basic');
aux = '';
for i=1:size(PhysParam,1),
    % Building the string command for the structure, last without ', '.
    if i<size(PhysParam,1),   
        aux = strcat(aux, char(39), pPhysNames(i), char(39), ', 0, '); 
    else
        aux = strcat(aux, char(39), pPhysNames(i), char(39), ', 0'); 
    end
end
aux = char(aux);
eval(strcat('pPhys = struct(', aux, ');'));
pPhys.pPhysNames = pPhysNames(:,1);
for i =1:size(PhysParam,1),
    pPhys.(char(pPhysNames(i))) = PhysParam(i,:);
end
clear PhysParam pPhysNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% SPECIAL SCRIPT FOR DENSITY CALCULATION AS f(T).
% Density of water (pW) calculated here for each T (ºC).
numGrd = size(R(idR).StM,2);
Tc_pW = [100000 958.4; 100 958.4; 80 971.8; 60 983.2; 40 992.2; 30 995.6502; 25 997.0479; 22 997.7735; 20 998.2071; 15 999.1026; 10 999.7026; 4 999.972; 0 999.8395; -100000 983.854];
X = Tc_pW(:,1);    Y = Tc_pW(:,2);    Tc = pPhys.T-273.15;
pPhys.pW = ones(1,numGrd)*1000;
for i=1:numGrd, 
    pPhys.pW(i) = interp1(X,Y,Tc(i)); 
end
%% %%% END OF SCRIPT FOR DENSITY CALCULATION AS f(T)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure with the THERMODYNAMIC PARAMETERS. (OLD LAKE MODEL VERSION)
% [ThParam,    pThNames]  = xlsread(FileModelXlsx, strcat('ThermoParam',    id));
% aux = '';
% numGrd = size(R(idR).StM,2);
% for i=1:size(ThParam,1),
%     % Building the string command for the structure, last without ', '.
%     if i<size(ThParam,1),   
%         aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:), '); 
%     else
%         aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:)'); 
%     end
% end
% aux = char(aux);
% eval(strcat('pTh = struct(', aux, ');'));
% pTh.pThNames = pThNames(:,1);
% pTh.Type = pThNames(:,size(pThNames,2));
% idnSp = prod(1*(char(pTh.Type)==repmat('nSp',size(pThNames,1),1)),2)==1;
% idG0f = sum([(prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('G0l',size(ThParam,1),1)),2)==1) ...
%     (prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('G0s',size(ThParam,1),1)),2)==1) ...
%     (prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('G0g',size(ThParam,1),1)),2)==1)],2);
% idH0f = sum([(prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('H0l',size(ThParam,1),1)),2)==1) ...
%     (prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('H0s',size(ThParam,1),1)),2)==1) ...
%     (prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('H0g',size(ThParam,1),1)),2)==1)],2);
% idchr = prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('chr',size(ThParam,1),1)),2)==1;
% pTh.spcNames= [];
% pTh.G0fM    = zeros(sum(idG0f),size(ThParam,2));   jG = 0;
% pTh.H0fM    = zeros(sum(idH0f),size(ThParam,2));   jH = 0;
% pTh.chrM    = zeros(sum(idchr),size(ThParam,2));   jch = 0;
% for i=1:size(pThNames,1),
%     if     idnSp(i)==1,
%         pTh.spcNames = [pTh.spcNames; pThNames(i,2:size(pThNames,2)-1)];
%     end
% end
% for i=1:size(ThParam,1),
%     if     idG0f(i)==1,
%         jG = jG+1;
%         pTh.G0fM(jG,:) = ThParam(i,:);
%     elseif idH0f(i)==1,
%         jH = jH+1;
%         pTh.H0fM(jH,:) = ThParam(i,:);
%     elseif idchr(i)==1,
%         jch = jch+1;
%         pTh.chrM(jch,:) = ThParam(i,:);
%     end        
% end
% % Temperature modified Gibbs free energies of formation according to Van't Hoff's equation.
% 
% G0fTM = zeros(size(pTh.G0fM,1),size(pTh.G0fM,2),numGrd);
% for j=1:numGrd;
% G0fTM(:,:,j) = pTh.H0fM + (pPhys.T(j)/298).* (pTh.G0fM -pTh.H0fM);
% end
% % pTh.GfM = pTh.G0fM; % If no temperature dependence yet.
% % Calculation of the acid constants from the Gibbs energy formation matrix provided in thM for the liquid states. The thM matrix provide Gf of 
% % [Uncharged non hydrated species, (hydrated if aplicable) fully protonated species, 1st deprotonated species, 2nd deprotonated species and 3rd deprotonated species].
% % Ka are calculated from exp{ (Gºf_prot - Gºf_deprot)/(R·T)}.
% % Submatrix containing only the Gf of the liquid state variables.
% idGf_L = prod(1*(char(pTh.Type(1:size(ThParam,1)))==repmat('G0l',size(ThParam,1),1)),2)==1;
% GfM_L = zeros(sum(idGf_L),size(ThParam,2));         jGl = 0;
% for i=1:size(ThParam,1),
%     if     idGf_L(i)==1,
%         jGl = jGl+1;
%         G0fTM_L(jGl,:,:) = G0fTM(i,:,:);
%     end
% end
% 
% % Calculating acid equilibrium constants (3) for all layers at each temperatures, 3D matrix.
% for j=1:numGrd,
%     for i=1:size(G0fTM_L,2)-2,
%         KaM(:,i,j) = (G0fTM_L(:,i+2,j)~=0).*(exp((G0fTM_L(:,i+1,j)-G0fTM_L(:,i+2,j))/(pOp.Rth.*pPhys.T(j))));
%     end
% end
% pTh.KaM = KaM;
% % Calculating hydration constant (1) for all layers at each temperatures, 2D matrix.
% for j=1:numGrd,
%     KdM(:,j) = (GfM_L(:,1)~=0).*(exp((GfM_L(:,2)-GfM_L(:,1)-GfM_L(size(GfM_L,1),2))./(pOp.Rth*pPhys.T(j))));
% end
% pTh.KdM = KdM;
% clear ThParam pThNames jG jH jch jGl;

%% Structure with the THERMODYNAMIC PARAMETERS.
[ThParam, pThNames]= xlsread(FileModelXlsx, strcat('ThermoParam', id), '', 'basic');
aux = '';
for i=1:size(ThParam,1),
    % Building the string command for the structure, last without ', '.
    if i<size(ThParam,1),   
        aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:), '); 
    else
        aux = strcat(aux, char(39), pThNames(i), char(39), ', ThParam(', num2str(i), ',:)'); 
    end
end
aux = char(aux);
eval(strcat('pTh = struct(', aux, ');'));
pTh.pThNames = pThNames(:,1);
pTh.Type = pThNames(:,size(pThNames,2));

%Identify variables in function of the string label using strcmp function
idnSp = strcmp(pTh.Type,'nSp');
% Gibbs calculated with temperature
idG0tf = strcmp(pTh.Type,'G0ts')| strcmp(pTh.Type,'G0tl')| strcmp(pTh.Type,'G0tb')| strcmp(pTh.Type,'G0tg');
% Gibbs at standards conditions
idG0f  = strcmp(pTh.Type,'G0s') | strcmp(pTh.Type,'G0l') | strcmp(pTh.Type,'G0b') | strcmp(pTh.Type,'G0g');
% Enthalpies.
idH0f = strcmp(pTh.Type,'H0s')| strcmp(pTh.Type,'H0l')| strcmp(pTh.Type,'H0b')| strcmp(pTh.Type,'H0g');
% Charge identifier
idchr = strcmp(pTh.Type,'chr'); 
% Acid-base and dehydration constant identifier
idKda = strcmp(pTh.Type,'Kda'); 
%Reshape to include Ka and Kd for aqueous-phase components only
idL = (strcmp(R.St.Phase,'L'));

%Preallocation for the matrices
pTh.spcNames = [];  % Species names matrix
% Assigns species names to each cell of the matrix
for i=1:size(pThNames,1),
    if     idnSp(i)==1,
        pTh.spcNames = [pTh.spcNames; pThNames(i,2:size(pThNames,2)-1)];
    end
end


% pTh.G0tfM = zeros(sum(idG0tf),size(ThParam,2));  jGt = 0; %Gibbs calculated with temperature
pTh.G0fM  = zeros(sum(idG0f),size(ThParam,2));   jG = 0;
pTh.H0fM  = zeros(sum(idH0f),size(ThParam,2));   jH = 0;
pTh.chrM  = zeros(sum(idchr),size(ThParam,2));   jch = 0;
pTh.KdaM  = zeros(sum(idKda),size(ThParam,2));   jKda = 0;



%Assigns values to each of the matrices in function of their identifier:
%jGt - Gibbs with Temperature correction
%jG - Gibbs formation value at standard conditions
%jH - Enthalpy formation value
%jch - Charge parameter
%jKda - Acid-base constant
for i=1:size(ThParam,1),
    if  idG0tf(i)==1,
        jGt = jGt+1;
        pTh.G0tfM(jGt,:) = ThParam(i,:);
    elseif     idG0f(i)==1,
        jG = jG+1;
        pTh.G0fM(jG,:) = ThParam(i,:);
    elseif idH0f(i)==1,
        jH = jH+1;
        pTh.H0fM(jH,:) = ThParam(i,:);
    elseif idchr(i)==1,
        jch = jch+1;
        pTh.chrM(jch,:) = ThParam(i,:);       
    elseif idKda(i)==1,     
        jKda = jKda+1;
        pTh.KdaM(jKda,:) = ThParam(i,:);
    end  
end

%Calculate KaM, KdM, and Gibbs Free energy for each layer
for j=1:numGrd,
    [pTh.KaM(:,:,j), pTh.KdM(:,:,j), pTh.G0ftM(:,:,j)] = DG_Ka_calc(pTh, pOp, pPhys.T(j), idL);
end

% Building a matrix for selecting the species required for Gibbs
% calculations. An extra row of zeros is added to take into account the
% Gibbs formation of the proton.
spcSelM = zeros(length(pTh.G0fM(:,1)) + 1, length(pTh.G0fM(1,:)));

% The flag will assign a value of 1 to the species used in the
% Gibbs calculation while the rest will be maintained as 0 
for i=1:length(spcThFlag), 
    spcSelM(i, spcThFlag(i)) = 1; 
end
pTh.spcSelM = spcSelM;
clear ThParam pThNames jG jH jch jGl spcNames spcSelM spcThFlag;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Structure with the KINETIC PARAMETERS.
[KinetParam,  pKtNames]   = xlsread(FileModelXlsx, strcat('KinetParam',  id), '', 'basic');
aux = '';
for i=1:length(KinetParam),
    % Building the string command for the structure, last without ', '.
    if i<length(KinetParam),   
        aux = strcat(aux, char(39), pKtNames(i), char(39), ', ',num2str(KinetParam(i)),', '); 
    else
        aux = strcat(aux, char(39), pKtNames(i), char(39), ', ',num2str(KinetParam(i))); 
    end
end
aux = char(aux);
eval(strcat('pKt = struct(', aux, ');'));
pKt.pKtNames = pKtNames(:,1);
clear KinetParam pKtNames;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output to the global variable R.
R(idR).pOp   = pOp;     clear pOp;
R(idR).pTh   = pTh;     clear pTh;
R(idR).pPhys = pPhys;   clear pPhys;
R(idR).pKt   = pKt;     clear pKt;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('model_loaded.mat');
load('model_loaded.mat');
fprintf('DONE!');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
