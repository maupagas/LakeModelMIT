% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae

function [rRc,rTr,rRcM,rTrM] = my_kinetics(idR, R)


% Reading all parameters and variables for this reactor.
St    = R(idR).St;          StM    = R(idR).StM;       
AlgSt = R(idR).AlgSt;       pOp   = R(idR).pOp;         
pKt    = R(idR).pKt;        pPhys = R(idR).pPhys;          
rTr   = R(idR).rTr;         rRc    = R(idR).rRc; 

% Compute number of grids and each variable at each phase
numGrd = size(StM,2);
numStG_S = sum(1*strcmp(St.Phase,'L'));     % Number of dissolved states.
numStG_B = sum(1*strcmp(St.Phase,'B'));     % Number of active biomass states.
numStG_X = sum(1*strcmp(St.Phase,'S'));     % Number of solid states.
numStG_G = sum(1*strcmp(St.Phase,'G'));     % Number of gas states.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CUSTOM AREA FOR THE MODEL BEING IMPLEMENTED
% The equations calculating the current rates of reaction must be written below.
% The variables refered must be taken from in the right structure they belong to.
% Fd:feeding, St:state variable, pOp:operative, pKt:kinetic
% Rates must be in the same order than in the stoichiometry matrix!!!

% %%%%%%%%%%%%%%%% TRANSPORT KINETICS %%%%%%%%%%%%
    % Diffussions function of temperature in Celsius at each layer.
    Tc = pPhys.T - 273.15;
    DiffSglu = pKt.sDiffSglu.*Tc + pKt.oDiffSglu;       DiffSac  = pKt.sDiffSac .*Tc + pKt.oDiffSac;
    DiffSch4 = pKt.sDiffSch4.*Tc + pKt.oDiffSch4;       DiffSic  = pKt.sDiffSic .*Tc + pKt.oDiffSic;
    DiffSh2  = pKt.sDiffSh2 .*Tc + pKt.oDiffSh2;        DiffSamn = pKt.sDiffSamn.*Tc + pKt.oDiffSamn;
    DiffSno2 = pKt.sDiffSno2.*Tc + pKt.oDiffSno2;       DiffSno3 = pKt.sDiffSno3.*Tc + pKt.oDiffSno3;
    DiffSn2  = pKt.sDiffSn2 .*Tc + pKt.oDiffSn2;        DiffSh2s = pKt.sDiffSh2s.*Tc + pKt.oDiffSh2s;
    DiffSso4 = pKt.sDiffSso4.*Tc + pKt.oDiffSso4;       DiffSfe2 = pKt.sDiffSfe2.*Tc + pKt.oDiffSfe2;
    DiffSfe3 = pKt.sDiffSfe3.*Tc + pKt.oDiffSfe3;       DiffSo2  = pKt.sDiffSo2 .*Tc + pKt.oDiffSo2;        
    DiffScat = pKt.sDiffScat.*Tc + pKt.oDiffScat;       DiffSan  = pKt.sDiffSan .*Tc + pKt.oDiffSan;
    DiffXb   = pKt.sDiffXb  .*Tc + pKt.oDiffXb;
    VsX_ = [pKt.VsXd; pKt.VsXfe; pKt.VsXi];
    VsG_ = pKt.VsG;
    
    % VERY IMPORTANT: THIS MUST REMAIN IN THE SAME ORDER THAN THE STATES IN
    % THE TRANSPORT MATRIX. ******(SHOULD BE CHECKED AT THE END #MPG)***********************
    DiffStL = [DiffSglu; DiffSac; DiffSch4; DiffSic; DiffSh2; DiffSamn; DiffSno2; DiffSno3; DiffSn2; DiffSh2s;...
        DiffSso4; DiffSfe2; DiffSfe3; DiffSo2; DiffScat; DiffSan];
    % Turbulence factors for increased diffussivity.
        for i=1:numStG_S,
    fmDiff(i,:) = pPhys.fmDiff;
        end
    DiffStL = fmDiff.* DiffStL;
    % Biomass brownian type of diffusion.
        for i=1:numStG_B,
    DiffXbB(i,:) = DiffXb;
        end
    % Solids settling velocities.
        for j=1:numGrd,
    VsX(:,j) = VsX_;
        end
    % Bubbles ascension velocities (NOT USED ANY LONGER)
        for j=1:numGrd,
    VsG(1:numStG_G,j) = VsG_;
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pOp.rTr == 1,      %Executes Transport Reactions if it is activated (1) or neglects it (0)
   
    % Transport submatrix for SOLUBLE species from UPPER layer.
    rTrMsT = zeros(numStG_S,numGrd);
    dCs= zeros(numStG_S, numGrd);   dz = zeros(numStG_S, numGrd);   dL = zeros(numStG_S, numGrd);
    for i=2:numGrd;
        dCs(:,i) = StM(1:numStG_S,i-1) - StM(1:numStG_S,i);
        dz(:,i)  = pPhys.z(i)-pPhys.z(i-1);    
        if i>1,    dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end
    end
    rTrMsT_ = (DiffStL./dL) .* (dCs./dz);
    rTrMsT(:,2:numGrd) = rTrMsT_(:,2:numGrd);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for SOLUBLE species from LOWER layer.
    rTrMsB = zeros(numStG_S,numGrd);
    dCs= zeros(numStG_S, numGrd);   dz = zeros(numStG_S, numGrd);   dL = zeros(numStG_S, numGrd);
    for i=1:numGrd-1;
        dCs(:,i) = StM(1:numStG_S,i+1) - StM(1:numStG_S,i);
        dz(:,i)  = pPhys.z(i+1)-pPhys.z(i);     
        if i>1,    dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end
    end
    rTrMsB_ = (DiffStL./dL) .* (dCs./dz);
    rTrMsB(:,1:numGrd-1) = rTrMsB_(:,1:numGrd-1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for active BIOMASS species from UPPER layer.
    rTrMbT = zeros(numStG_B,numGrd);
    dCb= zeros(numStG_B, numGrd);   dz = zeros(numStG_B, numGrd);   dL = zeros(numStG_B, numGrd);
    for i=2:numGrd;
        dCb(:,i) = StM(numStG_S+1:numStG_S+numStG_B,i-1) - StM(numStG_S+1:numStG_S+numStG_B,i);
        dz(:,i)  = pPhys.z(i)-pPhys.z(i-1);     
        if i>1,    dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end    
    end
    rTrMbT_ = (DiffXbB./dL) .* (dCb./dz);
    rTrMbT(:,2:numGrd) = rTrMbT_(:,2:numGrd);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for active BIOMASS species from LOWER layer.
    rTrMbB = zeros(numStG_B,numGrd);
    dCb= zeros(numStG_B, numGrd);   dz = zeros(numStG_B, numGrd);   dL = zeros(numStG_B, numGrd);
    for i=1:numGrd-1;
        dCb(:,i)= StM(numStG_S+1:numStG_S+numStG_B,i+1) - StM(numStG_S+1:numStG_S+numStG_B,i);
        dz(:,i) = pPhys.z(i+1)-pPhys.z(i);      
        if i>1,    dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end      
    end
    rTrMbB_ = (DiffXbB./dL) .* (dCb./dz);
    rTrMbB(:,1:numGrd-1) = rTrMbB_(:,1:numGrd-1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for settling SOLID species from UPPER layer.
    rTrMxT = zeros(numStG_X,numGrd);
    Cx= zeros(numStG_X, numGrd);   dL = zeros(numStG_X, numGrd);
    for i=2:numGrd;
        Cx(:,i)  = StM(numStG_S+numStG_B+1:numStG_S+numStG_B+numStG_X,i-1); 
        if i>1,    dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end
    end
    rTrMxT_ = (VsX ./ dL).*Cx;
    rTrMxT(:,2:numGrd) = rTrMxT_(:,2:numGrd);
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for active settling SOLID species from (to since it is <0) LOWER layer.
    rTrMxB = zeros(numStG_X,numGrd);
    Cx= zeros(numStG_X, numGrd);   dL = zeros(numStG_X, numGrd);
    for i=1:numGrd-1;
        Cx(:,i)= StM(numStG_S+numStG_B+1:numStG_S+numStG_B+numStG_X,i);   
        if i>1,   dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end      
    end
    rTrMxB_ = -(VsX ./ dL).*Cx;
    rTrMxB(:,1:numGrd-1) = rTrMxB_(:,1:numGrd-1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for rising BUBBLE GAS species from LOWER layer.
    rTrMgB = zeros(numStG_G,numGrd);
    Cg= zeros(numStG_G, numGrd);   dL = zeros(numStG_G, numGrd);
    for i=1:numGrd-1;
        Cg(:,i)= StM(numStG_S+numStG_B+numStG_X+1:numStG_S+numStG_B+numStG_X+numStG_G,i+1);   if i>1,    
        dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end      
    end
    rTrMgB_ = (VsG ./ dL).*Cg;
    rTrMgB(:,1:numGrd-1) = rTrMgB_(:,1:numGrd-1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transfer rates of SOLUBLE species TO GAS PHASE. (IS THIS ONLY UPPER LAYER?? #MPG)

    %Converting Matrix speciation Data into vectors 
    GgasT = (AlgSt.Vgas .* pPhys.P) ./ (pOp.Rg .* pPhys.T);   % molG/Lliq
    rTrSch4G = pPhys.kLa.*(St.Sch4    - (St.Gch4./GgasT) .* pPhys.P * pKt.Kh_ch4);
    rTrSco2G = pPhys.kLa.*(AlgSt.Sco2 - (St.Gco2./GgasT) .* pPhys.P * pKt.Kh_co2);    %(1) as subindex for being surface concentration of CO2
    rTrSh2G  = pPhys.kLa.*(St.Sh2     - (St.Gn2 ./GgasT) .* pPhys.P * pKt.Kh_n2);
    rTrSn2G  = pPhys.kLa.*(St.Sn2     - (St.Gh2 ./GgasT) .* pPhys.P * pKt.Kh_n2);
    rTrSsh2G = pPhys.kLa.*(AlgSt.Ssh2       - (St.Gsh2./GgasT) .* pPhys.P * pKt.Kh_sh2);    %(1) as subindex for being surface concentration of H2S
    rTrSo2G  = pPhys.kLa.*(St.So2     - (St.Go2 ./GgasT) .* pPhys.P * pKt.Kh_o2);
    rTrMsG = [rTrSch4G; rTrSco2G; rTrSh2G; rTrSn2G; rTrSsh2G; rTrSo2G];
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transport submatrix for rising BUBBLE GAS species from (to since <0) UPPER layer.
    rTrMgT = zeros(numStG_G,numGrd);
    Cg= zeros(numStG_G, numGrd);   dL = zeros(numStG_G, numGrd);
    for i=1:numGrd;
        Cg(:,i)  = StM(numStG_S+numStG_B+numStG_X+1:numStG_S+numStG_B+numStG_X+numStG_G,i); 
        if i>1, dL(:,i)  = pPhys.L(i)- pPhys.L(i-1);    else     dL(:,i)  = pPhys.L(i); end
    end
    rTrMgT_ = (VsG ./ dL).*Cg;
    rTrMgT(:,1:numGrd) = rTrMgT_(:,1:numGrd);
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gas bubbles removal rate total EQUALISED to transfer rate by using a factor.
    f_rTrMgT = sum(rTrMsG)./sum(rTrMgT);
    for i=1:numStG_G,
        rTrMgT(i,:) = f_rTrMgT.*rTrMgT(i,:);
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complete transport matrix built from all submatrixes.
    rTrM = [rTrMsT; rTrMsB; rTrMbT; rTrMbB; rTrMxT; rTrMxB; rTrMgT; rTrMgB; rTrMsG];

else
    rTrM = zeros(length(rTr.TrNames), numGrd);
end
%Giving values matrix of each state variable. 
% for i=1:length(rTr.TrNames),
%     eval(strcat('rTr.', char(rTr.TrNames(i)), ' = [', char(num2str(rTrM(i,:))), '];'));
% end

%Alternative to eval (MUCH FASTER #MPG)
for i = 1:length(rTr.TrNames), rTr.(char(rTr.TrNames{i})) = rTrM(i,:); end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% END OF TRANSPORT KINETICS %%%%%%%%%%%%

% %%%%%%%%%%%%%%%% REACTION KINETICS %%%%%%%%%%%%
% Expanding Ks into a vector for all layers. (No need to do this unless Ks
% is going to be function of something else #MPG)
% Ks for electron donors.

if pOp.rRc == 1,   %Calculates kinetics if they are activated in Excel (1). Otherwise, ignores it.
    
    KsNcyb  = pKt.KsNcyb;         KsDfer  = pKt.KsDfer;
    KsDhglu = pKt.KsDhglu;        KsDhac  = pKt.KsDhac;
    KsDaob  = pKt.KsDaob;         KsDnob  = pKt.KsDnob;
    KsDdnglu= pKt.KsDdnglu;       KsDdnac = pKt.KsDdnac;
    KsDdnh2 = pKt.KsDdnh2;
    KsDsrd  = pKt.KsDsrd;         KsDsrdh2= pKt.KsDsrdh2;
    KsDsox  = pKt.KsDsox;         KsDson  = pKt.KsDson;         
    KsDfeox = pKt.KsDfeox;        KsDfeon = pKt.KsDfeon;
    KsDferd = pKt.KsDferd;        KsDferdh2=pKt.KsDferdh2;
    KsDacm  = pKt.KsDacm;         KsDh2m  = pKt.KsDh2m;
    KsDmto  = pKt.KsDmto;         KsDmts  = pKt.KsDmts;
    % Ks for electron acceptors.
    KsCcyb  = pKt.KsCcyb;         KsAfer  = pKt.KsAfer;
    KsAhglu = pKt.KsAhglu;        KsAhac  = pKt.KsAhac;
    KsAaob  = pKt.KsAaob;         KsAnob  = pKt.KsAnob;
    KsAdnglu= pKt.KsAdnglu;       KsAdnac = pKt.KsAdnac;
    KsAdnh2 = pKt.KsAdnh2;
    KsAsrd  = pKt.KsAsrd;         KsAsrdh2= pKt.KsAsrdh2;
    KsAsox  = pKt.KsAsox;         KsAson  = pKt.KsAson;
    KsAfeox = pKt.KsAfeox;        KsAfeon = pKt.KsAfeon;
    KsAferd = pKt.KsAferd;        KsAferdh2=pKt.KsAferdh2;
    KsAacm  = pKt.KsAacm;         KsAh2m  = pKt.KsAh2m;
    KsAmto  = pKt.KsAmto;         KsAmts  = pKt.KsAmts;
    % Ks for nirogen source and inorganic carbon.
    KsAMNs = pKt.KsAMNs;            KsICs = pKt.KsICs;

    % Decay and hydrolysis constants may become temperature dependent for each layer.
    kdXb = pKt.kdXb * pPhys.T ./ pPhys.T;   khXd = pKt.khXd * pPhys.T ./ pPhys.T;

    % Monod terms for eD, eA, N source and C source (if applicable, otherwise 1) for all metabolisms.
    % Monods for electron donors.
    Md_cyb  = ones(1,numGrd);                   Md_fer  = St.Sglu ./ (KsDfer  + St.Sglu);
    Md_hglu = St.Sglu ./ (KsDhglu + St.Sglu);   Md_hac  = St.Sac  ./ (KsDhac  + St.Sac);
    Md_aob  = St.Samn ./ (KsDaob  + St.Samn);   Md_nob  = St.Sno2 ./ (KsDnob  + St.Sno2);
    Md_dnglu= St.Sglu ./ (KsDdnglu+ St.Sglu);   Md_dnac = St.Sac  ./ (KsDdnac + St.Sac);
    Md_dnh2 = St.Sh2  ./ (KsDdnh2 + St.Sh2);
    Md_srd  = St.Sac  ./ (KsDsrd  + St.Sac);    Md_srdh2= St.Sh2  ./ (KsDsrdh2+ St.Sh2);
    Md_sox  = St.Sh2s ./ (KsDsox  + St.Sh2s);   Md_son  = St.Sh2s ./ (KsDson  + St.Sh2s);   
    Md_feox = St.Sfe2 ./ (KsDfeox + St.Sfe2);   Md_feon = St.Sfe2 ./ (KsDfeon + St.Sfe2);   
    Md_ferd = St.Sac  ./ (KsDferd + St.Sac);    Md_ferdh2=St.Sh2  ./ (KsDferdh2+St.Sh2);
    Md_acm  = St.Sac  ./ (KsDacm  + St.Sac);    Md_h2m  = St.Sh2  ./ (KsDh2m  + St.Sh2);
    Md_mto  = St.Sch4 ./ (KsDmto  + St.Sch4);   Md_mts  = St.Sch4 ./ (KsDmts  + St.Sch4);
    % Monods for electron acceptors.
    Ma_cyb  = ones(1,numGrd);                   Ma_fer  = ones(1,numGrd);
    Ma_hglu = St.So2  ./ (KsAhglu + St.So2);    Ma_hac  = St.So2  ./ (KsAhac  + St.So2);
    Ma_aob  = St.So2  ./ (KsAaob  + St.So2);    Ma_nob  = St.So2  ./ (KsAnob  + St.So2);
    Ma_dnglu= St.Sno3 ./ (KsAdnglu+ St.Sno3);   Ma_dnac = St.Sno3 ./ (KsAdnac + St.Sno3);
    Ma_dnh2 = St.Sno3 ./ (KsAdnh2 + St.Sno3);
    Ma_srd  = St.Sso4 ./ (KsAsrd  + St.Sso4);   Ma_srdh2= St.Sso4 ./ (KsAsrdh2+ St.Sso4);
    Ma_sox  = St.So2  ./ (KsAsox  + St.So2);    Ma_son  = St.Sno3 ./ (KsAson  + St.Sno3);
    Ma_feox = St.So2  ./ (KsAfeox + St.So2);    Ma_feon = St.Sno3 ./ (KsAfeon + St.Sno3);   
    Ma_ferd = St.Xfeoh;                         Ma_ferdh2=St.Xfeoh;
    Ma_acm  = ones(1,numGrd);                   Ma_h2m  = St.Sic  ./ (KsAh2m  + St.Sic);
    Ma_mto  = St.So2  ./ (KsAmto  + St.So2);    Ma_mts  = St.Sso4 ./ (KsAmts  + St.Sso4);
    % Monods for nitrogen source (some have N already as eD or eA)
    Mn_cyb  = St.Samn ./ (KsNcyb + St.Samn);    Mn_fer  = St.Samn ./ (KsAMNs + St.Samn);
    Mn_hglu = St.Samn ./ (KsAMNs + St.Samn);    Mn_hac  = St.Samn ./ (KsAMNs + St.Samn);
    Mn_aob  = ones(1,numGrd);                   Mn_nob  = ones(1,numGrd);
    Mn_dnglu= ones(1,numGrd);                   Mn_dnac = ones(1,numGrd);
    Mn_dnh2 = ones(1,numGrd);
    Mn_srd  = St.Samn ./ (KsAMNs + St.Samn);    Mn_srdh2= St.Samn ./ (KsAMNs + St.Samn);
    Mn_sox  = St.Samn ./ (KsAMNs + St.Samn);    Mn_son  = ones(1,numGrd);
    Mn_feox = St.Samn ./ (KsAMNs + St.Samn);    Mn_feon = ones(1,numGrd);
    Mn_ferd = St.Samn ./ (KsAMNs + St.Samn);    Mn_ferdh2=St.Samn ./ (KsAMNs + St.Samn);
    Mn_acm  = St.Samn ./ (KsAMNs + St.Samn);    Mn_h2m  = St.Samn ./ (KsAMNs + St.Samn);
    Mn_mto  = St.Samn ./ (KsAMNs + St.Samn);    Mn_mts  = St.Samn ./ (KsAMNs + St.Samn);
    % Monods for autotrophic carbon source if no C is present in either eD or eA.
    Mc_cyb  = St.Sic ./ (KsCcyb + St.Sic);      Mc_fer  = ones(1,numGrd);
    Mc_hglu = ones(1,numGrd);                   Mc_hac  = ones(1,numGrd);
    Mc_aob  = St.Sic ./ (KsICs + St.Sic);       Mc_nob  = St.Sic ./ (KsICs + St.Sic);
    Mc_dnglu= ones(1,numGrd);                   Mc_dnac = ones(1,numGrd);
    Mc_dnh2 = ones(1,numGrd);
    Mc_srd  = ones(1,numGrd);                   Mc_srdh2= ones(1,numGrd);
    Mc_sox  = St.Sic ./ (KsICs + St.Sic);       Mc_son  = St.Sic ./ (KsICs + St.Sic);
    Mc_feox = St.Sic ./ (KsICs + St.Sic);       Mc_feon = St.Sic ./ (KsICs + St.Sic);
    Mc_ferd = ones(1,numGrd);                   Mc_ferdh2=ones(1,numGrd);
    Mc_acm  = ones(1,numGrd);                   Mc_h2m  = ones(1,numGrd);
    Mc_mto  = ones(1,numGrd);                   Mc_mts  = ones(1,numGrd);
    % Specific inhibition terms if needed they can be expanded here (e.g. by pH)
    % INHIBITIONS
    pH = AlgSt.pH;      Ki_nh3  = pKt.Ki_nh3;  %Is this necessary ?? (#MPG)
    % General pH inhibition.
    Iph = ones(size(pH)); % All ones.
    Iph = Iph + (pH <= pKt.pHlo)   .* ((exp(-3*((pH-pKt.pHlo)   /(pKt.pHlo   -pKt.pHll)).^2))-1);
    Iph = Iph + (pH >= pKt.pHuo)   .* ((exp(-3*((pH-pKt.pHuo)   /(pKt.pHuo   -pKt.pHul)).^2))-1);

    %Iph    = 1 + (pH < pKt.pHul)   .* ((exp(-3*((pH-pKt.pHul)   /(pKt.pHul   -pKt.pHll)).^2))-1);

    % Aceticlastic methanogens pH inhibition.
    Iph_ac = 1 + (pH < pKt.pHul_ac).* ((exp(-3*((pH-pKt.pHul_ac)/(pKt.pHul_ac-pKt.pHll_ac)).^2))-1);

    % Hydrogenotroph Methanogens pH inhibition.
    Iph_h2 = 1 + (pH < pKt.pHul_h2).* ((exp(-3*((pH-pKt.pHul_h2)/(pKt.pHul_h2-pKt.pHll_h2)).^2))-1);

    % Aceticlastic  methanogens free amonia inhibition.
    Inh3_xacm = Ki_nh3 ./ (Ki_nh3 + AlgSt.Snh3);

    % Competitive terms when one biomass groups targets multiple substrates.
    Ic_hglu = (St.Sglu ./(St.Sglu  + St.Sac));
    Ic_hac  = (St.Sac  ./(St.Sglu  + St.Sac));
    Ic_dnglu= (St.Sglu ./(St.Sglu  + St.Sac + St.Sh2));
    Ic_dnac = (St.Sac  ./(St.Sglu  + St.Sac + St.Sh2));
    Ic_dnh2 = (St.Sh2  ./(St.Sglu  + St.Sac + St.Sh2));
    Ic_srd  = (St.Sac  ./(St.Sac   + St.Sh2));
    Ic_srdh2= (St.Sh2  ./(St.Sac   + St.Sh2));
    Ic_ferd = (St.Sac  ./(St.Sac   + St.Sh2));
    Ic_ferdh2=(St.Sh2  ./(St.Sac   + St.Sh2));

    % Specific inhibitions for each microbial group.
    Icyb  = Iph.*pPhys.phv;         Ifer   = Iph;
    Ihglu = Iph.*Ic_hglu;           Ihac   = Iph.*Ic_hac;
    Iaob  = Iph;                    Inob   = Iph;
    Idnglu= Iph.*Ic_dnglu;          Idnac  = Iph.*Ic_dnac;
    Idnh2 = Iph.*Ic_dnh2;
    Isrd  = Iph.*Ic_srd;            Isrdh2 = Iph.*Ic_srdh2;
    Isox  = Iph;                    Ison   = Iph;
    Ifeox = Iph;                    Ifeon  = Iph;
    Iferd = Iph.*Ic_ferd;           Iferdh2= Iph.*Ic_ferdh2;
    Iacm  = Iph_ac.*Inh3_xacm;      Ih2m   = Iph.*Iph_h2;
    Imto  = Iph;                    Imts   = Iph;

    % Metabolism reaction kinetics.
    rRc.rMcyb  = pKt.qmcyb  .* Md_cyb  .* Ma_cyb  .* Mn_cyb  .* Mc_cyb  .* Icyb  .* St.Xcyb;
    rRc.rMfer  = pKt.qmfer  .* Md_fer  .* Ma_fer  .* Mn_fer  .* Mc_fer  .* Ifer  .* St.Xfer;
    rRc.rMhglu = pKt.qmhglu .* Md_hglu .* Ma_hglu .* Mn_hglu .* Mc_hglu .* Ihglu .* St.Xhet;    % Competing term for eD used.
    rRc.rMhac  = pKt.qmhac  .* Md_hac  .* Ma_hac  .* Mn_hac  .* Mc_hac  .* Ihac  .* St.Xhet;    % Competing term for eD used.
    rRc.rMaob  = pKt.qmaob  .* Md_aob  .* Ma_aob  .* Mn_aob  .* Mc_aob  .* Iaob  .* St.Xaob;
    rRc.rMnob  = pKt.qmnob  .* Md_nob  .* Ma_nob  .* Mn_nob  .* Mc_nob  .* Inob  .* St.Xnob;
    rRc.rMdnglu= pKt.qmdnglu.* Md_dnglu.* Ma_dnglu.* Mn_dnglu.* Mc_dnglu.* Idnglu.* St.Xdn;     % Competing term for eD used.
    rRc.rMdnac = pKt.qmdnac .* Md_dnac .* Ma_dnac .* Mn_dnac .* Mc_dnac .* Idnac .* St.Xdn;     % Competing term for eD used.
    rRc.rMdnh2 = pKt.qmdnh2 .* Md_dnh2 .* Ma_dnh2 .* Mn_dnh2 .* Mc_dnh2 .* Idnh2 .* St.Xdn;     % Competing term for eD used.
    rRc.rMsrd  = pKt.qmsrd  .* Md_srd  .* Ma_srd  .* Mn_srd  .* Mc_srd  .* Isrd  .* St.Xsrd;    % Competing term for eD used.
    rRc.rMsrdh2= pKt.qmsrdh2.* Md_srdh2.* Ma_srdh2.* Mn_srdh2.* Mc_srdh2.* Isrdh2.* St.Xsrd;    % Competing term for eD used.
    rRc.rMsox  = pKt.qmsox  .* Md_sox  .* Ma_sox  .* Mn_sox  .* Mc_sox  .* Isox  .* St.Xsox;
    rRc.rMson  = pKt.qmson  .* Md_son  .* Ma_son  .* Mn_son  .* Mc_son  .* Ison  .* St.Xson;
    rRc.rMfeox = pKt.qmfeox .* Md_feox .* Ma_feox .* Mn_feox .* Mc_feox .* Ifeox .* St.Xfeox;
    rRc.rMfeon = pKt.qmfeon .* Md_feon .* Ma_feon .* Mn_feon .* Mc_feon .* Ifeon .* St.Xfeon;
    rRc.rMferd = pKt.qmferd .* Md_ferd .* Ma_ferd .* Mn_ferd .* Mc_ferd .* Iferd .* St.Xferd;   % Competing term for eD used.
    rRc.rMferdh2=pKt.qmferdh2.*Md_ferdh2.*Ma_ferdh2.*Mn_ferdh2.*Mc_ferdh2.*Iferdh2.*St.Xferd;   
    rRc.rMacm  = pKt.qmacm  .* Md_acm  .* Ma_acm  .* Mn_acm  .* Mc_acm  .* Iacm  .* St.Xacm;
    rRc.rMh2m  = pKt.qmh2m  .* Md_h2m  .* Ma_h2m  .* Mn_h2m  .* Mc_h2m  .* Ih2m  .* St.Xh2m;
    rRc.rMmto  = pKt.qmmto  .* Md_mto  .* Ma_mto  .* Mn_mto  .* Mc_mto  .* Imto  .* St.Xmto;
    rRc.rMmts  = pKt.qmmts  .* Md_mts  .* Ma_mts  .* Mn_mts  .* Mc_mts  .* Imts  .* St.Xmts;

    % Precipitation rates. Set up such that does not exceed a max Fe3+ in solution. [M/h]
    rRc.rPrFe3 = pKt.Kprfe3 * ( (St.Sfe3>=pKt.Sfe3max).*St.Sfe3 + (St.Sfe3<pKt.Sfe3max).*St.Xfeoh  ) .* ((St.Sfe3./pKt.Sfe3max) - 1);

    % Biomass decay kinetics.
    rRc.rXcyb  = kdXb .* St.Xcyb;
    rRc.rXfer  = kdXb .* St.Xfer;
    rRc.rXhet  = kdXb .* St.Xhet;
    rRc.rXaob  = kdXb .* St.Xaob;
    rRc.rXnob  = kdXb .* St.Xnob;
    rRc.rXdn   = kdXb .* St.Xdn;
    rRc.rXsrd  = kdXb .* St.Xsrd;
    rRc.rXsox  = kdXb .* St.Xsox;
    rRc.rXson  = kdXb .* St.Xson;
    rRc.rXfeox = kdXb .* St.Xfeox;
    rRc.rXfeon = kdXb .* St.Xfeon;
    rRc.rXferd = kdXb .* St.Xferd;
    rRc.rXacm  = kdXb .* St.Xacm;
    rRc.rXh2m  = kdXb .* St.Xh2m;
    rRc.rXmto  = kdXb .* St.Xmto;
    rRc.rXmts  = kdXb .* St.Xmts;
    % Hydrolysis rate of dead biomass.
    rRc.rXd    = khXd .* St.Xd;


% %%%%%%%%%%%%%%%% END OF REACTION KINETICS %%%%%%%%%%%%

    % Matrix of reaction rates for each layer is created from named variables.
    rRcM = zeros(length(rRc.RcNames),numGrd);
    for i=1:length(rRc.RcNames),        
        rRcM(i,:) = rRc.(char(rRc.RcNames(i)));    
    end
    
else
    rRcM = zeros(length(rRc.RcNames),numGrd);
end

% END OF THE CUSTOM AREA
%% DOUBLE-CHECKING ORDER OF WRITTEN KINETICS IN BOTH EXCEL AND MATLAB 
%IT IS DONE ONLY AT t=0, AT FIRST TIME MY_KINETICS IS RUN
if R(idR).t == 0 && (pOp.rRc == 1 && pOp.rTr == 1),
      
    %%%Reaction kinetics 
    %Obtaining the order written in this script
    rRc_Matlab = fieldnames(rRc);
    rRc_Matlab = rRc_Matlab(2:length(rRc_Matlab));
    %Obtaining the order written in Excel
    rRc_Excel = rRc.RcNames;

    %Comparison of both vectors
    if any(strcmp(rRc_Matlab, rRc_Excel) == 0)
        error('REACTION KINETICS are not written in a proper order. Please, double-check and re-execute the model...\n')
    end

    %Transport kinetics
    %Obtaining the order written in this script
    rTr_Matlab = fieldnames(rTr);
    rTr_Matlab = rTr_Matlab(2:length(rTr_Matlab));
    %Obtaining the order located in excel
    rTr_Excel = rTr.TrNames;

    %Comparison of both vectors
    if any(strcmp(rTr_Matlab, rTr_Excel) == 0)
        error('TRANSPORT KINETICS are not written in a proper order. Please, double-check and re-execute the model...\n')
    end

end 

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
