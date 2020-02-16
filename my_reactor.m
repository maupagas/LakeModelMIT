% EnvBioProM Lab
% envbioprom.labs.masdar.ac.ae
% Chemical & Environmental Engineering
% Masdar Institute
% PO Box 54224 Abu Dhabi
% United Arab Emirates
% Contact: Jorge Rodríguez 
% jrodriguez@masdar.ac.ae
function [dStM_dt] = my_reactor(idR, R)
%% Preliminary calculations for gas balance.
% Calculating the total gas flow needed for the gas phase mass balances.
% Calculating the molar flow of gas generated [mol-gas/timeU].
% Sum of the generation term for gas transfer plus the gas flow fed.

% Moles of gas in the influent.
% NinT = (R(idR).Fd.Qgas_in*R(idR).pOp.Pgas)/(R(idR).pOp.Rg*R(idR).pOp.T); % molGas/h
NinT= 0;

% Moles of gas generated.  
ngas_tr = sum(bsxfun(@times,(char(R(idR).St.Phase) == 'G'), R(idR).RtM)) * R(idR).pOp.Vliq * R(idR).pOp.Vliq; % molGas/h

% Moles of water vapour.
nh2o_tr = (R(idR).pOp.Ph2o * (ngas_tr + NinT)) / (R(idR).pOp.Pgas - R(idR).pOp.Ph2o); % molVapour/h

% Moles of gas to get out.
NoutT = NinT + ngas_tr + nh2o_tr; % molGas/h

% Calculation of the gas flow.
Qgas = R(idR).pOp.Rg * R(idR).pPhys.T(1) * NoutT .* (NoutT > 0)  ./ R(idR).pOp.Pgas; % L/h%*(NoutT > 0)   (%Temperature in the top-most layer)
if NoutT < 0; R(idR).RtM(strcmp(R(idR).St.StNames, 'Gn2')) = -ngas_tr/R(idR).pOp.Vliq; end   %Why N2?

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Preliminary calculations for solids balance.
% Calculation of the dilution rate.
% D = R(idR).Fd.Qinf / R(idR).pOp.Vliq;    
HRT = R(idR).pOp.HRT;   D = 1/HRT;

% SRT might be specified different than HRT if solid retention is used.  
% If it is (-1) then SRT = HRT no solid retention otherwise the excel value given.
% SRT might be specified different than HRT if solid retention is used.  
    if     R(idR).pOp.SRT >=  0,          SRT = R(idR).pOp.SRT;    % If it is a positive value, there is solids retention
    elseif R(idR).pOp.SRT == -1,          SRT = HRT;               % If it is (1), then SRT = HRT
    elseif R(idR).pOp.SRT == -2,                                   % If it is (-2), there is solids retention that varies with time
         instSRT=0;
         for i=2:size(R(idR).SRTprofile,1),
             if R(idR).t >=  R(idR).SRTprofile(i-1,1) && R(idR).t < R(idR).SRTprofile(i,1) || R(idR).t > R(idR).SRTprofile(i,1)
                 instSRT=R(idR).SRTprofile(i,2);
             end
         end
         SRT=instSRT;
    else
        warning('Incorrect Choice Number for SRT mode')
    end

%% CALCULATION OF DERIVATIVES FOR ALL THE STATES.
dStM_dt = zeros(size(R(idR).StM));
dnGM_dt = zeros(size(R(idR).StM));
% For each state variable (Xi) the mass balance is applied according to the phase.
% If unchanged effluent concentration is state.
% InfV = R(idR).Fd.FdV(3:length(R(idR).Fd.FdV));              % kmol/m3 - First two elements are Q and Qgas
EffV = R(idR).StM;


% Depending on the type of state variable, a different differential equation may apply.
% Calculating total solids concentration from previuos time step. It
% calculates solids by adding up solid phases 'S' (i.e. precipitations) and
% biomass 'B'
% Xtot = sum((char(R(idR).St.Phase) == 'S' | char(R(idR).St.Phase) == 'B').* EffV);

% Depending on the type of state variable, a different differential equation may apply.
for i=1:R(idR).Dim.numStG,
    Phase = char(R(idR).St.Phase(i));
    switch Phase,
        case 'L'        % If it is a dissolved species.
            dStM_dt(i,:) = R(idR).RtM(i,:) + D * (R(idR).FlwM(i,:)-R(idR).StM(i,:));
        case 'B'        % If it is a biomass species.
            dStM_dt(i,:) = R(idR).RtM(i,:);
        case 'S'        % If it is a solids species.
            dStM_dt(i,:) = R(idR).RtM(i,:);
        case 'G'        % If it is a gas species.
            dStM_dt(i,:) = R(idR).RtM(i,:);
        otherwise
            dStM_dt(i,:) = 0;
    end
end


% for i=1:R(idR).Dim.numSt,
%     Phase = char(R(idR).St.Phase(i));
%     switch Phase,
%         case 'L'        % soluble phase state variable.
%             % In liquid phase the outlet of the reactor is the current concentration.
%             dStM_dt(i) = D*(InfV(i) - EffV(i)) + R(idR).RtM(i);
% 
%         case 'S'        %  solid phase state variable.
%             % In solid phase the outlet is defined by the SRT/HRT ratio.
%                       if  HRT==Inf || HRT==0 || SRT==-2,     
%                           EffV(i) = R(idR).StM(i);            
%                       else
%                           % To limit the maximum amount of solids (Xtmax) in the reactor
%                           if  R(idR).pOp.Xtmax ~= 0,
%                               if Xtot < R(idR).pOp.Xtmax*1.01,          % if Xtot is within Xtmax, the concentrations are coded to be unaffected
%                                   fx_out = 1 + 0.05*Xtot / abs(R(idR).pOp.Xtmax - Xtot);
%                               else        % if Xtot exceeds Xmax, solids are forced to exit the reactor
%                                   fx_out=SRT/HRT;
%                               end
%                               EffV(i) = fx_out*(R(idR).StM(i) * HRT / SRT);
%                           else
%                               % if there is no Xtmax limitation
%                               EffV(i) = R(idR).StM(i) * HRT / SRT;  
%                           end
%                       end
%                       % Mass balance equation for solids in a CSTR(idR).
%                       dStM_dt(i) = D*(InfV(i) - EffV(i)) + R(idR).RtM(i);
%        
%         case 'B'        %  biomass state variable.
%             % In solid phase the outlet is defined by the SRT/HRT ratio.
%                       if HRT==Inf || HRT==0 || SRT==-2,     
%                           EffV(i) = R(idR).StM(i);            
%                       else
%                           % To limit the maximum amount of solids (Xtmax) in the reactor
%                           if R(idR).pOp.Xtmax~=0
%                               if Xtot<R(idR).pOp.Xtmax*1.01          % if Xtot is within Xtmax, the concentrations are coded to be unaffected
%                                   fx_out = 1 + 0.05*Xtot / abs(R(idR).pOp.Xtmax - Xtot);
%                               else        % if Xtot exceeds Xmax, solids are forced to exit the reactor
%                                   fx_out=SRT/HRT;
%                               end
%                               EffV(i) = fx_out*(R(idR).StM(i) * HRT / SRT);
%                           else
%                               % if there is no Xtmax limitation
%                               EffV(i) = R(idR).StM(i) * HRT / SRT;  
%                           end
%                       end
%                       % Mass balance equation for solids in a CSTR(idR).
%                       dStM_dt(i) = D*(InfV(i) - EffV(i)) + R(idR).RtM(i);
%         
%         case 'G'        %  gas phase state variable. In the model, they are calculated as bar
%             % For gas phase the mass balance considers the Qgas calculated.
%             % Reactor units are in pressure units (bar) so it needs to be
%             % converted to moles dividing by R * T. Moles can be used to
%             % check the balance and will be stored in dnGM_dt.
%             InfV(i) = (InfV(i)*R(idR).Fd.Qgas_in)/(R(idR).pOp.Rg*R(idR).pOp.T);                       %    moles/h 
%             EffV(i) = EffV(i)*Qgas/(R(idR).pOp.Rg*R(idR).pOp.T);                                      %    moles/h
%            
%             dnGM_dt(i)=  (InfV(i)-EffV(i))+ (R(idR).RtM(i)*R(idR).pOp.Vliq);                          %    moles/h
%             dStM_dt(i)=  dnGM_dt(i)*R(idR).pOp.Rg*R(idR).pOp.T/R(idR).pOp.Vgas;                       %    bar/h
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cancelling derivatives of states defined as non dynamic blocked.
dStM_dt = bsxfun(@times, dStM_dt, R(idR).St.Dyn);
