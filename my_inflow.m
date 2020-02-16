% AWMC-MFC Model implementation by
%
% Dr. Jorge Rodríguez Rodríguez
% Advanced Wastewater Management Centre (AWMC)
% The University of Queensland
% L4/Gehrmann Bld. (60), Research dRd.
% Brisbane QLD 4072
% AUSTRALIA
%
% Contact: jorger@usc.es
function [Fd] = my_feeding(t, u, idR, Conn, mFd, R);

% Feed flow and concentrations are given in a matrix as well as the time vector.
TM = mFd.Time;
pFd = mFd;
% From the Excel feeding program mFd, depending on the current time value appropriate feed composition is assigned.
for i=2:length(TM),
    if t >= TM(i-1) & t < TM(i) | t > TM(i)
        for j=1:length(pFd.FdNames),
            pFd.(char(pFd.FdNames(j))) = mFd.(char(pFd.FdNames(j)))(i);
        end
    end
end
% Assigning the Excel programmed feed to actual feed.
Fd = pFd;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specific code for Cation addition by pH controller.
%     if Conn==1,
%         Q_add = u(1);
%         Scat_add = u(2);
%         % Update feed flow.
%         Fd.Q_inf = pFd.Q_inf + Q_add;
%         % Update all non gas influent concentrations, correcting Scat afterwards.
%         for j=(1+2):length(pFd.FdNames)-5,
%             Fd.(char(Fd.FdNames(j))) = pFd.(char(Fd.FdNames(j))) * (pFd.Q_inf/Fd.Q_inf);
%         end
%         Fd.Scat_inf = (pFd.Q_inf * pFd.Scat_inf + Q_add*Scat_add)/Fd.Q_inf;
%         % In case of zero flow batch process.
%         if pFd.Q_inf==0,
%             Fd = pFd;
%             Fd.Q_inf = Q_add;
%             Fd.Scat_inf = Scat_add;
%         end
%     end
    % END code for Cation addition by pH controller.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The feeding flows and concentrations are copied in a vector.
FdV = zeros(length(Fd.FdNames),1);
for i=1:length(FdV),
    FdV(i) = Fd.(char(Fd.FdNames(i)));
end
Fd.FdV = FdV;


% % If the feeding comes from other reactor.
% if Conn==1,
%     Fd=mFd;
%     for j=1:length(pFd.FdNames),
%         Fd.(char(Fd.FdNames(j))) = u(j);
%     end
% else
%     % Feed flow and concentrations are given in a matrix as well as the time vector.
%     % To avoid errors the first element of the vector of times is put to zero independently of the given value.
%     TM = pFd.Time;
%     Fd = pFd;
%     % To make easier the assignation, an extra element is added to TM vector with a very high value.
%     for i=2:length(TM),
%         if t >= TM(i-1) & t < TM(i) | t > TM(i)
%             for j=1:length(Fd.FdNames),
%                 Fd.(char(Fd.FdNames(j))) = pFd.(char(Fd.FdNames(j)))(i);
%             end
%         end
%     end
% end
% 
% % Feeding concentrations in a vector.
% FdV = zeros(R.Dim.numSt,1);
% for i=1:length(FdV),
%     FdV(i) = Fd.(char(Fd.FdNames(i+(length(Fd.FdNames)-R.Dim.numSt))));
% end    
% 
% % Output to the global R.
% R(idR).Fd = Fd;
% R(idR).FdV = FdV;