%Script to calculate Gibbs as a function of temperature, acid-base
%constants matrix (KaM) and Dehydration constants (KdM)

function [KaM, KdM, G0ftM] = DG_Ka_calc(pTh, pOp, T, idL)

%Temperature modified Gibbs free energies of formation according to Van't Hoff's equation.
G0ftM = pTh.H0fM + (T/298.15)* (pTh.G0fM -pTh.H0fM);

%Replace dehydration of water to 0
G0ftM(length(G0ftM(:,1)),1)=0;

%Pre-allocate KaM
KaM = zeros(length(G0ftM(:,1)),length(G0ftM(1,:))-1);

%Calculate Ka Matrix
for i = 1:length(G0ftM(1,:))-2
    KaM(:,i) = ((G0ftM(:,i+2)~=0)&(G0ftM(:,i+1)~=0)).* exp(-(G0ftM(:,i+2)-G0ftM(:,i+1))./(pOp.Rth * T));
end

%Calculate Kd Matrix (actually is a vector)
KdM = (G0ftM(:,1)~=0) .* exp (-(G0ftM(:,2)-G0ftM(:,1)-G0ftM(length(G0ftM(:,1)),2))./(pOp.Rth * T));

% %Reshape to include Ka and Kd for aqueous-phase components only
% idL = (strcmp(St.Phase,'L'));
%Selection only the state variables that are on liquid phase
KaM = [KaM(idL==1,:); KaM(length(KaM(:,1)),:)];
KdM = [KdM(idL==1,:); KdM(length(KdM(:,1)),:)];

    





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