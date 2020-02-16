% 3D Plot function example syntax:
%  plotLake3D(tout, allOut,{'Sno3'});
function [] = plotLake3D(tout, allOut, varName)
global R;
varNames = [R.St.StNames;R.AlgSt.AlgStNames];
% Find id of the name chosen.
id = find(strcmp(varName,varNames));

% Building 3D matrix for 3D plots.
fprintf('>> BUILDING 3D MATRIX WITH ALL RESULTS AND 3D PLOTTING, PLEASE WAIT...\n');
R.tOutM = tout;

for i=1:size(allOut,1),
    allOutMt = zeros(R.Dim.numGrd, R.Dim.NOutG);
    for j=1:R.Dim.numGrd,
        allOutMt(j,:) = [allOut(i,(j-1)*R.Dim.NOutG+1:j*R.Dim.NOutG)];
    end
    R.simOutM(i,:,:) = allOutMt;
end

[T, zL]=meshgrid(R.tOutM, R.pPhys.z);
%surf(T,zL,R.simOutM(:,:,id)');
mesh(T,zL,R.simOutM(:,:,id)');
axis([0 tout(length(tout)) R.pPhys.z(1) R.pPhys.z(length(R.pPhys.z)) 0 max(max(R.simOutM(:,:,id)))]);
if isequal(varName,{'pH'})==1,
    axis([0 tout(length(tout)) R.pPhys.z(1) R.pPhys.z(length(R.pPhys.z)) min(min(R.simOutM(:,:,id))) max(max(R.simOutM(:,:,id)))]);
end
xlabel('Time (h)');
ylabel('Depth (m)');
zlabel(varName);
fprintf('....DONE!!\n');