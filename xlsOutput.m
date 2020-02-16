% Jorge Rodríguez
% Masdar Institute
% PO Box 54224
% Abu Dhabi
% United Arab Emirates
% jrodriguez@masdar.ac.ae

% Collecting the output results from the MFC model including the biofilm and writting in the 
% simResults.xls file from Excel.
global R;
% Copying template as empty simResults file.
% WARNING THE PREVIOUYS simResults FILE WILL BE OVERWRITEN.
copyfile('temp simResults.xls', 'simResults.xls');
fprintf('>> NOTE THAT THE PREVIOUS simResults.xls FILE (IF ANY) HAS BEEN DELETED!!...\n');
fprintf('>> SAVING FINAL SIMULATION RESULTS INTO simResults.xls FILE, PLEASE WAIT...\n');
nCell='A2';
% Writting the variables 
nSheet = 'OutputFinal';

% Build matrix with results and layers.
idFinal = size(allOut,1);
for j=1:R.Dim.numGrd,
    allOutMt(j,:) = [allOut(idFinal,(j-1)*R.Dim.NOutG+1:j*R.Dim.NOutG)];
end
% Adding depth coordinates to first columns.
dta = [R.pPhys.z' R.StM' R.AlgStM' R.rRcM' R.rTrM' ];
eval(strcat('xlswrite(', char(39),'simResults.xls', char(39), ', [dta], nSheet, nCell);'));
fprintf('... DATA SAVED...\n');

% Building 3D matrix for 3D plots.
fprintf('>> BUILDING 3D MATRIX WITH ALL RESULTS, PLEASE WAIT...\n');
R.tOutM = tout;
R.simOutM = zeros(length(R.tOutM), R.Dim.numGrd, R.Dim.NOutG);
for i=1:size(allOut,1),
    allOutMt = zeros(R.Dim.numGrd, R.Dim.NOutG);
    for j=1:R.Dim.numGrd,
        allOutMt(j,:) = [allOut(i,(j-1)*R.Dim.NOutG+1:j*R.Dim.NOutG)];
    end
    R.simOutM(i,:,:) = allOutMt;
end
save simResult.mat
fprintf('... 3D MATRIX BUILT AND SAVED IN simResults.mat!!\n');