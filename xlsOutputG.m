% SAVES THE SIMULATION RESULTS IN THE EXCEL FILE 'simResults.xls',
% created from copy of 'temp simResults.xls'. 
% Any existing 'simResults.xls' will be overwriten.
global R;
% Copying template as empty simResults file.
% THE PREVIOUS simResults FILE (IF ANY) WILL BE OVERWRITEN.
% File='C:\Users\mpaton\Documents\MATLAB\ADM1 Study\myTelemac\simID\simResults.xls';
FileName = 'simResultsG.xlsx';
File = strcat(pwd,'\', FileName);
fprintf('>> NOTE THAT THE PREVIOUS simResults.xlsx FILE (IF ANY) HAS BEEN DELETED!!...\n');
fprintf('>> Collecting results into Excel file, please wait...\n');
copyfile('temp simResultsG.xlsx', File);

%Activates Excel server (I needed to add the file location with its full path)
Excel = actxserver ('Excel.Application');
if ~exist(File,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(File,1);
    ExcelWorkbook.Close(false);
end
invoke(Excel.Workbooks,'Open',File);

%Convert allOut into a Grid 3D Matrix Results
allOutG = zeros(length(tout), R.Dim.NOutG, R.Dim.numGrd);
for i = 1:R.Dim.numGrd, allOutG(:,:,i)=allOut(:,(1+((i-1)*R.Dim.NOutG):i*R.Dim.NOutG));  end

%Writing the output data in the excel
nCellHeader = 'A1';
nCell='A2';
nSheet = cell(1, length(allOutG(1,1,:)));
xlsHeaderNames = ['t (h)'; R.St.StNames; R.AlgSt.AlgStNames; R.rRc.RcNames; R.rTr.TrNames]';


%Writing the names of the output Sheets
for j=1:length(allOutG(1,1,:))
    nSheet{j} = strcat('OutputLayer_', num2str(j-1), 'm');
end

%Writes all at once in every Excel Sheet
for j=1:length(allOutG(1,1,:))
    xlswrite1(FileName, xlsHeaderNames,        nSheet{j}, nCellHeader);
    xlswrite1(FileName, [tout allOutG(:,:,j)], nSheet{j}, nCell);
end
fprintf('... simulation results succesfully stored in the file "simResults.xls" \n');

% %Preallocate variables
% pKtV = zeros(length(R.pKt.pKtNames), 1); 
% pOpV = zeros(length(R.pOp.pOpNames), 1); 
% 
% %Kinetic Parameters update
% for j=1:length(R.pKt.pKtNames),
%     pKtV(j) = R.pKt.(char(R.pKt.pKtNames(j)));
% end
%     nSheet = 'KinetParam';   nCellPar = 'B1';
%     xlswrite1(FileName, pKtV, nSheet, nCellPar);
%         
% %Operational Parameters update
% for j=1:length(R.pOp.pOpNames),
%     pOpV(j) = R.pOp.(char(R.pOp.pOpNames(j)));
% end
%     nSheet = 'OperatParam';   nCellPar = 'B1';
%     xlswrite1(FileName, pOpV, nSheet, nCellPar);
%     clearvars pOpV pKtV nCellPar nCell j File FileName nSheet
%     
% fprintf('... Operational and Kinetic parameters succesfully stored in the file "simResults.xls" >>\n');

invoke(Excel.ActiveWorkbook,'Save');
Excel.Quit
Excel.delete 
clear Excel