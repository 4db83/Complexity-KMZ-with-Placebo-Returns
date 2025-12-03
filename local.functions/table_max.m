function tbl_out = table_max(table_in,NewVarName)
% returns the max
% first turn off warnings
warning('off', 'MATLAB:table:ModifiedVarnamesRows2Vars'); 

SR_m = table_in;
SR_m2 = rows2vars(SR_m); SR_m2(:,'OriginalVariableNames') = []; 
SR_m2.Properties.RowNames = table_in.Properties.VariableNames; 
max_SR_m2 = rows2vars(max(SR_m2)); max_SR_m2.Properties.RowNames = max_SR_m2.OriginalVariableNames;
max_SR_m2(:,'OriginalVariableNames') = [];
max_SR_m2 = renamevars(max_SR_m2,'Var1',NewVarName);
tbl_out = max_SR_m2;