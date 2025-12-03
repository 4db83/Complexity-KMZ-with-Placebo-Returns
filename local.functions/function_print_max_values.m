function [p_kmz, z_kmz, max_kmz, p_0, z_0, max_0] = function_print_max_values(tbl_kmz, tbl_0, PTYPE, Plist, trnwin, PRINT_COMP, MIN, decml)
% helper function to print the max values and generate p,z max output for kmz and 0.

SetDefaultValue(6,'PRINT_COMP',0)
SetDefaultValue(7,'MIN'       ,0);  % default is to compute the max value
SetDefaultValue(8,'decml'     ,6)

if MIN == 1
  tbl_kmz = -1.*tbl_kmz;
  tbl_0   = -1.*tbl_0;
end

[p_kmz, z_kmz] = find(max(max(tbl_kmz.Variables)) == tbl_kmz.Variables);
[p_0,   z_0]   = find(max(max(tbl_0.Variables))   == tbl_0.Variables);
max_kmz = max(max(tbl_kmz.Variables));
max_0   = max(max(tbl_0.Variables)); sep('=')

length_Pz_kmz = length( char(tbl_kmz.Row(p_kmz))) + length( char(tbl_kmz.Properties.VariableNames(z_kmz)) );
length_Pz_0   = length( char(tbl_0.Row(p_0)))     + length( char(tbl_0.Properties.VariableNames(z_0)) );   

% length(char(PTYPE))
ADJST = "\t";
if length(char(PTYPE)) < 5;  ADJST = "\t\t"   ; end
if length(char(PTYPE)) < 3;  ADJST = "\t\t\t" ; end

MIN_or_MAX = "Max ";
if MIN == 1
  MIN_or_MAX = "Min ";
end

if MIN ==1
  max_kmz = -max_kmz;
  max_0   = -max_0;
end

fprintf(MIN_or_MAX + PTYPE + "(KMZ)" + ADJST + " = % 4.6f ",            max_kmz  )
fprintf(' RFF(%s,',             char(tbl_kmz.Row(                     p_kmz)))
fprintf(' %s). '  ,             char(tbl_kmz.Properties.VariableNames(z_kmz))); 

if length_Pz_kmz < 10
  fprintf('\t\t\t Complexity c = %.4g.\n',  Plist(p_kmz)/trnwin);
elseif length_Pz_kmz > 12
  fprintf('\t Complexity c = %.4g.\n',      Plist(p_kmz)/trnwin);
else 
  fprintf('\t\t Complexity c = %.4g.\n',    Plist(p_kmz)/trnwin);
end

fprintf(MIN_or_MAX + PTYPE + "(0)  "   + ADJST + " = % 4.6f ",           max_0    )
fprintf(' RFF(%s,',             char(tbl_0.Row(                       p_0))  )
fprintf(' %s). '  ,             char(tbl_0.Properties.VariableNames(  z_0))  ); 
if length_Pz_0 < 9
  fprintf('\t\t\t Complexity c = %.4g.\n',  Plist(p_0)/trnwin);
elseif length_Pz_0 > 12
  fprintf('\t Complexity c = %.4g.\n',      Plist(p_0)/trnwin);
else
  fprintf('\t\t Complexity c = %.4g.\n',    Plist(p_0)/trnwin);
end

% Print this to double check output
if PRINT_COMP == 1
  sep
  print_table(sortrows(table_max(tbl_kmz, PTYPE + "(KMZ)") , PTYPE + "(KMZ)",'descend'), decml, 1,0,2);sep
  print_table(sortrows(table_max(tbl_0,   PTYPE + "(0)"  ) , PTYPE + "(0)"  ,'descend'), decml, 1,0,2);sep;sep('')
end

end


