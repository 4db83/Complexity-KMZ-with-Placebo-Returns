function plot2pdf(file_name, PaperPosition, KeepEPS)

SetDefaultValue(2,'PaperPosition' ,'auto')
SetDefaultValue(3,'KeepEPS'       , 0)

set(gcf,'PaperPositionMode', PaperPosition);

% file_name = strcat( '../graphics/', 'YY=11', num2str(MA6));
% check if file_name includes eps extension

if length(file_name) < 5
 I_eps = 0;
else
  I_eps = strcmp(file_name(end-3:end),'.eps');
end

if I_eps
  file_name = file_name;
else
  file_name = strcat(file_name,'.eps');
end

% print 2 eps first
print(gcf, '-depsc2', file_name );
% check if epstopdf exists, if so, convert to PDF and delete eps
status = system(['epstopdf',' ',[file_name]]);
if status
  disp('    NO epstopdf on your system: --> file not converted to PDF. Work with EPS')
else
  if ~KeepEPS
    delete(file_name)
  end
end

end

% delete eps