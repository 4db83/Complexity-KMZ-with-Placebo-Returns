function show(data_in, decimals, leftSpace)
SetDefaultValue(2,'decimals'  , 4)
SetDefaultValue(3,'leftSpace' , 6)

if iscell(data_in)
  [nR, ~] = size(data_in);
  indx = (1:nR)';
   for k = 1:nR
    fprintf( '%6.6g ', indx(k,1));                  % First column: 0 decimal places
    fprintf( ' % s', char(data_in(k,1)));           % First column: 0 decimal places
    fprintf('\n');
   end 
else
  [nR, nC] = size(data_in);
  if nR == 1; data_in = data_in'; end

  if nargin > 2; leftSpace0 = leftSpace; end
  % first column space fix this
  first_column_space  = 6;
  leftSpaceInteger    = 7;

  kk = 4;

  if decimals > 2;  leftSpace = kk+7;  end
  if decimals > 3;  leftSpace = kk+8;  end
  if decimals > 4;  leftSpace = kk+10; end
  if decimals > 5;  leftSpace = kk+11; end
  if decimals > 6;  leftSpace = kk+12; end
  if decimals > 7;  leftSpace = kk+13; end
  if decimals > 8;  leftSpace = kk+14; end
  if decimals > 9;  leftSpace = kk+16; end
  if decimals > 10; leftSpace = kk+18; end
  if decimals > 14; leftSpace = kk+20; end

  % use passed through value if nargin > 2
  if nargin > 2; leftSpace = leftSpace0; end

  [Nr,Nc] = size(data_in);
  indx = (1:Nr)';
  indxS = num2str(indx);

  if iscell(data_in)
    aout = [cellstr(indxS) data_in];
  elseif istable(data_in)
    if ~isempty(data_in.Row)
      row_names = data_in.Row;
      aout = [cellstr(indxS) row_names cellstr(num2str(data_in.Variables))];
    else
      aout = [cellstr(indxS) cellstr(num2str(data_in.Variables))];
    end
  else
    aout = [data_in];
    % aout = [cellstr(indxS) cellstr(num2str(data_in,'%2.4f'))];
  end

  % find integer columns of data_in, if any
  isIntegerColumn = all(mod(aout,1) == 0, 1);
  bigDecimalsL     = repmat(leftSpace,Nc,1).*~isIntegerColumn';
  bigDecimalsL(bigDecimalsL==0) = leftSpaceInteger;
  bigDecimalsR     = repmat(decimals,Nc,1) .*~isIntegerColumn';

  % ['%' num2str(leftSpace) '.' num2str(bigDecimals(1)) 'f']

  % IsG = strcmp(decimals,'%9.9g');
  IsG = (decimals == 0);

  for k = 1:Nr
    fprintf( strcat(' %', num2str(first_column_space), '.0f'), indx(k,1));           % First column: 0 decimal places
    for cc = 1:Nc
      % strcat('%11.',num2str(bigDecimals(cc)),'f')
      % fprintf( [' %' num2str(leftSpace) '.' num2str(bigDecimalsR(cc)) 'f'], aout(k,cc) );       % Remaining columns: 4 decimal places
      % fprintf( [' %' num2str(bigDecimalsL(cc)) '.' num2str(bigDecimalsR(cc)) 'f'], aout(k,cc) );       % Remaining columns: 4 decimal places
      if IsG==1
        fprintf(' % 14.14g    ', aout(k,cc) );       % Remaining columns: 4 decimal places
      else
        fprintf( [' %' num2str(bigDecimalsL(cc)) '.' num2str(bigDecimalsR(cc)) 'f'], aout(k,cc) );       % Remaining columns: 4 decimal places
      end
    end
    fprintf('\n');                          % New line
  end
  % return nothing, just print to screen
  % decimals
  % if decimals == 0
  %   disp(aout)
  % else
  %   lst(aout,decimals)
  % end
  % info.fmt = {'%4.0f','%12.2f'}
  % % myprint0(aout,info)
  % myprint(aout)
  % for k = 1:length(aout)
  %   fprintf('% d ', aout(k,1) )
  %   lst( aout(k,2:end), decimals)
  % end

end

function info = myprint0(y,info)
% PURPOSE: print an (nobs x nvar) matrix in formatted form
% CALL AS: myprint(y,info) or myprint(y)
% ---------------------------------------------------------------------------------------------
% USAGE:     mprint(x,info) 
% where: x         = (nobs x nvar) matrix (or vector) to be printed
%        info      = a structure containing printing options
%        info.begr = beginning row to print,    (default = 1)
%        info.endr = ending row to print,       (default = nobs)
%        info.begc = beginning column to print, (default = 1
%        info.endc = ending column to print,    (default = nvar)        
%				 info.sep  = string, ie., '-','=','+','*' etc. type of sep lines to be printed
%				 info.width  = width of sep lines to be printed
%        info.cnames = an (nvar x 1) string vector of names for columns (optional)
%                      e.g. info.cnames = strvcat('col1','col2');
%                      (default = no column headings)
%        info.rnames = an (nobs+1 x 1) string vector of names for rows (optional)
%                      e.g. info.rnames = strvcat('Rows','row1','row2');
%                      (default = no row labels)
%        info.fmt    = a format string, e.g., '%12.6f' or '%12d' (default = %10.4f)
%                      or an (nvar x 1) string containing formats
%                      e.g., info.fmt=strvcat('%12.6f','%12.2f','%12d'); for nvar = 3
%        info.fid    = file-id for printing results to a file
%                      (defaults to the MATLAB command window)
%                      e.g. fid = fopen('file.out','w'); 
%        info.rflag  = 1 for row #'s printed, 0 for no row #'s (default = 0) 
%        info.width  = # of columns before wrapping occurs (default = 80)                                                  
% ---------------------------------------------------------------------------------------------
% e.g.   info.cnames = {'ColName1','ColName2','ColName3'}; etc
%        info.rnames = {'RowName1','RowName2'}; etc
%				 info.fmt		 = '%4.6f'  (Default is )
%        mprint(y,info), prints entire matrix, column and row headings
%
% NOTES: - defaults are used for info-elements not specified
%        - default wrapping occurs at 80 columns, which varies depending on the
%          format you use, e.g. %10.2f will wrap after 8 columns
% ---------------------------------------------------------------------------------------------
% SEE ALSO: tsprint, mprint_d, lprint
% ---------------------------------------------------------------------------------------------
[ry,cy] = size(y);

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

% header line default length and style
cwidth	= 120;
lsep		= '-';

if nargin < 2
	info.fmt = '%14.6f';
end

% set default row and column names
default_col_names = cellstr([repmat('Col(',cy,1) num2str((1:cy)') repmat(')',cy,1)]);
default_row_names = cellstr([repmat('Row(',ry,1) num2str((1:ry)') repmat(')',ry,1)]);

if ~isfield(info,'cnames')
	info.cnames = default_col_names;
end

if ~isfield(info,'rnames')
	info.rnames = default_row_names;
end

if isempty(info.rnames)
	info.rnames = default_row_names;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ColHeaders_in = info.cnames;
if ischar(ColHeaders_in)
	ColHeaders = ColHeaders_in;
else
	ColHeaders = char(ColHeaders_in);
end
info.cnames = ColHeaders;

RowHeaders_in = info.rnames;
if size(RowHeaders_in,1)==1
	RowHeaders_in = RowHeaders_in';
end

LR = length(RowHeaders_in);

if LR == ry
  % disp("this")
% 	RowHeaders = char(['{Variable}'; RowHeaders_in]);
	% RowHeaders = char([' '; RowHeaders_in]);
  % RowHeaders = char([RowHeaders_in]);
  RowHeaders = RowHeaders_in;
elseif LR == ry+1
	RowHeaders = char([RowHeaders_in]);
else
	 fprintf(' No. of Data Rows (plus header): %d ~= No. of Row Names: %d \n', [ry+1 LR]);
	 error('Wrong # rnames in myprint'); 
end

info.rnames = RowHeaders

% setup defaults
fid = 1; rflag = 0; cflag = 0; rnum = 0; nfmts = 1; 
[nobs nvars] = size(y);
begr = 1; endr = nobs; begc = 1; endc = nvars; 
fmt = '%14.6f';

if nargin == 1
% rely on defaults
elseif nargin == 2
  if ~isstruct(info)
    error('myprint: you must supply the options as a structure variable'); 
  end
fields = fieldnames(info);
nf = length(fields);
for k=1:nf
    if strcmp(fields{k},'fmt')
        fmts = info.fmt; 
  [nfmts junk] = size(fmts);
  if nfmts == nvars
   fmt = fmts;
  elseif nfmts == 1
   fmt = fmts;
  else
   error('myprint: wrong # of formats in string -- need nvar');
  end
    elseif strcmp(fields{k},'fid')
        fid = info.fid;
    elseif strcmp(fields{k},'begc');
        begc = info.begc;
    elseif strcmp(fields{k},'begr');
        begr = info.begr;
    elseif strcmp(fields{k},'endc');
        endc = info.endc;
    elseif strcmp(fields{k},'endr');
        endr = info.endr;
    elseif strcmp(fields{k},'width');
				cwidth = info.width;
		elseif strcmp(fields{k},'sep');
				lsep = info.sep;
    elseif strcmp(fields{k},'cnames');
        cnames = info.cnames;
        cflag = 1;
    elseif strcmp(fields{k},'rnames');
        rnames = info.rnames;
        rflag = 1;
    elseif strcmp(fields{k},'rflag');
        rnum = info.rflag;
    end
end

else
error('Wrong # of arguments to myprint');
   
end % end of if-elseif input checking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	CONTROL LINE STYLE PRINTING UNDER THE HEADERS

sa	= fmt(:,2:end);
saI = strfind(sa,'.');							% find decimal placement
N1_fmts = str2num(sa(1:saI-1));			% before	decimals together
N2_fmts = str2num(sa(saI+1:end-1));	% after		decimals together
dm_addon = 2;
l0_ = repmat('-',1,ry + 2);

if isempty(N1_fmts)
	N1_fmts = 10;
end

if isempty(N2_fmts)
	N2_fmts = 8;
end

if N1_fmts > 10;
	l0_ = repmat('-',1,N1_fmts + 2);
end

T_fmts = N2_fmts + dm_addon + 6;

%LL_ = [l0_ repmat('-',1,T_fmts*cy)];

LL_ = [repmat(lsep,1,cwidth)];

% see if the user supplied row names and set rnum
% correct her mistake if she did this
if rflag == 1
rnum = 0;
end

% parse formats
if nfmts == 1
   f1 = strtok(fmt,'%');
   f2 = strtok(f1,'.'); 
    if strcmp(f1,f2)
     f2 = strtok(f2,'d');
     dflag = 1;
     fflag = 0;
    else
     tmp1			= strtok(fmt,'f');
     [tmp2,a] = strtok(tmp1,'.');
     tmp1 = tmp1(2:length(tmp1));
     tmp2 = tmp2(2:length(tmp2));
     opoint = num2str(str2num(tmp1) - str2num(tmp2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		 % decimal = opoint(1,length(opoint));
		 % decimal = opoint(3:length(opoint));
		 % I HAVE CHANGED THIS BELOW
		 decimal = a(2:length(a));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     f2 = strtok(f2,'f');
     fflag = 1;
     dflag = 0;
    end
   f2 = str2num(f2);
   nwide = floor(cwidth/f2); % 80 columns divided by format
   nvar = endc-begc+1;
   nsets = ceil(nvar/nwide);
else %  wrapping in this case is based on widest format in the list
	
	nwidev		= zeros(nfmts,1);
	nsetsv		= zeros(nfmts,1);
	f2v				= zeros(nfmts,1);
	dflagv		= zeros(nfmts,1);
	fflagv		= zeros(nfmts,1);
	decimalv	= zeros(nfmts,1);

  for ii=1:nfmts
   f1 = strtok(fmt(ii,:),'%');
   f2 = strtok(f1,'.');
    if strcmp(f1,f2)
     f2 = strtok(f2,'d');
     dflagv(ii,1) = 1;
     fflagv(ii,1) = 0;     
    else
     tmp1 = strtok(fmt(ii,:),'f');
     tmp2 = strtok(tmp1,'.');
     tmp1 = tmp1(2:length(tmp1));
     tmp2 = tmp2(2:length(tmp2));
     opoint = num2str(str2num(tmp1) - str2num(tmp2));
     decimalv(ii,1) = opoint(1,length(opoint));     
     f2 = strtok(f2,'f');
     fflagv(ii,1) = 1;
     dflagv(ii,1) = 0;     
    end
   f2v(ii,1) = str2num(f2);
   nwidev(ii,1) = floor(cwidth/f2v(ii,1)); % cwidth columns divided by format
   nvar = endc-begc+1;
   nsetsv(ii,1) = ceil(nvar/nwidev(ii,1));   
	end
nsets = min(nsetsv); 
nwide = max(nwidev);
end 

% if we have row and column labels
% adjust variable labels and column heading strings
% to match the width of the printing format

if rnum == 1
dstr = 'Obs#';
end

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % THIS CONTROLS THE PRINT FROMAT STRING: BEFORE THAT THERE IS SOME FUNCKY TRIMMING AT WORK
%  ffmt = info.fmt;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if cflag == 1 % we have column headings
 [vsize nsize] = size(cnames); % error check cnames argument
 if vsize ~= nvars; 
	 fprintf(' No. of Data Column: %d ~= No. of Column Names: %d \n', [nvars vsize ]);
	 error('Wrong # cnames in myprint'); 
 end    
 if nfmts == 1 % case of only 1 format string
  nmax = max(f2,nsize); % build format strings 
                        % based on widest format              
  sfmt = ['%', num2str(nmax)];
  sfmt = [sfmt,'s ']; 
  ffmt = ['%', num2str(nmax)];
	
   if dflag == 1
			ffmt = [ffmt,'d '];
   elseif fflag == 1
			ffmt = [ffmt,'.'];
			ffmt = [ffmt,decimal];
			ffmt = [ffmt,'f '];
   end
 else % we have multiple format strings, process each
 sfmtv = []; fmtv = [];
  for ii=1:nfmts % find and parse multiple formats
  nmax = max(f2v(ii,:),nsize); % build format strings 
                        % based on widest format              
  sfmtv{ii} = ['%', num2str(nmax)];
  sfmtv{ii} = [sfmtv{ii},'s ']; 
  ffmtv{ii} = ['%', num2str(nmax)];
   if dflagv(ii,1) == 1
		ffmtv{ii} = [ffmtv{ii},'d '];
   elseif fflagv(ii,1) == 1
   ffmtv{ii} = [ffmtv{ii},'.'];
   ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];    
   ffmtv{ii} = [ffmtv{ii},'f '];
   end
  end % end of for ii loop
 end % end of if-else
elseif cflag == 0 % we have no column headings
 if nfmts == 1 % case of only 1 format string
  nmax = f2; % augment format string with a space (the hard way) 
  ffmt = ['%', num2str(nmax)];
   if dflag == 1
   ffmt = [ffmt,'d '];
   elseif fflag == 1
   ffmt = [ffmt,'.'];
   ffmt = [ffmt,decimal];
   ffmt = [ffmt,'f '];
   end
 else % we have multiple format strings, process each
 sfmtv = []; fmtv = [];
  for ii=1:nfmts % find and parse multiple formats
  nmax = f2v(ii,:); % augment format strings with a space 
  ffmtv{ii} = ['%', num2str(nmax)];
   if dflagv(ii,1) == 1
   ffmtv{ii} = [ffmtv{ii},'d '];
   elseif fflagv(ii,1) == 1
   ffmtv{ii} = [ffmtv{ii},'.'];
   ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];    
   ffmtv{ii} = [ffmtv{ii},'f '];
   end
  end % end of for ii loop
 end % end of if-else    
end % end of if-elseif cflag == 0,1

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if rflag == 1 % we have row labels
 [vsize nsize] = size(rnames); % error check cnames argument
 if vsize ~= nobs+1; 
	 fprintf(' No. of Data Rows (plus header): %d ~= No. of Row Names: %d \n', [nobs+1 vsize ]);
	 error('Wrong # rnames in myprint'); 
 end  
 rfmt = ['%', num2str(nsize)]; 
 rfmt = [rfmt,'s']; 
end % end of if rflag == 1

if (rflag == 0 & cflag == 0)
    ffmt = fmt;
end

% print matrix
for j=1:nsets;
 if nfmts == 1 % print row header and column headers
 if rnum == 1;fprintf(fid,'%5s',dstr);     
     elseif rflag == 1    
  fprintf(fid,rfmt,rnames(1,:));
     end  
     if cflag == 1
    for k = (j-1)*nwide+begc:j*nwide+begc-1
  if k <= endc
% find version #; 
  %[version,junk] = version; vers = str2num(version);
   %if vers == 5.2
   fprintf(fid,sfmt,strjust(cnames(k,:),'right'));
% 	 fprintf(fid,sfmt,strjust(cnames(i,:),alignment_));
	 
   %else
   %fprintf(fid,sfmt,strjust(cnames(i,:)));
   %end
  end
 end
     end
  % fprintf(fid,'\n------------------------------------------------------------\n');
  fprintf(fid,['\n' LL_ '\n']);
 else % we have multiple formats
 if rnum == 1;fprintf(fid,'%5s',dstr);     
    elseif rflag == 1   
 fprintf(fid,rfmt,rnames(1,:));
    end
    if cflag == 1
   for k = (j-1)*nwide+begc:j*nwide+begc-1
  if k <= endc
% find version #; 
  %[version,junk] = version; vers = str2num(version);
   %if vers == 5.2
   fprintf(fid,sfmtv{k},strjust(cnames(k,:),'right'));
   %else
   %fprintf(fid,sfmtv{i},strjust(cnames(i,:)));
   %end
  end
   end
    end
 % fprintf(fid,'\n------------------------------------------------------------\n');
 fprintf(fid,['\n' LL_ '\n']);
 end % end of if-else nfmts
 
 % ffmt = info.fmt;
 
 for K = begr:endr % print row labels and numbers in matrix
  if rnum == 1; fprintf(fid,'%5d',K);
        elseif rflag == 1        
  fprintf(fid,rfmt,rnames(K+1,:));
        end
  for L = (j-1)*nwide+begc:j*nwide+begc-1
   if L <= endc
    if nfmts == 1
      % if isnan(y(K,L)) 
      %   % ffmt
      %   fprintf(repmat(' ',17+1,1));
      %   % ffmt
      % else
        fprintf(fid,ffmt,y(K,L));
      % end
    else
    fprintf(fid,ffmtv{L},y(K,L));
    end
   end
  end % end of for l
  fprintf(fid,'\n');
 end % end of for k
%%%%%  fprintf(fid,'\n');			% i have uncommented this to leave no space after last line.
%%%%%  Daniel Buncic.
end % end of for j

end %EOF







end