function varargout = addlegend(LH,legNames,Anchor,FntSize,INTerPret,LineLength,Buffer,LineThickness,NumCols)%
% --------------------------------------------------------------------------------------------------
% CALL AS:  addlegend(LH,legNames, ... optional ) or
%           addlegend([],legNames, ... optional ) MUST INCLUDE []
% --------------------------------------------------------------------------------------------------
% QUICKLY ADD LEGEND TO PLOT WITH STANDARD SETTINGS

% to stop writing Latex font warnings to screen
% [lastMsg, lastId] = lastwarn;  % Get the last warning message and its ID; warning('off', lastId); 
% warning off

warning('off', 'MATLAB:legend:IgnoringExtraEntries') 
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')

% if ~isa(LH, 'matlab.graphics.chart.primitive.Line')
%   legNames = LH;
%   LH = [];
% end
% 
% ~(isa(LH, 'matlab.graphics.chart.primitive.Line') || isa(LH,'double'))
%   legNames = LH;
%   LH = [];

% set some default values
FNS0 = get(gca,'FontSize');
SetDefaultValue(3,'Anchor'        ,  3     );     
SetDefaultValue(4,'FntSize'       ,  FNS0);
SetDefaultValue(5,'INTerPret'     ,  'LaTex'  );
SetDefaultValue(6,'LineLength'    ,  4/5 );  
SetDefaultValue(7,'Buffer'        ,  [0 0]  );
SetDefaultValue(8,'LineThickness' ,  1.5 );
SetDefaultValue(9,'NumCols'       ,  1  );

% convert to legnames to string if not already
if ~iscell(legNames); legNames = cellstr(legNames); end

if isempty(LH)
  h_leg =	legendflex( legNames, ...
  						        'fontsize'	  , FntSize, ...
  						        'anchor'		  , Anchor*ones(1,2), ... % 3=top Right, 2=top Center, 1=top Left, 5=bottom Right, 7=bottom Left, 
  						        'buffer'		  ,	Buffer, ... 
                      'xscale'      , LineLength, ...      % length of the lines
                      'ncol'        , NumCols   , ...
  						        'Interpreter' , INTerPret);
else
  h_leg =	legendflex( LH, ...
                      legNames, ...
  						        'fontsize'	  , FntSize, ...
  						        'anchor'		  , Anchor*ones(1,2), ... % 3=top Right, 2=top Center, 1=top Left, 5=bottom Right, 7=bottom Left, 
  						        'buffer'		  ,	Buffer, ... 
                      'xscale'      , LineLength, ...      % length of the lines
                      'ncol'        , NumCols   , ...
  						        'Interpreter' , INTerPret);
end

% --------------------------------------------------------------------------------------------------
% set the linewidth of the legend lines
% LineThickness = 1.5;
legend_lines = findobj(h_leg, 'type', 'line');
set(legend_lines, 'LineWidth', LineThickness);

% varible output arguments
if nargout > 0
	varargout{1} = h_leg;
end

warning('on', 'MATLAB:legend:IgnoringExtraEntries') 








end










































