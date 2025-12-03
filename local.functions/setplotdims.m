function varargout = setplotdims(figdim,axLineWidth) 
%{F: specifies plots dimensions, fontsize and fontname
%===============================================================================
% Makes it easy to view strings in excel which are actually numbers without 
% converting them to numbers`
%-------------------------------------------------------------------------------
% 	USGAGE:	(1) setplot([height], Fontsize)
%           (2) setplot([top-positon height], Fontsize)
%						(3) setplot([top-positon width height], Fontsize) (rightpos at 0)
%						(4) setplot([leftpos rightpos height width], Fontsize) 
%           [leftpos rightpos] is optional and will only be called if one 
%           needs to center the plot. Default font size is 10, fontname is Palation.
%-------------------------------------------------------------------------------
% 	INPUT : 
%	  figdim		  =  (1x2) matrix with [height width] dimension of plot
%   axLineWidth	=  scalar, linewidth.
% 	OUTPUT:       
%	  zero arguments: adjusted plot.
% OLD:
% function varargout = setplot(figdim,fntsize,y_digits,axLineWidth,fighandle) 
%===============================================================================
% 	NOTES :   Always call it at the end of the plotting commands, after hold off;
%-------------------------------------------------------------------------------
% Created :		27.07.2013.
% Modified:		03.11.2023.
% Copyleft:		Daniel Buncic.
%------------------------------------------------------------------------------%}

widh_def = .9;
% SetDefaultValue(1, 'figdim'			, [widh_def .2]);
SetDefaultValue(1, 'figdim',      get(gca,'Position'));
SetDefaultValue(2, 'axLineWidth', 5/5);

xd = size(figdim);
if ~(max(xd) == 1)
	Big1 = figdim > 1;
	if any(Big1) > 0
		figdim = figdim./10;
	end
end

% added the font name conversion now to print in TNR as opposed to Palation!
FName = 'Times New Roman';
top_	= .45;
left_ = .05;
ax    = gca;

% if nargin < 2
  % if only one input argument is given
  if max(xd) == 2
    % sets width and height, takes the centered location [.1 .2]
%     set(fighandle,'Position',[[left_ top_] figdim],'FontSize',fntsize,'FontName',FName);
    set(ax,'Position',[left_ figdim(1) widh_def figdim(2)],'FontName',FName);
  elseif max(xd) == 3
    % sets width and height as well as the location
    set(ax,'Position',[left_ figdim(1:3)],'FontName',FName);
	elseif max(xd) == 4
    % sets width and height as well as the location
    set(ax,'Position',figdim,'FontName',FName);
  elseif max(xd) == 1
    % only scalar input means that only the font size is to be set.
%     fntsize = figdim;
%     set(fighandle,'FontSize',fntsize,'FontName',FName);
		set(ax,'Position',[[left_ top_ widh_def] figdim],'FontName',FName);
  else 
    disp('Error in plot diminsion');
	end

% MAKE LINEWIDTH OF BOX
ax.LineWidth = axLineWidth;

% plotting standards, can be removed
% tickshrink(2/3);
box on; grid on;
set(gca,'GridLineStyle',':','GridAlpha',1/3)
plot_position = get(gca,'Position');
% set(gca, 'Layer', 'bottom')

if nargout > 0
  varargout{1} = plot_position(2);
end


