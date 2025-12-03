function varargout = addgrid(LineWidth, Alpha, LineStyle)
% function varargout = addgrid(GridLayer,Alpha,LineStyle)
% return only grid axis handle

box on; grid on;

% LineWidth0 = get(gca,'LineWidth');
LW0 = 4/5; % set(hg,'LineWidth') = 4/3;

Alpha0 = 1/5; LS0 = ':';
SetDefaultValue(1,'LineWidth' , LW0)
SetDefaultValue(2,'Alpha'     , Alpha0)
SetDefaultValue(3,'LineStyle' , LS0)
% SetDefaultValue(1,'GridLayer' , 0)

set(gca,'LineWidth', LineWidth)

hg = gca;
hg.GridLineStyle  = LineStyle;
% check if property exists for older versions of matlab
% if isprop(hg,"GridLineWidth"); hg.GridLineWidth = LineWidth; end
hg.GridAlpha = Alpha;
% hg.Layer     = 'bottom';
% if GridLayer == 1; hg.Layer = 'top'; end

for k = 1:nargout
  varargout{k} = hg;
end

% SHRINK TICKS IF NEEDED --> better to do manually
% tickshrink(.5)

% AXIS LAYER ON TOP OR BOTTOM --> better to leave under the plots, looks ugly otherwise
% layer axis on tp  
% hg.Layer = 'top';