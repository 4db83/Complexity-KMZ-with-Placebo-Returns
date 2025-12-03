function [] = addrowspace(FNS)
% function add row-space between figures when using nexttile by increasing the font of a last figure entry
% simply call as: addrowspace or addrowspace(40) or so, where 40 is the fontsize.  
% may need to be careful to contain enough rows in the tiledlayout definition to have enough columns.
% 
%
% NOTE: this does not work with setoutside ticks, because this creates another axis and only one
% axis gets pushed up.


if nargin < 1 || isempty(nargin)
  FNS = 20;
end

nexttile
setyticklabels([],0,FNS)
set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide only the axes lines
axis off;