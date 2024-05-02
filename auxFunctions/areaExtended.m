function [hP, hN, ahP] = areaExtended(x, varargin) 

Y = varargin{1} ;

% Positive components. 
fh(1) = figure(get(gcf,'Number')); 
hP    = area(x,Y.*(Y>0), 'Facecolor', 'flat', 'Facealpha', 1, 'linestyle', 'none', 'EdgeColor', 'none', 'linewidth', .01) ; 
ahP   = gca; 

% Negative components. 
fh(2) = figure; 
hN    = area(x, Y.*(Y<0), 'Facecolor', 'flat', 'Facealpha', 1, 'linestyle', 'none', 'EdgeColor', 'none', 'linewidth', .01) ; 


end 