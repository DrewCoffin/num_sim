function [x,y,TH,R] = get_grid(th,r)
%th is the theta coordinate
%r is the radial coordinate
[TH,R] = meshgrid(th,r); %establish grid
[x,y] = pol2cart(TH,R);  %convert to x-y coordinates
end

