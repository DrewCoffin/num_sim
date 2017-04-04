function [phi] = boundary2(phi,nth,nr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
phi(:,1) = 0;
phi(nr,:) = 1.0;
phi(1,:) = 0.0;
phi(:,nth) = 0.0;
phi(nr/2 - 2:nr/2+2,nth/2 -2:nth/2 +2) = 2.0;
end

