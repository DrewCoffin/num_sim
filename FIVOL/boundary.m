function [phi] = boundary(phi,phix,nth,nr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
phi(:,1) = 0
phi(nr,:) = phix(nr,:);
phi(1,:) = phix(1,:);
phi(:,nth) = phix(:,nth);
end

