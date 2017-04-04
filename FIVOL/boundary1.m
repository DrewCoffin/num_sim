function [phi] = boundary1(phi,phix,nth,nr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
phi(:,1) = 0
phi(nr,:) = 2.0;
phi(1,:) = 1.0;
phi(:,nth) = 0.5;
end

