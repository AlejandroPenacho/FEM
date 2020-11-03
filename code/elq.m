function [Qe]=elq(le,q,qt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble element load vector
% File name: elq.m
% 
% le [m]	Element length
% q  [N/m]	Element distributed load (constant in the lab)
% qt  [N]	Element distributed torque (constant in the lab)
% Qe is returned - element load vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qe = (le/2)*[q, -q*le/6, qt, q, q*le/6, qt]; %Element load vector using the equation from page 38 as well as eq. 79
