function [Me]=elm(le,m,J0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble element mass matrix
% File name: elM.m
% 
% le [m]	Element length
% m  [kg]	Element mass per unit length
% Me is returned - element mass matrix
%
% Make sure the mass matrix is symmetric!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Me =  eye(6) * [le*m/2; (m*le^3)/24; J0 * le/2; le*m/2; (m*le^3)/24; J0 * le/2];
    
end

