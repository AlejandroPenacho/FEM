function [Ke]=elk(le,EI,GJ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble element stiffness matrix
% File name: elk.m
% 
% le [m]	Element length
% EI [Nm2]	Element bending stiffness (constant in the lab)
% GJ [Nm2]	Element torsional stiffness (constant in the lab)
%
% Ke is returned - element stiffness matrix
%
% Make sure the stiffness matrix is symmetric!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

