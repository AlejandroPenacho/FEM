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

Ke = zeros(6,6);

%Making the element stiffness matrix
Ke(1,1:6) = (6*EI/le^2)*[2/le, -1, 0, -2/le, -1, 0];
Ke(2,2:6) = (2*EI/le)*[2, 0, 3/le, 1, 0];
Ke(3,3:6) = (GJ/le)*[1, 0, 0, -1];
Ke(4,4:6) = (6*EI/le^2)*[2/le, 1, 0];
Ke(5,5:6) = (2*EI/le)*[2, 0];
Ke(2,6)   = GJ/le;

%Making it symmetric
Ke = (Ke+Ke') - eye(size(Ke,1)).*diag(Ke);
















