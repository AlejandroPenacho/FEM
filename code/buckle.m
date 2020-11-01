function [pb,ub]=buckle(Ks,Ksigmas,nnode,node_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve beam buckling equation
% File name: buckle.m
% 
% Ks		Structural stiffness matrix
% Ksigmas	Structural inital stiffness matrix
% nnode		Number of nodes
% node_z	Nodal x-coordinates
%
% pb		Matrix with eigenvalues (load factors) in diagonal
% ub		Corresponding matrix of eigenvectors (buckling modes)
% 	(Column i of ub shows the buckling mode for buckling load (i,i) in pb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate eigenvalues and eigenvectors

% Split bending and twist modes into separate vectors

% Normalise deflections, rotations and twist for plotting purposes (without risking to mix up signs or divide by zero)

% Plot buckling modes

% Present the buckling loads






