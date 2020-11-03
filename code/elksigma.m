function [Kesigma]=elksigma(le,P,I0,A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble element initial stress stiffness matrix
% File name: elksigma.m
% 
% le [m]	Element length
% P  [N]	"Tensile" buckling load
% Kesigma is returned - element initial stress matrix
%
% Make sure the initial stress matrix is symmetric!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kesigma = zeros(6,6);

%Making the element intitial stress matrix
%All elements apart from row 3 and 6 taken from equation 48. Rows 3 and 6
%taken from equation 77
Kesigma(1,1:6) = (P/(30*le))*[36, -3*le, 0, -36, -3*le, 0];
Kesigma(2,2:6) = (P/(30*le))*[4*le^2, 0, 3*le, -le^2, 0];
Kesigma(3,3:6) = (I0*P/(A*le))*[1, 0, 0, -1];
Kesigma(4,4:6) = (P/(30*le))*[36, 3*le, 0];
Kesigma(5,5:6) = (P/(30*le))*[4*le^2, 0];
Kesigma(6,6)   = (I0*P/(A*le));

%Making it symmetric
Kesigma = (Kesigma+Kesigma') - eye(size(Kesigma,1)).*diag(Kesigma);