function [defl,teta,fi,umax,tetamax,fimax,RF]=bending(Ks,Qs,K,Q,nnode,node_z,S,q,T,E,I,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate deformed beam bending and torsion shape and plot results
% File name: bending.m
%
% Ks		Structural stiffness matrix
% Qs		Structural load vector
% K		System stiffness matrix
% Q		System load vector
% nnode		Number of nodes
% node_z	Nodal z-coordinates
%
% defl		Deflection vector of size nnodes
% teta		Rotation vector of size nnodes
% fi		Twist vector of size nnodes
% umax		Maximum deflection
% tetamax	Maximum rotation
% fimax		Maximum twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve system of equations
uReducedVector = linsolve(Ks,Qs);

% Present displacements at the free end

umax    = uReducedVector(end-2);
tetamax = uReducedVector(end-1);
fimax   = uReducedVector(end);

% Present reaction forces

RF = K * [0;0;0; uReducedVector] - Q;

% Create result vector containing deflections, rotations and twist

uVector = [0;0;0; uReducedVector];

% Split deflections, rotations and twist into separate vectors
defl = uVector(1:3:end);
teta = uVector(2:3:end);
fi   = uVector(3:3:end);

% Normalise deflections, rotations and twist and plot results

normalizedDefl = defl / umax;
normalizedTeta = teta / tetamax;
normalizedFi = fi / fimax;

%Calculating analytical results to compare values


x = [1:-1/(length(node_z)-1):0];
for i = 1:length(x)
    %Deflection
    if x(i) > 0
        defl_S(i) = S*L^3/(6*E*I) * (x(i)^3-3*x(i)+2);
        defl_T(i) = T*L^2/(2*E*I) * (x(i)^2 - 2*x(i)+1);
        
        
        %Teta
        teta_S(i) = S*L^2/(2*E*I) * (1-x(i)^2);
        teta_T(i) = T*L/(E*I) * (1-x(i));
       
        
    else
        defl_S(i) = S*L^3/(6*E*I) * (-1 + 3*(1-x(i)));
        defl_T(i) = T*L^2/(2*E*I) * 2*(1-x(i)-0.5);
        teta_S(i) = S*L^2/(2*E*I);
        teta_T(i) = T*L/(E*I);
    end
    defl_q(i) = q*L^3/(24*E*I) * (x(i)^4-4*x(i)+3);
    teta_q(i) = q*L^2/(6*E*I) * (1-x(i)^3);
    
end

defl_analytical  = defl_S + defl_T + defl_q;
teta_analytical = teta_S + teta_T + teta_q;
normalizedDef_an = defl_analytical/max(defl_analytical);
normalizedTeta_an = teta_analytical/max(teta_analytical);


figure (1)

subplot(3,1,1)
title("Deflection")
scatter(node_z, defl)
xlabel("z")
ylabel("u/u_{max}")
grid minor

subplot(3,1,2)
title("Teta")
scatter(node_z, teta)
xlabel("z")
ylabel("\theta/\theta_{max}")
grid minor

subplot(3,1,3)
title("Fi")
scatter(node_z, normalizedFi)
xlabel("z")
ylabel("\phi/\phi_{max}")
grid minor

figure (2)

subplot(2,1,1)
title("Deflection")
scatter(node_z, defl_analytical)
xlabel("z")
ylabel("u/u_{max}")
grid minor

subplot(2,1,2)
title("Teta")
scatter(node_z, teta_analytical)
xlabel("z")
ylabel("\theta/\theta_{max}")
grid minor

end
