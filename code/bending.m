function [defl,teta,fi,umax,tetamax,fimax,RF]=bending(Ks,Qs,K,Q,nnode,node_z,S,q,T,E,I,L,G,J,qt);

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

umax    = max(abs(uReducedVector(1:3:end)));
tetamax = max(abs(uReducedVector(2:3:end)));
fimax   = max(abs(uReducedVector(3:3:end)));

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


fID = fopen("output/results.txt","a");
fprintf(fID, "Bending results\n\n");
fprintf(fID, "Max. displacement:\t%.3f mm\n", umax*1000 );
fprintf(fID, "Max. bending angle:\t%.3f °\n", tetamax*180/pi );
fprintf(fID, "Max. torsion angle:\t%.3f °\n", fimax*180/pi );

fprintf(fID, "\nForces at clampped end:\n");
fprintf(fID, "\t Shear force: %.3f N\n", RF(1));
fprintf(fID, "\t Moment force: %.3f Nm\n", RF(2));
fprintf(fID, "\t Torque: %.3f Nm\n", RF(3));

fclose(fID);
%Calculating analytical results to compare values


% x = [1:-1/(length(node_z)-1):0];
x = linspace(0,1,300);
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
    defl_q(i) = q*L^4/(24*E*I) * (x(i)^4-4*x(i)+3);
    teta_q(i) = q*L^3/(6*E*I) * (1-x(i)^3);
    phi(i) = L/(G*J)*(T+qt*L)*x(i) - qt*L^2/(2*G*J)*x(i)^2;
end

defl_analytical  = defl_S + defl_q; %defl_T
teta_analytical = teta_S + teta_q; %teta_T
normalizedDef_an = defl_analytical/umax;
normalizedTeta_an = teta_analytical/tetamax;
normalizedPhi_an = phi/fimax;

figure (1)
lineInstances = FEMplot(node_z, uVector, umax, tetamax, fimax);



% subplot(3,1,1)
% title("Deflection")
% scatter(node_z, normalizedDefl)
% xlabel("z")
% ylabel("u/u_{max}")
% grid minor
% 
% subplot(3,1,2)
% title("Teta")
% scatter(node_z, normalizedTeta)
% xlabel("z")
% ylabel("\theta/\theta_{max}")
% grid minor
% 
% subplot(3,1,3)
% title("Fi")
% scatter(node_z, normalizedFi)
% xlabel("z")
% ylabel("\phi/\phi_{max}")
% grid minor

subplot(3,1,1)
title("Deflection")
hold on
h = plot(x*L, flip(normalizedDef_an), '--','LineWidth', 1.5, 'color', [0.4660, 0.6740, 0.1880]);
hold off
legend([lineInstances{1}, h], ["FEM", "Analytical"],'location' ,'northwest')

subplot(3,1,2)
title("Teta")
hold on
h = plot(x*L, flip(normalizedTeta_an), '--','LineWidth', 1.5,'color', [0.4660, 0.6740, 0.1880]);
hold off
legend([lineInstances{2}, h], ["FEM", "Analytical"],'location' ,'northwest')

subplot(3,1,3)
title("Phi")
hold on
h = plot(x*L, normalizedPhi_an, '--','LineWidth', 1.5,'color', [0.4660, 0.6740, 0.1880]);
hold off
legend([lineInstances{2}, h], ["FEM", "Analytical"],'location' ,'northwest')

end


function lineInstances = FEMplot(node_z, uVector, umax, tetamax, fimax)


    lineInstances = cell(3,1);
    % Functions N are necessary to get the shape of the beam between nodes

    Nfun = @(x) [ 1      - 3*x.^2 + 2*x.^3 ; ...
                 -x      + 2*x.^2 - x.^3   ; ...
                  3*x.^2  - 2*x.^3         ; ...
                  x.^2    - x.^3          ];

    Bfun = @(x) [ -6*x + 6*x.^2             ; ...
                  -1   + 4*x    - 3*x.^2    ; ...
                  6*x  - 6*x.^2             ; ...
                  2*x  - 3*x.^2             ];
              
    Fifun = @(x) [1-x; x];
    
    subplot(3,1,1)
    hold on
    for i=1:length(node_z)
        if i<length(node_z)
            h = node_z(i+1) - node_z(i);
            z = linspace(node_z(i), node_z(i+1), 20);
            lineInstances{1} = plot(z, Nfun(linspace(0,1,20))' * ([uVector(3*(i-1)+[1,2]); uVector(3*i+[1,2])].*[1;h;1;h])/umax,'LineWidth', 2, 'color', [0, 0.4470, 0.7410]);
        end
        scatter(node_z(i), uVector(3*(i-1)+1)/umax,'filled','MarkerFaceColor', [0, 0.4470, 0.7410])
    end
    hold off
    xlabel("z")
    ylabel("u/u_{max}")
    grid minor
    
    
    subplot(3,1,2)
    hold on
    for i=1:length(node_z)
        if i<length(node_z)
            h = node_z(i+1) - node_z(i);
            z = linspace(node_z(i), node_z(i+1), 20);
            lineInstances{2} = plot(z, Bfun(linspace(0,1,20))' * ([uVector(3*(i-1)+[1,2]); uVector(3*i+[1,2])].*[1/h;1;1/h;1])/tetamax, 'LineWidth', 2, 'color', [0, 0.4470, 0.7410]);
        end
        scatter(node_z(i), -uVector(3*(i-1)+2)/tetamax,'filled','MarkerFaceColor', [0, 0.4470, 0.7410])
    end
    hold off
    xlabel("z")
    ylabel("\theta/\theta_{max}")
    grid minor
    
    
    subplot(3,1,3)
    hold on
    for i=1:length(node_z)
        if i<length(node_z)
            z = linspace(node_z(i), node_z(i+1), 20);
            lineInstances{3} = plot(z, Fifun(linspace(0,1,20))' * ([uVector(3*(i-1)+3); uVector(3*i+3)])/fimax,'LineWidth', 2, 'color', [0, 0.4470, 0.7410]);
        end
        scatter(node_z(i), uVector(3*(i-1)+3)/fimax,'filled','MarkerFaceColor', [0, 0.4470, 0.7410])
    end
    hold off
    xlabel("z")
    ylabel("\phi/\phi_{max}")
    grid minor

end
