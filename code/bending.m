function [defl,teta,fi,umax,tetamax,fimax]=bending(Ks,Qs,K,Q,nnode,node_z);

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
    
    [umax, tetamax, fimax] = uReducedVector(end-2:end);

    % Present reaction forces
    
    RF = K * [0;0;0; uReducedVector];

    % Create result vector containing deflections, rotations and twist
    
    uVector = [0;0;0; uReducedVector];

    % Split deflections, rotations and twist into separate vectors
    defl = uVector(1:3:end);
    teta = uVector(2:3:end);
    fi   = uVector(3:3:end);

    % Normalise deflections, rotations and twist and plot results
    
    normalizedDefl = def / umax;
    normalizedTeta = teta / tetamax;
    normalizedFi = fi / fimax;
    
    figure
    
    subplot(3,1,1)
    title("Deflection")
    scatter(node_z, normalizedDefl)
    xlabel("z")
    ylabel("u/u_max")
    minor grid
  
    subplot(3,1,2)
    title("Teta")
    scatter(node_z, normalizedTeta)
    xlabel("z")
    ylabel("u/u_max")
    minor grid
    
    subplot(3,1,3)
    title("Fi")
    scatter(node_z, normalizedFi)
    xlabel("z")
    ylabel("u/u_max")
    minor grid
end
