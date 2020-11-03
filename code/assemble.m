function [K,Q,M,Ksigma]=assemble(le,EI,GJ,I0,A,J0,q,qt,S,T,m,P,ndof,nelem);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble system stiffness matrix, load vector, mass matrix (not necessary)
    % and initial stress matrix
    % File name: assemble.m
    %
    % K		System stiffness matrix
    % Q		System load vector
    % M		System mass matrix
    % Ksigma       System initial stress matrix
    %
    % le		element length [m]
    % EI		element bending stiffness [Nm2]
    % GJ		element torsional stiffness [Nm2]
    % I0		element polar moment of inertia [m4]
    % A		element cross-section area [m2]
    % J0		element mass moment of inertia [kgm]
    % q		element transverse pressure load [N/m]
    % qt		element torsion pressure load [Nm/m]
    % S		transverse point load [N]
    % T		local torque [Nm]
    % m		element mass per unit length [kg/m]
    % P		applied buckling load [N], positive if tensile
    % ndof		number of degrees of freedom
    % nelem		number of elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K=zeros(ndof);
    Q=zeros(ndof,1);
    M=zeros(ndof);
    Ksigma=zeros(ndof);
    
    %% Assemble of stiffness matrix
    
    % For every element, the element stiffness matrix is computed and added
    % to the total stiffness matrix of the structure
    for iElement = 1:nelem
        
        % For the case of a simple beam, the nodes at the tips of each
        % element are easily computed. Nodes are explicitly calculated so
        % changes can be done easily to this process.
        
        firstNode = iElement;
        secondNode = iElement + 1;
        
        
        firstNodeBaseIndex = (firstNode-1) * 3;
        secondNodeBaseIndex = (secondNode-1) * 3;
        
        [Ke]=elk(le,EI,GJ);
        [Kesigma]=elksigma(le,P,I0,A);
        [Qe]=elq(le,q,qt);              
        
        deltaK = zeros(ndof);
        deltaKsigma = zeros(ndof);
        deltaQload = zeros(ndof,1);
        
        for i=1:6
            for j=1:6
                if i<=3
                    baseI = firstNodeBaseIndex +i;
                else
                    baseI = secondNodeBaseIndex +i - 3;
                end
                if j<=3
                    baseJ = firstNodeBaseIndex + j;
                else
                    baseJ = secondNodeBaseIndex +j - 3;
                end                
                
                deltaK(baseI, baseJ) = Ke(i,j);
                deltaKsigma(baseI, baseJ) = Kesigma(i,j);
                deltaQload(baseI) = Qe(i);
            end
        end
        K = K + deltaK;
        Ksigma = Ksigma + deltaKsigma;
        Q = Q + deltaQload;
    end




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add concentrated loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

