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
    [ub,pb] = eig(Ks, -Ksigmas);

    [ndof, ~] = size(Ks);

    % Split bending and twist modes into separate vectors
    
    % The "vector" variables are matrices with the vectors of the shapes of
    % the different modes, while the "lambda" variables are the eigenvalues
    % of those vectors.

    bendingVectors = zeros(ndof);
    bendingLambdas = zeros(ndof, 1);
    nBendingModes = 0;

    twistVectors = zeros(ndof);
    twistLambdas = zeros(ndof, 1);
    nTwistModes = 0;

    % To determine whether a mode is related to bending or twist, the value
    % of the maximum value of phi is checked. If it is below a certain
    % threhold, it is considered to be zero, so it is a bending mode.
    
    threshold = 1E-10;

    for i=1:ndof
        if max(abs(ub(3:3:end,i))) < threshold
            nBendingModes = nBendingModes + 1;
            bendingVectors(:,nBendingModes) = ub(:,i);
            bendingLambdas(nBendingModes) = pb(i,i);
        else
            nTwistModes = nTwistModes + 1;
            twistVectors(:,nTwistModes) = ub(:,i);
            twistLambdas(nTwistModes) = pb(i,i);        
        end
    end
    
    % Vectors are trimmed

    bendingVectors = bendingVectors(:,1:nBendingModes);
    bendingLambdas = bendingLambdas(1:nBendingModes);
    twistVectors = twistVectors(:,1:nTwistModes);
    twistLambdas = twistLambdas(1:nTwistModes);
    
    
    % Normalise deflections, rotations and twist for plotting purposes (without risking to mix up signs or divide by zero)

    if nnode == 2
        nBendingModesPlotting = 1;
    else
        nBendingModesPlotting = 2;   
    end
    
    %% Bending Modes
    
    
    bendingLegend = cell(nBendingModesPlotting,1);
    
    for i=1:nBendingModesPlotting
        bendingLegend{i} = sprintf("N_{z}= %.2f kN", bendingLambdas(i)/1000);
    end
    
    % Max displacements and rotations are obtained
    
    maxDisplacementVector = max(abs(bendingVectors(1:3:end,:)));
    maxRotationVector = max(abs(bendingVectors(2:3:end,:)));
    
    figure

    sgtitle("Bending modes")
    
    colorCell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980]};
    for i=1:nBendingModesPlotting
        lineInstances(i) = FEMBucklePlot(node_z, [0;0;0; bendingVectors(:,i)], ...
                                        maxDisplacementVector(i), maxRotationVector(i), colorCell{i});
    end
    
    subplot(2,1,1)
    legend(lineInstances, bendingLegend, 'location' ,'northwest')
    
%     subplot(2,1,1)
%     title("Displacement")
%     
%     hold on
%     for i = 1:nBendingModesPlotting
%         scatter(node_z, [0; bendingVectors(1:3:end,i)/maxDisplacementVector(i)], 'filled')
%     end
%     hold off
%     grid minor
%     legend(bendingLegend, "location", "northwest")
%     xlabel("z")
%     ylabel("u/u_{max}")
%     
%     
%     subplot(2,1,2)
%     title("Rotation")
%     hold on
%     for i = 1:nBendingModesPlotting
%         scatter(node_z, [0; bendingVectors(2:3:end,i)/maxRotationVector(i)], 'filled')
%     end
%     hold off
%     grid minor
%     xlabel("z")
%     ylabel("\theta/\theta_{max}")
    
    
    % Results are written to a text file
    
    fID = fopen("output/results.txt", "a");
    
    fprintf(fID, "\n\n\nBending buckling modes\n\n");
    
    fprintf(fID, "Mode nº");
    for i=1:nBendingModes
        fprintf(fID, "\t\t\t%d", i);
    end
    
    fprintf(fID, "\nLoad (kN)");
    for i=1:nBendingModes
        fprintf(fID, "\t\t%3f", bendingLambdas(i)/1000);
    end    
    
    
    
    %% Twist Modes
    
    fprintf(fID, "\n\n\nTwist buckling\n\n");
    if (max(twistLambdas)-min(twistLambdas))/max(twistLambdas) < 10^-10
        fprintf(fID, "All eigenvalues for twist buckling coincide\n");
        fprintf(fID, "Load (kN): %.3f", twistLambdas(1)/1000);
    end
      
    fclose(fID);
end


function lineInstances = FEMBucklePlot(node_z, uVector, umax, tetamax, color)


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

    subplot(2,1,1)
    hold on
    for i=1:length(node_z)
        if i<length(node_z)
            h = node_z(i+1) - node_z(i);
            z = linspace(node_z(i), node_z(i+1), 20);
            lineInstances = plot(z, Nfun(linspace(0,1,20))' * ([uVector(3*(i-1)+[1,2]); uVector(3*i+[1,2])].*[1;h;1;h])/umax,'LineWidth', 2, 'color', color);
        end
        scatter(node_z(i), uVector(3*(i-1)+1)/umax,'filled','MarkerFaceColor', color)
    end
    hold off
    xlabel("z")
    ylabel("u/u_{max}")
    grid on
    
    
    subplot(2,1,2)
    hold on
    for i=1:length(node_z)
        if i<length(node_z)
            h = node_z(i+1) - node_z(i);
            z = linspace(node_z(i), node_z(i+1), 20);
            plot(z, Bfun(linspace(0,1,20))' * ([uVector(3*(i-1)+[1,2]); uVector(3*i+[1,2])].*[1/h;1;1/h;1])/tetamax, 'LineWidth', 2, 'color', color);
        end
        scatter(node_z(i), -uVector(3*(i-1)+2)/tetamax,'filled','MarkerFaceColor', color)
    end
    hold off
    xlabel("z")
    ylabel("\theta/\theta_{max}")
    grid on
    
end


% Plot buckling modes

% Present the buckling loads






