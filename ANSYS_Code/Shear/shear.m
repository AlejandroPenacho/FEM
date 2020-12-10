clc; clear

%% Extract data
rawPos = [];

fID = fopen("nPos.txt");
fgetl(fID);
fgetl(fID);
fgetl(fID);
fgetl(fID);

while ~feof(fID)
    fgetl(fID);
    rawPos = [rawPos; fscanf(fID,"%d %f %f %f")];
end

fclose(fID);

nNodes = length(rawPos)/4;

positionArray = zeros(nNodes, 3);
positionArray(:,1) = rawPos(1:4:end);
positionArray(:,2) = rawPos(2:4:end);
positionArray(:,3) = rawPos(3:4:end);


%% Stress data

zRawData = [];

fID = fopen("stress.txt");
for i=1:3
    fgetl(fID);
end

done = false;

while ~done
    for i=1:10
        fgetl(fID);
    end
    
    newData = fscanf(fID,"%d %f %f %f %f %f %f");
    
    if newData
        zRawData = [zRawData; newData];
    else
        done = true;
    end
    
end

fclose(fID);

nZData = length(zRawData)/7;

zData = zeros(nZData, 7);

for i=1:7
    zData(:,i) = zRawData(i:7:end);
end

% Assign stresses

finalData = zeros(nNodes, 5);

finalData(:,1:3) = positionArray;

for i = 1:nNodes
    stressValues = zData(zData(:,1) == finalData(i,1),:);
    
    if length(stressValues(:,1)) == 4
        finalData(i,4) = -mean(stressValues(2:3,6));
        finalData(i,5) = -mean(stressValues(2:3,7));        
    else
        finalData(i,4) = mean(stressValues(:,6));
        finalData(i,5) = mean(stressValues(:,7));
    end
end



%% Sort the nodes

firstFlange = finalData(positionArray(:,3) == 0 & positionArray(:,2) ~= 0,:);
secondFlange = finalData(positionArray(:,2) == 0, :); 
thirdFlange = finalData(positionArray(:,3) == 40 & positionArray(:,2) ~= 0, :);

[~, firstFlangeIndex] = sort(firstFlange(:,2), 'descend');
[~, secondFlangeIndex] = sort(secondFlange(:,3));
[~, thirdFlangeIndex] = sort(thirdFlange(:,2));

finalData(:,1:3) = [firstFlange(firstFlangeIndex,1:3); 
             secondFlange(secondFlangeIndex,1:3); 
             thirdFlange(thirdFlangeIndex,1:3)];

finalData(:,4) = [firstFlange(firstFlangeIndex,5); 
             -secondFlange(secondFlangeIndex,4); 
             -thirdFlange(thirdFlangeIndex,5)];

finalData = finalData(~isnan(finalData(:,4)),:);

nNodes = length(finalData(:,1));
         
%% Analyitical results

S = 100;
h = 40;
t_f = 3;
t_w = 1;
b = 40;
I_xx = 101333.333;

m_f = h*t_f / (2*I_xx) * S / 3;

w_coeff = (S/(2*I_xx)) * [-t_w, h*t_w, h*b*t_f];


         
%% Plot results
nPoints = 100;

% figure
% scatter3(-(finalData(:,2)-40), zeros(nNodes,1), finalData(:,3), "filled")
% hold on
% h(1) = plot3(-(finalData(:,2)-40), finalData(:,4), finalData(:,3), "LineWidth", 2);
% 
% h(2) = plot3(linspace(0,40,nPoints), linspace(0,40,nPoints) * m_f, zeros(nPoints,1), "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
% plot3(linspace(0,40,nPoints), linspace(0,40,nPoints) * m_f, 40 + zeros(nPoints,1), "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
% plot3(zeros(nPoints,1)+40,...
%     w_coeff(1) * linspace(0,40,nPoints).^2 + w_coeff(2) * linspace(0,40,nPoints) + w_coeff(3) * ones(1, nPoints), ...
%     linspace(0,40,nPoints), ...
%       "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);

figure
scatter3(finalData(:,2), zeros(nNodes,1), finalData(:,3), "filled")
hold on
h(1) = plot3(finalData(:,2), finalData(:,4), finalData(:,3), "LineWidth", 2);

h(2) = plot3(linspace(40,0,nPoints), linspace(0,40,nPoints) * m_f, zeros(nPoints,1), "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
plot3(linspace(40,0,nPoints), linspace(0,40,nPoints) * m_f, 40 + zeros(nPoints,1), "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
plot3(zeros(nPoints,1), ...
    w_coeff(1) * linspace(0,40,nPoints).^2 + w_coeff(2) * linspace(0,40,nPoints) + w_coeff(3) * ones(1, nPoints), ...
    linspace(0,40,nPoints), ...
      "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
  view(135,20);
hold off

title("Shear stress distribution at z=L/2")
xlabel("X (mm)")
ylabel("q (N/mm^2)")
zlabel("Y (mm)")

legend(h, "ANSYS", "Analytical")
