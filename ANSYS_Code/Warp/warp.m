clc; clear

%% Extract data

fID = fopen("nPos.txt");
for i=1:5
    fgetl(fID);
end

posData = fscanf(fID,"%d %f %f %f");
fclose(fID);

nNodes = length(posData)/4;

dataArray = zeros(nNodes, 3);
dataArray(:,1) = posData(2:4:end);
dataArray(:,2) = posData(3:4:end);


fID = fopen("zDisp.txt");
for i=1:11
    fgetl(fID);
end

zData = fscanf(fID,"%d %f");
fclose(fID);

dataArray(:,3) = zData(2:2:end);



%% Substract rotation

fID = fopen("xRot.txt");
for i=1:11
    fgetl(fID);
end

rotData = fscanf(fID,"%d %f");
fclose(fID);

rotation = mean(rotData(2:2:end));

rotation = 0.75 * rotation;

dataArray(:,3) = dataArray(:,3) - rotation * (dataArray(:,2)-20);

%% Sort the nodes

firstFlange = dataArray(dataArray(:,2) == 0 & dataArray(:,1) ~= 0,:);
secondFlange = dataArray(dataArray(:,1) == 0, :); 
thirdFlange = dataArray(dataArray(:,2) == 40 & dataArray(:,1) ~= 0, :);

[~, firstFlangeIndex] = sort(firstFlange(:,1), 'descend');
[~, secondFlangeIndex] = sort(secondFlange(:,2));
[~, thirdFlangeIndex] = sort(thirdFlange(:,1));

dataArray = [firstFlange(firstFlangeIndex,:); 
             secondFlange(secondFlangeIndex,:); 
             thirdFlange(thirdFlangeIndex,:)];

         
%% Analyitical results

T = 100 * 360/19;
J = 2200/3;
G = 26923;

limitW1 = 20 * (360/19) * T/(G*J);

limitW2 = -20 * (360/19 - 40) * T/(G*J);

         
%% Plot results
figure
scatter3(dataArray(:,1), zeros(nNodes,1), dataArray(:,2), "filled")
hold on
h(1) = plot3(dataArray(:,1), dataArray(:,3), dataArray(:,2), "LineWidth", 2);

h(2) = plot3([0,40], [-limitW1, limitW2], [0,0], "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
plot3([0,0], [-limitW1, limitW1], [0,40], "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
plot3([0,40], [limitW1, -limitW2], [40,40], "color", [0.4940, 0.1840, 0.5560], "LineWidth", 2);
view(135,20);
hold off

title("Warping distribution at z=L")
xlabel("X (mm)")
ylabel("w (mm)")
zlabel("Y (mm)")


legend(h, "ANSYS", "Analytical")
