clc; clear

%% Extract data1

fID = fopen("data1.txt");
for i=1:6
    fgetl(fID);
end

posData = fscanf(fID,"%f %f %f");
fclose(fID);

nNodes = length(posData)/3;

dataArray1 = zeros(nNodes, 2);
dataArray1(:,1) = posData(2:3:end);
dataArray1(:,2) = posData(3:3:end);


%% Extract data2

fID = fopen("data2.txt");
for i=1:6
    fgetl(fID);
end

posData = fscanf(fID,"%f %f %f");
fclose(fID);

nNodes = length(posData)/3;

dataArray2 = zeros(nNodes, 2);
dataArray2(:,1) = posData(2:3:end);
dataArray2(:,2) = posData(3:3:end);


%% Analytical
S = 100;
h = 40;
Ixx = 304000/3;

maxStress = @(z) S*h/(2*Ixx) * (1000-z);

%% PLOT

figure
hold on
plot(dataArray1(:,1), dataArray1(:,2), "LineWidth", 2);
plot(dataArray2(:,1), dataArray2(:,2), "LineWidth", 2);
plot(0:0.1:1000, maxStress(0:0.1:1000), "LineWidth", 2);
yline(0);
hold off
grid minor

legend("ANSYS (X=0)", "ANSYS (X=40)", "Analytical")