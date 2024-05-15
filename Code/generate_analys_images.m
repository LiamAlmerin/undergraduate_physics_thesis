% Import initial conditions
initialConditionsTable = readtable('initial_conditions.csv');
% Set up variables
diffCoef = initialConditionsTable.Diffusion_Coefficient;
gridSize = initialConditionsTable.Grid_Size;
delT = initialConditionsTable.Delta_t;
cellDiam = initialConditionsTable.Cell_Diameter;
growthRate = initialConditionsTable.Growth_Rate;
proteinRadiusOver = initialConditionsTable.Protein_Radius_Overlapping;
proteinChanceIter = initialConditionsTable.Protein_Chance_Per_Iteration;
clear initialConditionsTable

% Import sim data
simDataTable = readtable('sim_data.csv');
% Set up arrays
mSizeVector = simDataTable.M_Size;
nSizeVector = simDataTable.N_Size;
addProteinVector = simDataTable.Added_Proteins;
%xAvgVector = simDataTable.X_Average;
%yAvgVector = simDataTable.Y_Average;
%xSqAvgVector = simDataTable.X_Square_Average;
%ySqAvgVector = simDataTable.Y_Square_Average;
%xDiffVector = simDataTable.X_Diffusion_Value;
%yDiffVector = simDataTable.Y_Diffusion_Value;
clear simDataTable

% Import div data
divDataTable = readtable('div_data.csv');
% Set up arrays
sizePreDivVector = divDataTable.Size_Pre_Division;
sizePostDivVector = divDataTable.Size_Post_Division;
timeBetweenDivVector = divDataTable.Time_Between_Divisions;
numDivVector = divDataTable.Number_of_Divisions;
clear divDataTable

% Import protein locations
proteinLocDataTable = readtable('protein_loc.csv');
% Set up arrays
iterNumVector = proteinLocDataTable.Iteration_Number;
mPosVector = proteinLocDataTable.M_Position;
nPosVector = proteinLocDataTable.N_Position;
clear proteinLocDataTable

% Precalculations

proteinMicroRadiusOver = proteinRadiusOver * gridSize; % Radius for overlapping protein in micrometers
iterCount = size(mSizeVector,1); % Total iterations done

% Create vector for the theta positions of the proteins
cellRad = cellDiam / 2; % Cell radius in micrometers
theta = gridSize/cellRad; % Grid size of cell in radians
thetaPosVector = nPosVector * theta; % Vector for the theta positions

% Set up m positions in micrometers
mMicroPosVector = mPosVector * gridSize;


if exist('asymp.csv', 'file') == 2 & exist('timeDivVector.txt', 'file') == 2 & exist('analysisVectors.csv', 'file') == 2
	% Import already made arrays
	expAsym = readtable('asymp.csv');
	expAsymFarther = expAsym.expAsymFarther;
	expAsymClosest = expAsym.expAsymClosest;
	clear expAsym
	
	timeDivVector = readmatrix('timeDivVector.txt');
	
	analysisVectors = readtable('analysisVectors.csv');

	timeIterVector = analysisVectors.timeIterVector;
	timeSecVector = analysisVectors.timeSecVector;
	numProteinsOverTime = analysisVectors.numProteinsOverTime;
	cellLenVector = analysisVectors.cellLenVector;
	cellCircumVector = analysisVectors.cellCircumVector;
	cellSurfAreaVector = analysisVectors.cellSurfAreaVector;
	averageDistanceOverTime = analysisVectors.averageDistanceOverTime;
	densityVec = analysisVectors.densityVec;
	clear analysisVectors

else
	x = -1
end
% Draw figue
figure('Name', 'Average Distance of Closest Proteins' , 'Position', [0, 0, 1920, 1080],'NumberTitle','off');
hold on

% Average distance over time graph
title('Average Distance of Closest Proteins')
xlabel('Time(in seconds)')
ylabel('Average Distance (in micrometers)')
plot(timeSecVector,averageDistanceOverTime, '-b'); % Solid
yline(expAsymFarther, '--b'); % Dashed
yline(expAsymClosest, ':b'); % Dotted
set(gca, 'FontSize', 15);
print('avgDistClosest.png','-dpng', '-r300');
hold off

% Number of proteins over time graph
figure('Name', 'Number of Proteins' , 'Position', [0, 0, 1920, 1080],'NumberTitle','off');
hold on

title('Protein Count')
xlabel('Time(in seconds)')
ylabel('Number of Proteins')
plot(timeSecVector, numProteinsOverTime, '-r');
xline(timeDivVector, '--k');
set(gca, 'FontSize', 15);
print('proteinCount.png','-dpng', '-r300');
hold off

% Density of proteins over time graph
figure('Name', 'Density of Proteins' , 'Position', [0, 0, 1920, 1080],'NumberTitle','off');
hold on

title('Protein Density')
xlabel('Time(in seconds)')
ylabel('Density of Proteins (in proteins per micrometers^2)')
plot(timeSecVector, densityVec, '-r');
xline(timeDivVector, '--k');
set(gca, 'FontSize', 15);
print('proteinDensity.png','-dpng', '-r300');
hold off

hold on

yyaxis left
title('Protein Density and Cell Length over Time')
xlabel('Time(in seconds)')
ylabel('Density of Proteins (in proteins per micrometers^2)', 'Color', 'r')
plot(timeSecVector, densityVec, '-r');

yyaxis right
ylabel('Cell Length (in micrometers)', 'Color', 'b')
plot(timeSecVector, cellLenVector, '-b');
xline(timeDivVector, '--k');

ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';

set(gca, 'FontSize', 15);
print('proteinDensityLen.png','-dpng', '-r300');
hold off

close all