% TODO: User Input

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


if exist('asymp.csv', 'file') == 2 & exist('timeDivVector.txt', 'file') == 2 & exist('analysisVectors.csv', 'file') == 2 & exist('analysisVectors_prev.csv', 'file') == 2 & exist('currIter.txt', 'file') == 2
	% Import already made arrays
	expAsym = readtable('asymp.csv');
	expAsymFarther = expAsym.expAsymFarther;
	expAsymClosest = expAsym.expAsymClosest;
	clear expAsym
	
	timeDivVector = readmatrix('timeDivVector.txt');
	
	analysisVectors = readtable('analysisVectors.csv');
	if size(analysisVectors.timeIterVector, 1)	~= iterCount
		analysisVectors = readtable('analysisVectors_prev.csv');
		if size(analysisVectors.timeIterVector, 1)	~= iterCount
			currIter = -1;
			writematrix(currIter);
		end
	end
	timeIterVector = analysisVectors.timeIterVector;
	timeSecVector = analysisVectors.timeSecVector;
	numProteinsOverTime = analysisVectors.numProteinsOverTime;
	cellLenVector = analysisVectors.cellLenVector;
	cellCircumVector = analysisVectors.cellCircumVector;
	cellSurfAreaVector = analysisVectors.cellSurfAreaVector;
	averageDistanceOverTime = analysisVectors.averageDistanceOverTime;
	densityVec = analysisVectors.densityVec;
	clear analysisVectors
	
	currIter = readmatrix('currIter.txt');
else
	% Bounds for optimal packing of proteins
	expAsymFarther = sqrt(2*(proteinMicroRadiusOver)^2); % Units are micrometers
		% Proteins can be farther, but then not optimal
		% Expected distance if proteins were in a pattern similar to the 5 on a die, dashed line
	expAsymClosest = proteinMicroRadiusOver; % Units are micrometers
		% Proteins cannot be closer than this 
		% Expected distance if proteins were in a cross pattern, dotted line
	writetable(table(expAsymFarther, expAsymClosest), 'asymp.csv');
	
	% Create array for time of divisions
	previousTime = 0;
	timeDivVector = NaN([(size(numDivVector,1) - 1) 1]);
	for j = 1:size(numDivVector,1)
		if numDivVector(j) > 0
			timeDivVector(j-1) = (previousTime + timeBetweenDivVector(j)) * delT;
			previousTime = previousTime + timeBetweenDivVector(j);
		end
	end
	writematrix(timeDivVector);
	
	% Set up time in seconds
	timeIterVector = (1:iterCount)' - 1;
	timeSecVector = timeIterVector * delT;

	% Set up size in micrometers
	cellLenVector = mSizeVector * gridSize;
	cellCircumVector = nSizeVector * gridSize;
	cellSurfAreaVector = 2 * pi * cellRad * cellLenVector;

	% Create needed empty arrays
	numProteinsOverTime = NaN([iterCount 1]); % Number of proteins over time
	averageDistanceOverTime = NaN([iterCount 1]); % Average distance between proteins over time
	densityVec = NaN([iterCount 1]); % Density of proteins over time
	
	writetable(table(timeIterVector, timeSecVector, numProteinsOverTime, cellLenVector, cellCircumVector, cellSurfAreaVector, averageDistanceOverTime, densityVec), 'analysisVectors.csv');
	
	currIter = 1;
	writematrix(currIter);
end
if currIter == -1
		% Bounds for optimal packing of proteins
	expAsymFarther = sqrt(2*(proteinMicroRadiusOver)^2); % Units are micrometers
		% Proteins can be farther, but then not optimal
		% Expected distance if proteins were in a pattern similar to the 5 on a die, dashed line
	expAsymClosest = proteinMicroRadiusOver; % Units are micrometers
		% Proteins cannot be closer than this 
		% Expected distance if proteins were in a cross pattern, dotted line
	writetable(table(expAsymFarther, expAsymClosest), 'asymp.csv');
	
	% Create array for time of divisions
	previousTime = 0;
	timeDivVector = NaN([(size(numDivVector,1) - 1) 1]);
	for j = 1:size(numDivVector,1)
		if numDivVector(j) > 0
			timeDivVector(j-1) = (previousTime + timeBetweenDivVector(j)) * delT;
			previousTime = previousTime + timeBetweenDivVector(j);
		end
	end
	writematrix(timeDivVector);
	
	% Set up time in seconds
	timeIterVector = (1:iterCount)' - 1;
	timeSecVector = timeIterVector * delT;

	% Set up size in micrometers
	cellLenVector = mSizeVector * gridSize;
	cellCircumVector = nSizeVector * gridSize;
	cellSurfAreaVector = 2 * pi * cellRad * cellLenVector;

	% Create needed empty arrays
	numProteinsOverTime = NaN([iterCount 1]); % Number of proteins over time
	averageDistanceOverTime = NaN([iterCount 1]); % Average distance between proteins over time
	densityVec = NaN([iterCount 1]); % Density of proteins over time
	
	writetable(table(timeIterVector, timeSecVector, numProteinsOverTime, cellLenVector, cellCircumVector, cellSurfAreaVector, averageDistanceOverTime, densityVec), 'analysisVectors.csv');
	
	currIter = 1;
	writematrix(currIter);
end
clear numDivVector

% Analyze

printTime = 1000;
if currIter > 1
	iPrintTime = -1;
else
	iPrintTime = 0;
end

for j = currIter:iterCount

	iPrintTime = iPrintTime + 1;

	averageDistance = 0;
	numProteins = 0;
	if iterNumVector(iterNumVector == j-1)
		numProteins = size(iterNumVector(iterNumVector == j-1), 1); % Number of proteins
		if numProteins > 1
			mIterPosVector = mMicroPosVector(iterNumVector == j-1); % M locations for this iteration
			thetaIterPosVector = thetaPosVector(iterNumVector == j-1);% Theta locations for this iteration
			for i = 1:numProteins
				hVector = NaN([(numProteins - 1) 1]); % Helper vector for calculating average distance over time
				for k = 1:numProteins
					if k ~= i
						if thetaIterPosVector(i) < thetaIterPosVector(k)
							phi_1 = thetaIterPosVector(k) - thetaIterPosVector(i);
							phi_2 = 2*pi + thetaIterPosVector(i) - thetaIterPosVector(k);
						else
							phi_1 = thetaIterPosVector(i) - thetaIterPosVector(k);
							phi_2 = 2*pi + thetaIterPosVector(k) - thetaIterPosVector(i);
						end
						mDist = mIterPosVector(i) - mIterPosVector(k);
						if phi_1 < phi_2
							if k < i
								hVector(k) = sqrt((phi_1 * cellRad)^2 + mDist^2);
							else
								hVector(k-1) = sqrt((phi_1 * cellRad)^2 + mDist^2);
							end
						else
							if k < i
								hVector(k) = sqrt((phi_2 * cellRad)^2 + mDist^2);
							else
								hVector(k-1) = sqrt((phi_2 * cellRad)^2 + mDist^2);
							end
						end
					end
				end
				% Average closest 4 distances, then add to the average distance
				hClosest = mink(hVector, 4);
				averageDistance = averageDistance + sum(hClosest)/size(hClosest,1);
			end
			% Average of the added averages of the closest 4 proteins over all proteins
			averageDistance = averageDistance / numProteins;
		else
			averageDistance = NaN;
		end
	else
		numProteins = 0; % Number of proteins
		averageDistance = NaN;
	end
	numProteinsOverTime(j) = numProteins;
	averageDistanceOverTime(j) = averageDistance;
	
	m = mSizeVector(j);
	n = nSizeVector(j);
	densityVec(j) = numProteins/(m * n * gridSize^2);
	
	if iPrintTime == printTime
		writetable(table(timeIterVector, timeSecVector, numProteinsOverTime, cellLenVector, cellCircumVector, cellSurfAreaVector, averageDistanceOverTime, densityVec), 'analysisVectors.csv');
		writetable(table(timeIterVector, timeSecVector, numProteinsOverTime, cellLenVector, cellCircumVector, cellSurfAreaVector, averageDistanceOverTime, densityVec), 'analysisVectors_prev.csv');
		if size(timeIterVector, 1)	~= iterCount
			currIter = -1;
			writematrix(currIter);
			error('Size of timeIterVector incorrect.')
		else
			currIter = j;
			writematrix(currIter);
		end
		iPrintTime = 0;
	end
end
if size(timeIterVector, 1)	~= iterCount
	currIter = -1;
	writematrix(currIter);
	error('Size of timeIterVector incorrect.')
writetable(table(timeIterVector, timeSecVector, numProteinsOverTime, cellLenVector, cellCircumVector, cellSurfAreaVector, averageDistanceOverTime, densityVec), 'analysisVectors.csv');

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

close all