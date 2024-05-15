
% Set up initial variables
matSizeVector = [101 101]; % The initial matrix size
initialProteins = 0; % The maximum number of proteins to insert before any diffusion, growth, or division.
diffusionCoefficient = 5*10^(-3); % Units are micrometers^2/seconds
cellDiam = 1; % Cell diameter in micrometers
proteinSizeInt = 10; % Interger "radial" size of the protein counted from center(inclusive) to horizontal or vertical edge. Units are the same as the matrix. Use intergers greater than or equal to 0
maxTime = 50000; % Time in iterations
printTime = 10000; % Number of iterations to do before outputing intermediate data and protein locations
growthRate = .0307/60; % Units might be micrometers/second
rateProtein = 1350; % Protein insertion rate per division
proteinChance = rateProtein/1350; % Protein chance is protein insertion rate divided by the expected average time in seconds for cell to divide
initProteinsRandom = 1; % If 1, randomizes the initial protein insertion points, else (suggest set to 0) inserts proteins at center of matrix
doGrowthAndDiv = 1; % If 1, avgEffInsertRatedoes do growth or division of cell, else (suggest set to 0) doesn't do growth or division of cell
doDiff = 1; % If 1, does do diffusion of proteins, else (suggest set to 0) doesn't do diffusion of proteins
continueRun = 1; % If 1, continue the previeous run

% Uncomment the line below in order to make sure proteins don't overlap, untested
%proteinSizeInt = 2*proteinSizeInt - 1; % Make sure proteins don't overlap, untested

% Calls the function diffuseProtein
diffuseProtein(matSizeVector, initialProteins, proteinChance, growthRate, diffusionCoefficient, cellDiam, proteinSizeInt, maxTime, printTime, doGrowthAndDiv, doDiff, initProteinsRandom, continueRun);
