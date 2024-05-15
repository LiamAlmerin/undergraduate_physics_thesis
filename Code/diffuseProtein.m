function diffuseProtein(matSizeVector, initialProteins, proteinChance, growthRate, diffusionCoefficient, cellDiam, proteinSizeInt, maxTime, printTime, doGrowthAndDiv, doDiff, initProteinsRandom, continueRun)

% Set up variables
if continueRun == 1
	% Import initial conditions
	initialConditionsTable = readtable('initial_conditions.csv');
	% Set up variables
	diffusionCoefficient = initialConditionsTable.Diffusion_Coefficient;
	gridSize = initialConditionsTable.Grid_Size;
	delT = initialConditionsTable.Delta_t;
	cellDiam = initialConditionsTable.Cell_Diameter;
	growthRate = initialConditionsTable.Growth_Rate;
	proteinSizeInt = initialConditionsTable.Protein_Radius_Overlapping;
	proteinChance = initialConditionsTable.Protein_Chance_Per_Iteration;
	clear initialConditionsTable

	% Import sim data
	simDataTable = readtable('sim_data.csv');
	% Set up arrays
	mSizeVector = simDataTable.M_Size;
	nSizeVector = simDataTable.N_Size;
	addedProteinsVector = simDataTable.Added_Proteins;
	clear simDataTable
	
	m = mSizeVector(end);
	n = nSizeVector(end);
	matSizeVector = [m n];
	
	% Set up time in seconds
	iterCount = size(mSizeVector,1); % Total iterations done
	timeIterVector = (1:iterCount)' - 1;
	timeSecVector = timeIterVector * delT;
	ts = timeSecVector(end);

	% Import div data
	divDataTable = readtable('div_data.csv');
	% Set up arrays
	divSizeVector = divDataTable.Size_Pre_Division;
	x0Vector = divDataTable.Size_Post_Division;
	divTimeVector = divDataTable.Time_Between_Divisions;
	divNumVector = divDataTable.Number_of_Divisions;
	clear divDataTable
	
	divNum = divNumVector(end);
	mi = x0Vector(end);
	
	neccesaryValuesTable = readtable('neccesary_values.csv');
	timeSinceLastDiv = neccesaryValuesTable.timeSinceLastDiv;
	proteinTagInt = neccesaryValuesTable.proteinTagInt;
	clear timeSinceLastDivTable
	
	% Import protein locations
	proteinLocDataTable = readtable('protein_loc.csv');
	% Set up arrays
	iterNumVector = proteinLocDataTable.Iteration_Number;
	mPosVector = proteinLocDataTable.M_Position;
	nPosVector = proteinLocDataTable.N_Position;
	proteinTagOverTimeVector = proteinLocDataTable.Protein_Tag;
	clear proteinLocDataTable
	
	iterNum = iterNumVector(end);
	xPosVector = nPosVector(iterNumVector == iterNumVector(end)).'; % X position of protein, related to the n value of the matrix
	yPosVector = mPosVector(iterNumVector == iterNumVector(end)).'; % Y position of protein, related to the m value of the matrix
	proteinTagVector = proteinTagOverTimeVector(iterNumVector == iterNumVector(end)).'; % Vector to 'tag' the proteins
	numProteins = size(proteinTagVector);
	addedProteins = 0; % Number of added proteins
	
	cellAgeMat = readmatrix('cell_age.csv');
	cellVolMat = readmatrix('cell_vol.csv');
	
	i_time = size(cellAgeMat, 2);
	i_cell = size(cellAgeMat, 1);
	i_printTime = 0;
	i_picture = 0;
	hillCoefficient = 12;
else
	xPosVector = []; % X position of protein, related to the n value of the matrix
	yPosVector = []; % Y position of protein, related to the m value of the matrix
	addedProteins = 0; % Number of added proteins
	proteinTagInt = 0; % Tag to know which protein
	proteinTagVector = []; % Vector to 'tag' the proteins
	ts = 0; % Current time in seconds
	i_time = 1; % Current cell iteration number, starts at 1
	iterNum = 0; % Current iteration number, starts at 0
	i_cell = 1; % Current cell number, starts at 1
	i_printTime = 0; % Current iteration number to be reset after reaching printTime
	i_picture = 0; % Current image number
	cellAgeMat(1,1) = 0; % Current age of cell
	divNum = 0; % Number of divisions
	divNumVector = [divNum];
	divSizeVector = [NaN]; % Length of cell at division 
	divTimeVector = [NaN]; % Time between divisions
	timeSinceLastDiv = 0; % Time since last division
	effInsertRate = []; % Effective Insertion Rate in proteins per second
	expInsertRate = []; % Expected Insertion Rate in proteins per second
	avgEffInsertRate = []; % Average Effective Insertion Rate in proteins per second
	iterNumVector = [];
	mPosVector = [];
	nPosVector = [];
	proteinTagOverTimeVector = []; % Vector hold all the tags for each iteration

	% Calculations required and more setting up variables
	cellCircum = cellDiam * pi; % Cell circumference in micrometers
	m = matSizeVector(1); % Matrix component that may change
	mi = m; % Initial m value
	n = matSizeVector(2); % Matrix component that remains unchanged
	mSizeVector = [matSizeVector(1)];
	nSizeVector = [matSizeVector(2)];
	gridSize = cellCircum / n; % The size of the matrix grid in micrometers
	x0Vector = [mi*gridSize]; % Initial cell size
	delT = (gridSize)^2 /(2 * diffusionCoefficient); % Time in seconds between each iteration
	proteinChance = proteinChance * delT; % Protein chance multiplied by delta t to make the chance per iteration
	hillCoefficient = 12;
	volume0 = pi * ((cellDiam/2)^2)*(gridSize * mi); % Initial cell volume in micrometers^3
	cellVolMat(1,1) = volume0; % Initial cell volume


	initConditionsTable = table(diffusionCoefficient, gridSize, delT, cellDiam, growthRate, proteinSizeInt, proteinChance, 'VariableNames', {'Diffusion_Coefficient', 'Grid_Size', 'Delta_t', 'Cell_Diameter', 'Growth_Rate', 'Protein_Radius_Overlapping', 'Protein_Chance_Per_Iteration'});
	writetable(initConditionsTable, 'initial_conditions.csv');

	% Insert initial proteins if initialProteins are set to greater than 0
	if initialProteins > 0
		if initProteinsRandom == 1 % Random insertion
			for i = 1:initialProteins
				proteinChanceInit = 1; % Separate proteinChance variable to force a 100% insertion rate
				
				% Calls insertProtein with proteinChance variable swapped with proteinChanceInit
				[xPosVector, yPosVector, addedProteins, proteinTagInt, proteinTagVector] = insertProtein(xPosVector, yPosVector, matSizeVector, proteinSizeInt, addedProteins, proteinChanceInit, proteinTagInt, proteinTagVector);
			end
		else % Insertion at center of matrix
			centerPosVector = [ceil(m/2) ceil(n/2)];
			for i = 1:initialProteins
				xPosVector = [xPosVector centerPosVector(1)];
				yPosVector = [yPosVector centerPosVector(2)];
				addedProteins = addedProteins + 1;
				
			end
		end
	end

	addedProteinsVector = [addedProteins]; % Vector of added proteins
	addedProteins = 0; % Reset added protein count

	numProteins = size(xPosVector,2); % Current total number of proteins
	if numProteins > 0
		for num = 1:numProteins
			iterNumVector = [iterNumVector ; 0];
			mPosVector = [mPosVector ; yPosVector(1, num)];
			nPosVector = [nPosVector ; xPosVector(1, num)];
			proteinTagOverTimeVector = [proteinTagOverTimeVector ; proteinTagVector(1, num)];
		end
	end
	numProteins = 0;% Reset number of proteins

end
rng('shuffle'); % Set seed before simulation
for totalTime = 1:maxTime
	% Increment time
    ts = ts + delT;
    timeSinceLastDiv = timeSinceLastDiv + 1;
    i_time = i_time + 1;
    i_printTime = i_printTime + 1;
    i_picture = i_picture + 1;
    cellAgeMat(i_cell, i_time) = cellAgeMat(i_cell, i_time-1) + 1;
	iterNum = iterNum + 1;
	
    if proteinChance > 0
		% Calls insertProtein
        [xPosVector, yPosVector, addedProteins, proteinTagInt, proteinTagVector] = insertProtein(xPosVector, yPosVector, matSizeVector, proteinSizeInt, addedProteins, proteinChance, proteinTagInt, proteinTagVector);
    end
	
	addedProteinsVector = [addedProteinsVector ; addedProteins];
	addedProteins = 0; % Reset added protein count
	
    if doGrowthAndDiv == 1
		% Calls growCell
        [yPosVector, matSizeVector, cellVolMat] = growCell(yPosVector, matSizeVector, gridSize, growthRate, cellDiam, cellVolMat, i_cell, i_time, delT);
    end
	
    if doDiff == 1
		% Calls diffuseVector twice
        [xPosVector, yPosVector] = diffuseVector(xPosVector, yPosVector, matSizeVector, proteinSizeInt);
		[xPosVector, yPosVector] = diffuseVector(xPosVector, yPosVector, matSizeVector, proteinSizeInt);

    end
    if doGrowthAndDiv == 1
		% Calls divCell
        [xPosVector, yPosVector, matSizeVector, cellAgeMat, cellVolMat, timeSinceLastDiv, divTimeVector, x0Vector, divSizeVector, mi, i_cell, divNum, proteinTagVector] = divCell(xPosVector, yPosVector, matSizeVector, cellAgeMat, cellVolMat, mi, x0Vector, divSizeVector, timeSinceLastDiv, divTimeVector, i_cell, i_time, delT, divNum, gridSize, hillCoefficient, growthRate, cellDiam, proteinTagVector);

	end
	divVectorSize = size(divNumVector, 1);
	if divNum > divNumVector(divVectorSize)
		divNumVector = [divNumVector ; divNum];
	end
	mSizeVector = [mSizeVector ; matSizeVector(1)];
	nSizeVector = [nSizeVector ; matSizeVector(2)];
	
	numProteins = size(xPosVector,2); % Current total number of proteins
	if numProteins > 0
		for num = 1:numProteins
			iterNumVector = [iterNumVector ; iterNum];
			mPosVector = [mPosVector ; yPosVector(1, num)];
			nPosVector = [nPosVector ; xPosVector(1, num)];
			proteinTagOverTimeVector = [proteinTagOverTimeVector ; proteinTagVector(1, num)];
		end
	end
	numProteins = 0;% Reset number of proteins
	
	if i_printTime == printTime
		% Write out sim_data.csv
		if doDiff == 1
			simDataTable = table(mSizeVector, nSizeVector, addedProteinsVector, 'VariableNames', {'M_Size', 'N_Size', 'Added_Proteins',});
			writetable(simDataTable, 'sim_data.csv');
		else
			simDataTable = table(mSizeVector, nSizeVector, 'VariableNames', {'M_Size', 'N_Size'});
			writetable(simDataTable, 'sim_data.csv');
		end
		
		% Write out div_data.csv
		divDataTable = table(divSizeVector, x0Vector, divTimeVector, divNumVector, 'VariableNames', {'Size_Pre_Division', 'Size_Post_Division', 'Time_Between_Divisions', 'Number_of_Divisions'});
		writetable(divDataTable, 'div_data.csv');
		writematrix(cellAgeMat, 'cell_age.csv');
		writematrix(cellVolMat, 'cell_vol.csv');
		neccesaryValuesTable = table(timeSinceLastDiv, proteinTagInt);
		writetable(neccesaryValuesTable, 'neccesary_values.csv')
		
		% Write out protein_loc.csv
		proteinLocTable = table(iterNumVector, mPosVector, nPosVector, proteinTagOverTimeVector, 'VariableNames', {'Iteration_Number', 'M_Position', 'N_Position', 'Protein_Tag'});
		writetable(proteinLocTable, 'protein_loc.csv');
		
		i_printTime = 0; % Resets iteration number dedicated for data output
	end
end

% Write out sim_data.csv
if doDiff == 1
simDataTable = table(mSizeVector, nSizeVector, addedProteinsVector, 'VariableNames', {'M_Size', 'N_Size', 'Added_Proteins',});
writetable(simDataTable, 'sim_data.csv');
else
simDataTable = table(mSizeVector, nSizeVector, 'VariableNames', {'M_Size', 'N_Size'});
writetable(simDataTable, 'sim_data.csv');
end

% Write out div_data.csv
divDataTable = table(divSizeVector, x0Vector, divTimeVector, divNumVector, 'VariableNames', {'Size_Pre_Division', 'Size_Post_Division', 'Time_Between_Divisions', 'Number_of_Divisions'});
writetable(divDataTable, 'div_data.csv');
cellTable = table(cellAgeMat, cellVolMat, 'VariableNames', {'cell_age', 'cell_volume'});
writetable(cellTable, 'cell_data.csv');
timeSinceLastDivTable = table(timeSinceLastDiv);
writetable(timeSinceLastDivTable, 'time_last_div.csv')

% Write out protein_loc.csv
proteinLocTable = table(iterNumVector, mPosVector, nPosVector, proteinTagOverTimeVector, 'VariableNames', {'Iteration_Number', 'M_Position', 'N_Position', 'Protein_Tag'});
writetable(proteinLocTable, 'protein_loc.csv');
