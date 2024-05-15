function [xPosVector, yPosVector, matSizeVector, cellAgeMat, cellVolMat, timeSinceLastDiv, divTimeVector, x0Vector, divSizeVector, mi, i_cell, divNum, proteinTagVector] = divCell(xPosVector, yPosVector, matSizeVector, cellAgeMat, cellVolMat, mi, x0Vector, divSizeVector, timeSinceLastDiv, divTimeVector, i_cell, i_time, delT, divNum, gridSize, hillCoefficient, growthRate, cellDiam, proteinTagVector)
currentTotalProteins = size(xPosVector,2); % Current total number of proteins
nh = hillCoefficient;
x = cellVolMat(i_cell, i_time)/(pi * ((cellDiam/2)^2));  % Current cell length (um) 
x0 = mi* gridSize; % Initial cell length (um) (Cell length after division)
age = cellAgeMat(i_cell, i_time)/60 *delT; % Current cell age at current time in minutes

% Equations to calculate the probablity to divide
hm = 7.7/(1+age/(47)); % Units must be micrometers
kmax = .36*age^2/((18^2+age^2)*60); % Units are per second
x_prev = cellVolMat(i_cell, i_time-1)/(pi * ((cellDiam/2)^2)); % Previous cell length (um) 
F = 1 - ((hm^nh+x0^nh)/(hm^nh+x^nh))^(kmax/(growthRate*nh));
Fp = 1 - ((hm^nh+x0^nh)/(hm^nh+x_prev^nh))^(kmax/(growthRate*nh));
prob_div = (F-Fp)/(1-Fp); % Probability to divide

z = rand; % Number to compare prob_div against
m = matSizeVector(1); % Matrix component that may change

if z <= prob_div % Check to see if cell should be divided
	mc = round(m/2, 0); % Find center of cell to cut
	mi = mc; % Set inital size of cell after cutting
	x0Vector = [x0Vector ; mi*gridSize]; % Calculate and record size of cell after division
	divSizeVector = [divSizeVector ; m * gridSize]; % Calculate the size of the cell before division

    xPosVectorNew = []; % X position of protein, related to the n value of the matrix
    yPosVectorNew = []; % Y position of protein, related to the m value of the matrix
    proteinTagVectorNew = [];
	
	% Remove proteins that extend the bounds of the new cell
	for i=1:currentTotalProteins
        if yPosVector(1,i) <= mc
            xPosVectorNew = [xPosVectorNew xPosVector(1,i)];
            yPosVectorNew = [yPosVectorNew yPosVector(1,i)];
			proteinTagVectorNew = [proteinTagVectorNew proteinTagVector(1,i)];
        end
    end
    xPosVector = xPosVectorNew;
    yPosVector = yPosVectorNew;
	proteinTagVector = proteinTagVectorNew;
	
	i_cell = i_cell + 1; % Increment cell number
	volume0 = pi * ((cellDiam/2)^2)*(gridSize * mc); % Cell volume in micrometers^3 after division
	cellVolMat(i_cell,i_time) = volume0; % Record cell volume after division
	cellAgeMat(i_cell, i_time) = 0; % Set cell age back to 0
	
	% Record the time since last division
	divTimeVector = [divTimeVector ; timeSinceLastDiv];
	
	divNum = divNum + 1; % Increment the number of divisions occured
	timeSinceLastDiv = 0; % Reset time since last division
    m = mc; % Set matrix component that may change size to new size of cell
else 
end
matSizeVector(1) = m; % Set matrix to the size of the cell, either after cutting or keep the current size if not cut
end