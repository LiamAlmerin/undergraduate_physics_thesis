function [xPosVector, yPosVector, addedProteins, proteinTagInt, proteinTagVector] = insertProtein(xPosVector, yPosVector, matSizeVector, proteinSizeInt, addedProteins, proteinChance, proteinTagInt, proteinTagVector)

m = matSizeVector(1); % Matrix component that may change
n = matSizeVector(2); % Matrix component that remains unchanged
y = rand; % Number to compare proteinChance against

notFit = 1; % Legacy, comment out first and check if no crash before deleting
proteinFitInt = 0; % Legacy, comment out first and check if no crash before deleting

currentTotalProteins = size(xPosVector,2); % Current total number of proteins

% Calls createBoxMat 
boxSizeMat = createBoxMat(xPosVector, yPosVector, matSizeVector, proteinSizeInt);

if y < proteinChance % Check to see if protein is to be inserted
	
	% Find all acceptable insertion points
	[i, j] =  find(boxSizeMat == 0); 
	
	e = size(i,1); % Number of possible locations
	if e > 0 % If possible locations are greater than 0, insert protein
		s = randi(e); % Choose random location
		k = j(s); % X location to insert
		l = i(s); % Y location to insert
		
		% Insert protein
		xPosVector = [xPosVector k];
		yPosVector = [yPosVector l];
		proteinTagVector = [proteinTagVector proteinTagInt];
		proteinTagInt = proteinTagInt + 1;
		addedProteins = addedProteins + 1; % Increase number of added proteins
		
		% Calls createBoxMat 
		boxSizeMat = createBoxMat(xPosVector, yPosVector, matSizeVector, proteinSizeInt);
	end
end
end
