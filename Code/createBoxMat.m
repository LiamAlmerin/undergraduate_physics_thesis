function [boxSizeMat] = createBoxMat(xPosVector, yPosVector, matSizeVector, proteinSizeInt)

m = matSizeVector(1); % Matrix component that may change
n = matSizeVector(2); % Matrix component that remains unchanged
currentTotalProteins = size(xPosVector,2); % Current total number of proteins
xCoordVector = []; % X positions that proteins can't go
yCoordVector = []; % Y position that proteins can't go

if proteinSizeInt > 0 %Check if size is greater than 0
	% Create a barrier calculated from center of hypothetical protein
	barrierVector = (1:1:proteinSizeInt) - 1; % 
	for b = 1:barrierVector(end)
		barrierVector = [barrierVector -(barrierVector(b + 1))];
	end
	
	% Locate the space each protein takes up 
	for l = 1:currentTotalProteins
		for j = 1:size(barrierVector,2)
			for i = 1:size(barrierVector,2)
				% Makes sure to respect boundary conditions
				if (xPosVector(1,l) + barrierVector(i)) >= 1 && (xPosVector(1,l) + barrierVector(i)) <= n && (yPosVector(1,l) + barrierVector(j)) >= 1 && (yPosVector(1,l) + barrierVector(j)) <= m
					xCoordVector = [xCoordVector (xPosVector(1,l) + barrierVector(i))];
					yCoordVector = [yCoordVector (yPosVector(1,l) + barrierVector(j))];
				elseif (xPosVector(1,l) + barrierVector(i)) < 1 && (yPosVector(1,l) + barrierVector(j)) >= 1 && (yPosVector(1,l) + barrierVector(j)) <= m
					xCoordVector = [xCoordVector n+(xPosVector(1,l) + barrierVector(i))];
					yCoordVector = [yCoordVector (yPosVector(1,l) + barrierVector(j))];
				elseif (xPosVector(1,l) + barrierVector(i)) > n && (yPosVector(1,l) + barrierVector(j)) >= 1 && (yPosVector(1,l) + barrierVector(j)) <= m
					xCoordVector = [xCoordVector (xPosVector(1,l) + barrierVector(i))-n];
					yCoordVector = [yCoordVector (yPosVector(1,l) + barrierVector(j))];
				end
			end
		end
	end
	
	% Prevent proteins from being too close to n (x) edge, retired
	%for g = 1:n
	%	for h = 1:(proteinSizeInt-1)
	%		if mod(h, 2) ~= 0
	%			xCoordVector = [xCoordVector g];
	%			yCoordVector = [yCoordVector ceil(h/2)];
	%		else
	%			xCoordVector = [xCoordVector g];
	%			yCoordVector = [yCoordVector (m+1-ceil((h-1)/2))];
	%		end
	%	end
	%end
	
	boxSizeMat = zeros(m,n); % Set up boxSizeMat
	boxSizeMat(sub2ind(size(boxSizeMat), yCoordVector, xCoordVector)) = 1; % Sets locations that proteins can't be to 1
    
else
    boxSizeMat = zeros(m,n); % If proteins have no size
end
end