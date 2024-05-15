function [diffBoxMat] = diffBoxMake(xPosVector, yPosVector, matSizeVector, proteinSizeInt)
m = matSizeVector(1); % Matrix component that may change
n = matSizeVector(2); % Matrix component that remains unchanged
currentTotalProteins = size(xPosVector,2); % Current total number of proteins
xCoordVector = []; % X positions that proteins can't go
yCoordVector = [];

if proteinSizeInt > 0 %Check if size is greater than 0
    if proteinSizeInt == 1 % The case where the protein size is 1
		% Create a barrier calculated from center of hypothetical protein
		barrierVector = (1:1:proteinSizeInt) - 1;
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
    else % The case where the protein size is not 1
		% Locate the places where a protein cannont move to, while making sure each protein doesn't inhibit itself
        outsideVector = (1:1:proteinSizeInt) - 1; % Size of the barrier
        insideVector = outsideVector(1:(size(outsideVector,2)-1)); % Size of the interior area to exclude
        maxVector = [-outsideVector(end),outsideVector(end)]; % Maximum size of the barrier for a hypothetical protein
        
		% Create a barrier calculated from center of hypothetical protein
		for b = 1:outsideVector(end) 
            outsideVector = [outsideVector -(outsideVector(b + 1))];
        end
		
		% Create the interior area to exclude for a hypothetical protein
        for b = 1:insideVector(end)
            insideVector = [insideVector -(insideVector(b + 1))];
        end
		
		% Calculate the location of the barrier
        for l = 1:currentTotalProteins
            for j = 1:size(maxVector,2)
                for i = 1:size(outsideVector,2)
					% Makes sure to respect boundary conditions
                    if (xPosVector(1,l) + outsideVector(i)) >= 1 && (xPosVector(1,l) + outsideVector(i)) <= n && (yPosVector(1,l) + maxVector(j)) >= 1 && (yPosVector(1,l) + maxVector(j)) <= m
                        xCoordVector = [xCoordVector (xPosVector(1,l) + outsideVector(i))];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + maxVector(j))];
                    elseif (xPosVector(1,l) + outsideVector(i)) < 1 && (yPosVector(1,l) + maxVector(j)) >= 1 && (yPosVector(1,l) + maxVector(j)) <= m
                        xCoordVector = [xCoordVector n+(xPosVector(1,l) + outsideVector(i))];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + maxVector(j))];
                    elseif (xPosVector(1,l) + outsideVector(i)) > n && (yPosVector(1,l) + maxVector(j)) >= 1 && (yPosVector(1,l) + maxVector(j)) <= m
                        xCoordVector = [xCoordVector (xPosVector(1,l) + outsideVector(i))-n];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + maxVector(j))];
                    end
                end
            end
        end
		
		% Calculate the location of the interior of the barrier
        for l = 1:currentTotalProteins
            for j = 1:size(insideVector,2)
                for i = 1:size(maxVector,2)
					% Makes sure to respect boundary conditions
                    if (xPosVector(1,l) + maxVector(i)) >= 1 && (xPosVector(1,l) + maxVector(i)) <= n && (yPosVector(1,l) + insideVector(j)) >= 1 && (yPosVector(1,l) + insideVector(j)) <= m
                        xCoordVector = [xCoordVector (xPosVector(1,l) + maxVector(i))];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + insideVector(j))];
                    elseif (xPosVector(1,l) + maxVector(i)) < 1 && (yPosVector(1,l) + insideVector(j)) >= 1 && (yPosVector(1,l) + insideVector(j)) <= m
                        xCoordVector = [xCoordVector n+(xPosVector(1,l) + maxVector(i))];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + insideVector(j))];
                    elseif (xPosVector(1,l) + maxVector(i)) > n && (yPosVector(1,l) + insideVector(j)) >= 1 && (yPosVector(1,l) + insideVector(j)) <= m
                        xCoordVector = [xCoordVector (xPosVector(1,l) + maxVector(i))-n];
                        yCoordVector = [yCoordVector (yPosVector(1,l) + insideVector(j))];
                    end
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
	
	diffBoxMat = zeros(m,n); % Set up diffBoxMat
	diffBoxMat(sub2ind(size(diffBoxMat), yCoordVector, xCoordVector)) = 1; % Sets locations that proteins can't be to 1
else
    diffBoxMat = zeros(m,n); % If proteins have no size
end
end