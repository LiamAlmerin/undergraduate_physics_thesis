function [yPosVector, matSizeVector, cellVolMat] = growCell(yPosVector, matSizeVector, gridSize, growthRate, cellDiam, cellVolMat, i_cell, i_time, delT)

m = matSizeVector(1); % Matrix component that may change
currentTotalProteins = size(yPosVector,2); % Current total number of proteins

cellVolMat(i_cell, i_time) = cellVolMat(i_cell,i_time-1)*exp(growthRate*delT); % Calculate and record grown cell volume
q = round(cellVolMat(i_cell, i_time)/(pi * ((cellDiam/2)^2)*gridSize), 0); % Calculate new matrix size

if q > m % If new matrix size is bigger than old matrix size
    matSizeVector(1) = q; % Set new matrix size
    if currentTotalProteins > 0 % If there is more than zero proteins, add lines randomly between them
        s = q-m; % Calculate the number of lines to add to matrix
        for i = 1:s % For each new line to add
			if rem(m,2) < 1
				y = rand;
				if y < .5
					r = floor((m + 1) / 2);
				else
					r = ceil((m + 1) / 2);
				end
			else
				r = (m + 1) / 2;
			end
            r = randi(m); % Find a random acceptable location
            for j = 1:currentTotalProteins
				% If the current protein location is greather than the new line that will be inserted, move the protein up one.
                if yPosVector(1,j) >= r
                    yPosVector(1,j) = yPosVector(1,j) +1;
                end
            end
			m = m + 1;
        end
    end
end

end