function [xPosVector, yPosVector] = diffuseVector(xPosVector, yPosVector, matSizeVector, proteinSizeInt)
currentTotalProteins = size(xPosVector,2); % Current total number of proteins
m = matSizeVector(1); % Matrix component that may change
n = matSizeVector(2); % Matrix component that remains unchanged

% If there are proteins, do diffusion
if currentTotalProteins > 0
    diffOrderVector = randperm(currentTotalProteins); % Randomize order of proteins to move
    for i = 1:currentTotalProteins
        probLoc = randperm(100,4) -1; % Randomize 4 values
		newLoc = find(probLoc == max(probLoc)); % Highest of the 4 values is where protein will diffuse to
        
		% If the proteins are larger than 0
		if proteinSizeInt > 0
			% Call diffBoxMake
			diffBoxMat = diffBoxMake(xPosVector, yPosVector, matSizeVector, proteinSizeInt); % Locations where proteins may and may not go based on size of protein
			
			% Check to see if it is valid for the protein to move in that location, if so move there, else check to see if it can move the opposite direction, if so move there, if it can't do either, don't move
			if newLoc == 1 
				if yPosVector(1,diffOrderVector(i)) == m
					if diffBoxMat(yPosVector(1,diffOrderVector(i))-1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) - 1;
					end
				else
					if diffBoxMat(yPosVector(1,diffOrderVector(i))+1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) + 1;
					elseif yPosVector(1,diffOrderVector(i)) ~= 1 && diffBoxMat(yPosVector(1,diffOrderVector(i))-1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) - 1;
					end
				end
			elseif newLoc == 2
				if xPosVector(1,diffOrderVector(i)) == 1
					if diffBoxMat(yPosVector(1,diffOrderVector(i)), n) == 0
						xPosVector(1,diffOrderVector(i)) = n;
					elseif diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))+1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) + 1;
					end
				else
					if diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))-1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) - 1;
					elseif xPosVector(1,diffOrderVector(i)) ~= n && diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))+1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) + 1;
					elseif xPosVector(1,diffOrderVector(i)) == n && diffBoxMat(yPosVector(1,diffOrderVector(i)), 1) == 0
						xPosVector(1,diffOrderVector(i)) = 1;
					end
				end
			elseif newLoc == 3
				if yPosVector(1,diffOrderVector(i)) == 1
					if diffBoxMat(yPosVector(1,diffOrderVector(i))+1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) + 1;
					end
				else
					if diffBoxMat(yPosVector(1,diffOrderVector(i))-1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) - 1;
					elseif yPosVector(1,diffOrderVector(i)) ~= m && diffBoxMat(yPosVector(1,diffOrderVector(i))+1, xPosVector(1,diffOrderVector(i))) == 0
						yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) + 1;
					end
				end
			elseif newLoc == 4
				if xPosVector(1,diffOrderVector(i)) == n
					if diffBoxMat(yPosVector(1,diffOrderVector(i)), 1) == 0
						xPosVector(1,diffOrderVector(i)) = 1;
					elseif diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))-1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) - 1;
					end
				else
					if diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))+1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) + 1;
					elseif xPosVector(1,diffOrderVector(i)) ~= 1 && diffBoxMat(yPosVector(1,diffOrderVector(i)), xPosVector(1,diffOrderVector(i))-1) == 0
						xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) - 1;
					elseif xPosVector(1,diffOrderVector(i)) == 1 && diffBoxMat(yPosVector(1,diffOrderVector(i)), n) == 0
						xPosVector(1,diffOrderVector(i)) = n;
					end
				end

			end
			
		% If the proteins are 0 or lower in size
        else
			
			% Check to see if it is valid for the protein to move in that location, if so move there, else check to see if it can move the opposite direction, if so move there, if it can't do either, don't move
			if newLoc == 1 
				if yPosVector(1,diffOrderVector(i)) == m
					yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) - 1;
				else
					yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) + 1;
				end
			elseif newLoc == 2
				if xPosVector(1,diffOrderVector(i)) == 1
					
					xPosVector(1,diffOrderVector(i)) = n;
				else
					xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) - 1;
				end
			elseif newLoc == 3
				if yPosVector(1,diffOrderVector(i)) == 1
					yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) + 1;
				else
					yPosVector(1,diffOrderVector(i)) = yPosVector(1,diffOrderVector(i)) - 1;
				end
			elseif newLoc == 4
				if xPosVector(1,diffOrderVector(i)) == n
					xPosVector(1,diffOrderVector(i)) = 1;
				else
					xPosVector(1,diffOrderVector(i)) = xPosVector(1,diffOrderVector(i)) + 1;
				end

			end
		end
        
    end
    
end
end