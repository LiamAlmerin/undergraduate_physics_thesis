function [xAvg, yAvg, xSqAvg, ySqAvg, xDif, yDif] = calcAvgDiff(matSizeVector, xPosVector, yPosVector, gridSize, diffusionCoefficient, i_time)
m = matSizeVector(1); % Matrix component that may change
n = matSizeVector(2); % Matrix component that remains unchanged
delT = (gridSize)^2 /(2 * diffusionCoefficient); % Time in seconds between each iteration
currentTotalProteins = size(xPosVector,2); % Current total number of proteins

% Calculate middle of matrix
mMid = (1 + m)/2;
nMid = (1 + n)/2;

% Set up variables
xNum = 0; % Setup for the added locations of the x positions
xSqNum = 0; % Setup for the added locations of the squared x positions
yNum = 0; % Setup for the added locations of the y positions
ySqNum = 0; % Setup for the added locations of the squared y positions

for i = 1:currentTotalProteins
    yNum = yNum + ((yPosVector(1,i) - mMid) * gridSize); % Added locations of the y positions
	ySqNum = ySqNum + (((yPosVector(1,i) - mMid) * gridSize)^2); % Added locations of the squared y positions
    xNum = xNum + ((xPosVector(1,i) - nMid) * gridSize); % Added locations of the x positions
	xSqNum = xSqNum + (((xPosVector(1,i) - nMid) * gridSize)^2); % Added locations of the squared x positions
end

% Calculated results
xAvg = xNum/currentTotalProteins; % Calculated average x value of proteins
xSqAvg = xSqNum/currentTotalProteins; % Calculated average x squared value of proteins

yAvg = yNum/currentTotalProteins; % Calculated average y value of proteins
ySqAvg = ySqNum/currentTotalProteins; % Calculated average y squared value of proteins

xDif = xSqAvg/(2*(delT*i_time)); % Calculated x diffusion value
yDif = ySqAvg/(2*(delT*i_time)); % Calculated y diffusion value
end