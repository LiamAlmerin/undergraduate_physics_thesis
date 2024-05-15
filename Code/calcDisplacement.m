
% Import inital conditions
initialConditionsTable = readtable('initial_conditions.csv');
gridSize = initialConditionsTable.Grid_Size;
delT = initialConditionsTable.Delta_t;
cellDiam = initialConditionsTable.Cell_Diameter;
clear initialConditionsTable

% Import analysis data
analysisTable = readtable('analysisVectors.csv');
timeSecVector = analysisTable.timeSecVector;
timeIterVector = analysisTable.timeIterVector;
numProteins = analysisTable.numProteinsOverTime;
clear analysisTable

iterCount = size(timeIterVector, 1);

divDataTable = readtable('div_data.csv');
numDivVector = divDataTable.Number_of_Divisions;
timeBetweenDivVector = divDataTable.Time_Between_Divisions;
clear divDataTable

previousTime = 0;
timeDivVector = NaN([(size(numDivVector,1) - 1) 1]);
for j = 1:size(numDivVector,1)
	if numDivVector(j) > 0
		timeDivVector(j-1) = (previousTime + timeBetweenDivVector(j));
		previousTime = previousTime + timeBetweenDivVector(j);
	end
end

% Import sim data
simDataTable = readtable('sim_data.csv');
proteinInsert = simDataTable.Added_Proteins;
clear simDataTable

% Import protein locations
proteinLocDataTable = readtable('protein_loc.csv');
iterNumVector = proteinLocDataTable.Iteration_Number;
mPosVector = proteinLocDataTable.M_Position;
nPosVector = proteinLocDataTable.N_Position;
proteinTagVector = proteinLocDataTable.Protein_Tag;
clear proteinLocDataTable

cellRad = cellDiam / 2; % Cell radius in micrometers
theta = gridSize/cellRad; % Grid size of cell in radians
thetaPosVector = nPosVector * theta; % Vector for the theta positions

% Set up arrays for output
nDispMean = NaN([iterCount 1]);
mDispMean = NaN([iterCount 1]);
nDispMeanSqrd = NaN([iterCount 1]);
mDispMeanSqrd = NaN([iterCount 1]);
nDispRtMeanSqrd = NaN([iterCount 1]);
mDispRtMeanSqrd = NaN([iterCount 1]);

nDispMeanMicro = NaN([iterCount 1]);
mDispMeanMicro = NaN([iterCount 1]);
nDispMeanSqrdMicro = NaN([iterCount 1]);
mDispMeanSqrdMicro = NaN([iterCount 1]);
nDispRtMeanSqrdMicro = NaN([iterCount 1]);
mDispRtMeanSqrdMicro = NaN([iterCount 1]);

% Arrays for loop
initPosTheta = [];
initPosM = [];
initTag = [];
curPosTheta = [];
curPosM = [];
curTag = [];

currentDiv = 1;
printTime = 1000;
printIter = 0;


for i = 1:iterCount
	printIter = printIter + 1;

	if numProteins(i) > 0
		% The current positions of the proteins
		curPosTheta = [thetaPosVector(iterNumVector == timeIterVector(i))];
		curPosM = [mPosVector(iterNumVector == timeIterVector(i))];
		curTag = [proteinTagVector(iterNumVector == timeIterVector(i))];
		
		% Appends array with the initial positions of a protein if one is inserted
		if proteinInsert(i) > 0
			if isempty(initTag) || isempty(initTag(initTag == curTag(end)))
				initPosTheta = [initPosTheta ; curPosTheta(end)];
				initPosM = [initPosM ; curPosM(end)];
				initTag = [initTag ; curTag(end)];
			end
		end
		
		% Removes the initial positions of all proteins that are removed
		initPosTheta = initPosTheta(ismember(initTag, curTag));
		initPosM = initPosM(ismember(initTag, curTag));
		initTag = curTag;
		currentDiv = currentDiv + 1;

		% Finds the distance between the current positions and the initial positions
		phi_1 = curPosTheta - initPosTheta;
		phi_2 = 2*pi + curPosTheta - initPosTheta;
		phiVector = [phi_1, phi_2];
		[~, index] = min(abs(phiVector), [], 2, 'linear');
		
		diffPosN = phiVector(index) / theta;
		diffPosM = (curPosM - initPosM);
		
		% Calculates the Mean, Squared Mean, and Root Mean Squared of the displacement
		nDispMean(i) = sum(diffPosN)/numProteins(i);
		mDispMean(i) = sum(diffPosM)/numProteins(i);
		nDispMeanSqrd(i) = sum((diffPosN).^2)/numProteins(i);
		mDispMeanSqrd(i) = sum((diffPosM ).^2)/numProteins(i);
		nDispRtMeanSqrd(i) = sqrt(nDispMeanSqrd(i));
		mDispRtMeanSqrd(i) = sqrt(mDispMeanSqrd(i));
		
		% Distance between the current positions and the initial positions in micrometers
		diffPosMicroN = diffPosN * gridSize;
		diffPosMicroM = diffPosM * gridSize;
		
		% Calculates the Mean, Squared Mean, and Root Mean Squared of the displacement in micrometers
		nDispMeanMicro(i) = sum(diffPosMicroN)/numProteins(i);
		mDispMeanMicro(i) = sum(diffPosMicroM)/numProteins(i);
		nDispMeanSqrdMicro(i) = sum((diffPosMicroN).^2)/numProteins(i);
		mDispMeanSqrdMicro(i) = sum((diffPosMicroM ).^2)/numProteins(i);
		nDispRtMeanSqrdMicro(i) = sqrt(nDispMeanSqrdMicro(i));
		mDispRtMeanSqrdMicro(i) = sqrt(mDispMeanSqrdMicro(i));
	
	else
		% Empties the initial position vectors when there are no proteins currently
		curPosTheta = [];
		curPosM = [];
		curTag = [];
		initPosTheta = [];
		initPosM = [];
		initTag = [];
	end
	
	if printIter == previousTime
		writetable(table(nDispMean, mDispMean, nDispMeanSqrd, mDispMeanSqrd, nDispRtMeanSqrd, mDispRtMeanSqrd, nDispMeanMicro, mDispMeanMicro, nDispMeanSqrdMicro, mDispMeanSqrdMicro, nDispRtMeanSqrdMicro, mDispRtMeanSqrdMicro), 'dispMeanVectors.csv');
		printIter = 0;
	end
end

writetable(table(nDispMean, mDispMean, nDispMeanSqrd, mDispMeanSqrd, nDispRtMeanSqrd, mDispRtMeanSqrd, nDispMeanMicro, mDispMeanMicro, nDispMeanSqrdMicro, mDispMeanSqrdMicro, nDispRtMeanSqrdMicro, mDispRtMeanSqrdMicro), 'dispMeanVectors.csv');