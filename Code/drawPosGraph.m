
iterCount = 500001; % Full length
%iterCount = 103359; % This corresponds to 10000 seconds

myColorMap = [0 0 0; 0 0.2 1; 1 1 0; 1 0 0];

% Import Variables
initialTable = readtable('initial_conditions.csv');
proteinSizeInt = ceil((initialTable.Protein_Radius_Overlapping + 1) /2);

clear initialTable

simDataTable = readtable('sim_data.csv');
mSizeVec = simDataTable.M_Size;
nSizeVec = simDataTable.N_Size;
mMax = max(mSizeVec);

clear simDataTable

analysisTable = readtable('analysisVectors.csv');
numProteinsVec = analysisTable.numProteinsOverTime;
timeSecVec = analysisTable.timeSecVector;

clear analysisTable

% Import protein locations
proteinLocDataTable = readtable('protein_loc.csv');
% Set up arrays
iterNumVector = proteinLocDataTable.Iteration_Number;
mPosVector = proteinLocDataTable.M_Position;
nPosVector = proteinLocDataTable.N_Position;

clear proteinLocDataTable

parfor i = 1:iterCount

	cellMatrix = -1.*ones(nSizeVec(i)+1, mMax+1); % Sets up matrix to be a consistant size to -1
	
	cellMatrix(1:nSizeVec(i), 1:mSizeVec(i)) = 0; % Sets positions where the cell exists to 0
	
	cellMatrix(nSizeVec(i)+1, mMax+1) = 1;
	cellMatrix(nSizeVec(i), mMax+1) = 2; 	% Values entered to force consistant colors
	
	% Sets positions where the proteins exist to 1 
	if numProteinsVec(i) > 0
		nIterPosVec = nPosVector(iterNumVector == i-1).';
		mIterPosVec = mPosVector(iterNumVector == i-1).';
		boxMat = createBoxMat(nIterPosVec, mIterPosVec, [mSizeVec(i), nSizeVec(i)], proteinSizeInt);
		outerMat = diffBoxMake(nIterPosVec, mIterPosVec, [mSizeVec(i), nSizeVec(i)], proteinSizeInt);
		posMat = padarray(boxMat.' + outerMat.', [1,mMax - mSizeVec(i) + 1], 0, 'post');
		cellMatrix = cellMatrix + posMat;
	end
	
	% Makes and prints out figure
	f = figure('Name', 'Image' , 'Position', [0, 0, 1920, 1080],'NumberTitle','off', 'WindowState', 'minimized');
	
	pcolor(cellMatrix);
	colormap(myColorMap);
	axis ij
	axis image
	
	titleStr = sprintf('Protein Assembly Positions at Time %8.2f Seconds', timeSecVec(i));
	title(titleStr);
	set(gca, 'Box', 'off');
	set(gca, 'TickDir', 'out');
	set(gca, 'FontSize', 15);
	f.PaperSize = [16 9];

	exportName = sprintf('anim/proteinAnim_%06d', (i-1));
	print(exportName, '-dpng', '-r0')
	% print(exportName, '-dpdf', '-bestfit')
	close all
	
	
end