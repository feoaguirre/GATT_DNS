%% Intro
% This is the main file for the DNS
% It will call all the other routines

function [flowHandles, info] = runDNS(caseFile, extraParameters) %#ok<*STOUT>

%% Define case file
if nargin == 0
    caseFile = 'parameters';
end

%% Define if simulation will actually be compiled and run or if just the preprocessing will be done
runSimulation = true;
compileCode = true;
plotDNSDomain = false;

%% Compiling parameters
forceRecompileAll = false;
displayCompiling = false;
optimizeCode = true;
debugger = false;
profiler = false;

%% Set folders for libraries
matlabDir = ''; % Leave empty for automatic directory
decompDir = '/usr/local/2decomp_fft';

%% Data logging
logAll = false; % Save all iterations to log or just when a flow is saved

%% Run parameters file
eval(caseFile);

if nargin == 2
	extraVars = unpackStruct(extraParameters);
    for i = 1:length(extraVars)
        eval([extraVars{i} ' = extraParameters.' extraVars{i} ';'])
    end
end

%% PARALELIZATION PARAMETERS
%% Check if 2D or 3D
tridimensional = mesh.z.n + mesh.z.buffer.i.n + mesh.z.buffer.f.n > 1; %#ok<*NODEF>
if ~tridimensional
    p_col = 1;
		[p_row,~] = get_p_row_col(caseFile,p_max,p_col);
elseif tridimensional AND exist(p_row,'var') AND exist(p_col,'var')
	p_row = p_row;
	p_col =p_col;
else
		[p_row,p_col] = get_p_row_col(caseFile,p_max);
end


%% Add source code path
addpath source
addpath source/boundaries

%% Initialize log file
if ~exist(caseName,'dir')
    mkdir(caseName);
end

fprintf('Parameters file: %s\nCase name: %s\n', caseFile, caseName);

if ~exist([caseName '/log.txt'],'file')
    logFile = fopen([caseName '/log.txt'],'w');
    fprintf(logFile,'Save number\tIteration\tSimulation time\tdt       \tCFL      \tU change\tV change\tW change\tR change\tE change');
    if isfield(mesh,'trackedPoints')
        for i = 1:size(mesh.trackedPoints,1)
            fprintf(logFile,'\tU%d            \tV%d            \tW%d            \tR%d            \tE%d            ',i,i,i,i,i);
        end
    end
    fprintf(logFile,'\n');
    fclose(logFile);
end

%% Run preprocessing routine or reload previous
if ~exist([caseName '/bin'],'dir')
    mkdir([caseName '/bin'])
end

if forceRecompileAll
    compileCode = true;
    delete([caseName '/bin/*.mod']);
	delete([caseName '/bin/*.o']);
end
    
fprintf('Running preprocessor\n')
preprocessing

fprintf('Mesh size: %d x %d x %d\n', mesh.nx, mesh.ny, mesh.nz)

%% Generate initial flow if needed
if genInitialFlow % If there is no previous run, generate new initial flow
    fprintf('Generating new initial flow\n');
     %Compute initial flow
    flow = generateInitialFlow(mesh,flowParameters,flowType.initial,boundary.insideWall);

    % Save initial flow to file
    flowToSave = flow;
    flowToSave.t = 0;
	for var = 'UVWRE'
		flowToSave.(var)(boundary.insideWall) = nan;
	end
    save([caseName '/flow_0000000000.mat'],'-struct','flowToSave','-v7.3');
	
	if isfield(flowType.initial,'meanFile')
		flowTypeTemp.initial.type = 'file';
		flowTypeTemp.initial.flowFile = flowType.initial.meanFile;
		if isfield(flowType.initial,'meshFile')
			flowTypeTemp.initial.meshFile = flowType.initial.meshFile;
        end
		
		meanFlow = generateInitialFlow(mesh,flowParameters,flowTypeTemp.initial,boundary.insideWall);

		% Save initial flow to file
		flowToSave = meanFlow;
		flowToSave.t = 0;
		for var = 'UVWRE'
			flowToSave.(var)(boundary.insideWall) = nan;
		end
		save([caseName '/meanflowSFD.mat'],'-struct','flowToSave','-v7.3');
	end
	

else
    fprintf('Resuming from file number %d\n',time.nStep)
end

% Copy parameters file to Fortran folder and write to log2.txt
logFile2 = fopen([caseName '/bin/log2.txt'],'a');
fprintf(logFile2,'DNS started at %s\n', datestr(now,'dd-mmm-yyyy HH:MM:SS'));
fprintf(logFile2,'Parameters file: %s\n', [caseFile '.m']);
fprintf(logFile2,'Starting flow file: flow_%.10d.mat\n\n', time.nStep);
if exist([caseName '/bin/parameters.m'],'file')
	[parametersDiffStatus,parametersDiff] = system(['diff ' caseName '/bin/parameters.m ' caseFile '.m']);
	if parametersDiffStatus
		fprintf(logFile2,'Parameters file was changed:\n%s\n', parametersDiff);
	end
end
fclose(logFile2);

copyfile([caseFile '.m'],[caseName '/bin/parameters.m']);

if nargin == 2
	save([caseName '/bin/extraParameters.mat'], '-struct', 'extraParameters')
end

% Compile code
if compileCode
	fprintf('Compiling code\n')
	compileFortran
end

%% Plot domain
if plotDNSDomain
    plotDomain %#ok<*UNRCH>
    drawnow
end

%% Call Fortran code
if runSimulation && ~debugger && ~profiler
    fprintf('Starting code\n')
    tic
    system(['cd ' caseName '/bin && mpirun -np ' num2str(p_row*p_col) ' main ' caseName]);
    toc
elseif runSimulation && debugger
    fprintf('Starting code with debugger\n')
    system(['cd ' caseName '/bin && mpirun -n ' num2str(p_row*p_col) ' xterm -sl 1000000 -fg white -bg black -hold -e gdb -ex run --args ./main ' caseName]);
elseif runSimulation && profiler
    fprintf('Starting code with profiler\n')
    system('export GMON_OUT_PREFIX=''gmon.out''');
    tic
    system(['cd ' caseName '/bin && mpirun -np ' num2str(p_row*p_col) ' main ' caseName]);
    toc
    system(['cd ' caseName '/bin && gprof -l main gmon.out > profile.txt']);
    system(['mv ' caseName '/bin/profile.txt .']);
end

%% Write to log2 file
logFile2 = fopen([caseName '/bin/log2.txt'],'a');
fprintf(logFile2,'DNS finished at %s\n\n', datestr(now,'dd-mmm-yyyy HH:MM:SS'));
fclose(logFile2);

%% Get outputs if needed
if nargout > 0
    allCaseFiles = dir(caseName);
    flowHandles = {};
    for i = 1:length(allCaseFiles)
        name = allCaseFiles(i).name;
        if length(name) == 19 && ~isempty(regexp(name,'flow_\d*.mat','once'))
            flowHandles{end+1} = matfile([caseName '/' name]); %#ok<AGROW>
        end
    end
end
if nargout == 2
    allVars = whos;
    for i = 1:length(allVars)
        eval(['info.' allVars(i).name ' = ' allVars(i).name ';']);
    end
end
end


function varList = unpackStruct(structure)
    varList = {};
    varNames = fieldnames(structure);
    for i = 1:length(varNames)
        varName = varNames{i};
        if isstruct(structure.(varName))
            varList2 = unpackStruct(structure.(varName));
            for j = 1:length(varList2)
                varList{end+1} = [varName '.' varList2{j}];
            end
        else
            varList{end+1} = varName;
        end
    end
end

function [p_row,pcow] = get_p_row_col(baseFile,p_max,varArgin)

    addpath source
    addpath source/boundaries
    
		if exist(varArgin(1),'var')
    
	    eval(baseFile)
	    p_row = p_max;
	    p_col = 1;



	    while p_col > 0
	        try
	            p_row
	            meshAddFixedPoints
	
	            % Run mesh generator
	            [mesh.X, mesh.x, mesh.nx] = generateMesh(domain.xi,domain.xf,mesh.x,'X');
	            [mesh.Y, mesh.y, mesh.ny] = generateMesh(domain.yi,domain.yf,mesh.y,'Y');
	            [mesh.Z, mesh.z, mesh.nz] = generateMesh(domain.zi,domain.zf,mesh.z,'Z');
	
	            %% Select boundary conditions
	            [boundary,mesh] = getBoundaryConditions(flowType,mesh,flowParameters,[numMethods.neumannOrder numMethods.neumann2Order]);
	
	            domainSlicesY = getDomainSlices(mesh.ny,p_row);
	            domainSlicesZ = getDomainSlices(mesh.nz,p_col);
	
	            [~] = initBoundaries(boundary,mesh,domainSlicesY,domainSlicesZ,p_row,p_col);
	            
	            return
	        catch
	            p_row = p_row - 1;
	        end
	    end

			else.  %% need to check if everything works correctly

			eval(baseFile)


maxProduct = 0;
    p_row = 1;
    p_col = p_max;
    for i = p_max:-1:1
        if mod(N, i) == 0
            j = N / i;
            if i * j <= N && i * j > maxProduct
                p_row = i;
                p_col = j;
                
							try
	            p_row
	            meshAddFixedPoints
	
	            % Run mesh generator
	            [mesh.X, mesh.x, mesh.nx] = generateMesh(domain.xi,domain.xf,mesh.x,'X');
	            [mesh.Y, mesh.y, mesh.ny] = generateMesh(domain.yi,domain.yf,mesh.y,'Y');
	            [mesh.Z, mesh.z, mesh.nz] = generateMesh(domain.zi,domain.zf,mesh.z,'Z');
	
	            %% Select boundary conditions
	            [boundary,mesh] = getBoundaryConditions(flowType,mesh,flowParameters,[numMethods.neumannOrder numMethods.neumann2Order]);
	
	            domainSlicesY = getDomainSlices(mesh.ny,p_row);
	            domainSlicesZ = getDomainSlices(mesh.nz,p_col);
	
	            [~] = initBoundaries(boundary,mesh,domainSlicesY,domainSlicesZ,p_row,p_col);
	            
	            return
				 end


            end
        end
    end



	        
end
