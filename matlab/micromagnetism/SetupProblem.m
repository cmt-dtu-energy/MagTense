function [ProblemSetupStruct, SimulationFileName] = SetupProblem(MySim)

%% Create a directory for the simulation
k = 1;
DirectoryFilename = ['Simulations\' char(datetime('today','Format','yyyy-MM-dd')) '_' MySim.SimulationName '_' num2str(k,'%04.0f')];
if exist(DirectoryFilename,'dir')
    DirectoryFilename_temp = DirectoryFilename;
    while (exist(DirectoryFilename_temp,'dir'))
        k = k+1;
        DirectoryFilename_temp = ['Simulations\' char(datetime('today','Format','yyyy-MM-dd')) '_' MySim.SimulationName '_' num2str(k,'%04.0f')];
    end
    DirectoryFilename = DirectoryFilename_temp;
end


MySim.FileName = [DirectoryFilename '\' MySim.SimulationName];

%% Call the default problem parameters
ProblemSetupStruct = DefaultProblemParameters();

if ProblemSetupStruct.SaveTheResult
    mkdir([DirectoryFilename]);
end

%% Evaluate all variables in the user provided struct MySim and put them in ProblemSetupStruct
user_defined_vars = fieldnames(MySim) ;
for i=1:numel(user_defined_vars)
    ProblemSetupStruct.(user_defined_vars{i}) = MySim.(user_defined_vars{i}) ;
end

%% Save the data
if ProblemSetupStruct.SaveTheResult
    save(MySim.FileName,'ProblemSetupStruct');
end
SimulationFileName = MySim.FileName;
end