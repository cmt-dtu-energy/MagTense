function [ProblemSetupStruct, SimulationFileName] = SetupProblem(MySim)

SimulationFileName = '';

%% Call the default problem parameters
ProblemSetupStruct = DefaultProblemParameters();

%% Evaluate all variables in the user provided struct MySim and put them in ProblemSetupStruct
user_defined_vars = fieldnames(MySim) ;
for i=1:numel(user_defined_vars)
    ProblemSetupStruct.(user_defined_vars{i}) = MySim.(user_defined_vars{i}) ;
end

%% Create a directory for the simulation and save the data
if ProblemSetupStruct.SaveTheResult
    k = 1;
    if ~isfield(MySim,'DirectoryFilename')
        DirectoryFilename = ['Simulations\' char(datetime('today','Format','yyyy-MM-dd')) '_' MySim.SimulationName '_' num2str(k,'%04.0f')];
        if exist(DirectoryFilename,'dir')
            DirectoryFilename_temp = DirectoryFilename;
            while (exist(DirectoryFilename_temp,'dir'))
                k = k+1;
                DirectoryFilename_temp = ['Simulations\' char(datetime('today','Format','yyyy-MM-dd')) '_' MySim.SimulationName '_' num2str(k,'%04.0f')];
            end
            DirectoryFilename = DirectoryFilename_temp;
        end
    end

    mkdir([DirectoryFilename]);
    MySim.FileName = [DirectoryFilename '\' MySim.SimulationName];

    save(MySim.FileName,'ProblemSetupStruct');
    
    SimulationFileName = MySim.FileName;
end

end