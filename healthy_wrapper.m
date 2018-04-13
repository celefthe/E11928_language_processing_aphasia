%% Wrapper script to run all analysis for healthy subjects
code_dir = pwd;
if ~strcmp(code_dir(end-3:end),'code')
    error('not running from code directory, CD to the right place then run the code')
end

healthy_first_level;
cd(code_dir);
healthy_second_level;
cd(code_dir);
healthy_behaviour;