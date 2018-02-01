function read_fMRI_logfile(varargin)

current = pwd;
% input check
if nargin == 0
    [f,p] = uigetfile('*.mat','select file');
    cd(p); load(f); cd(current)
elseif nargin == 1
    data = varargin{1};
else
    error('to many arguments')    
end

% which session
session = input('which session? ');
if session == 1
    onset = data.onsets1;
    duration =  data.durations1;
    conditions = data.conditions1;
    perf =  data.perf1;
    percent = data.perf_per_condition1';
    rt = data.rt1';
elseif session == 2
    onset = data.onsets2;
    duration =  data.durations2;
    conditions = data.conditions2;
    perf =  data.perf2;
    percent = data.perf_per_condition2';
    rt = data.rt2';
end

% outputs for SPM
for i= 1:4
    index = find(conditions == i);
    condition_onset{i} = onset(index);
    condition_duration{i} = duration(index);
    test = perf(index);
    if sum(test) == 30
        errors_onset{i} = [];
        errors_duration{i} = [];
    else
        test  = single(test == 0) + single(isnan(test));
        errors_onset{i} = condition_onset{i}(find(onset(index).*test));
        errors_duration{i} =  condition_duration{i}(find(duration(index).*test));
    end
end
d = uigetdir(current,'save spm info in?'); cd(d);
save condition_onset condition_onset
save condition_duration condition_duration
save errors_onset errors_onset
save errors_duration  errors_duration

% behavioural outputs
cd (current)
accuracy = zeros(1,4);
RT = zeros(1,4);
for i=1:4
     index = find(conditions == i);
     accuracy(:,i) = nanmean(percent(index));
     RT(:,i) = nanmean(rt(index));
end
d = uigetdir(current,'save behavioural info in?'); cd(d);
save accuracy.txt accuracy -ascii
save RT.txt RT -ascii

