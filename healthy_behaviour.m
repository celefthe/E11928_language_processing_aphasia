%% Behavioural data analysis of the fMRI task

% Total scores, partial scores and reaction times, d prime
% A drift diffusion model is also computed in python using the HDDM library
%
% Sujects had to perform a word-picture verification task, with 1
% congruent condition and 3 incongruent conditions, either phonologically
% related, semantically related, or unrelated, that is for a given image, 
% there are 4 possible associated words anf subjects indicate if it
% mataches or not. Scoring of the test itself is on 30, when all 4 cases of
% an item are correct. The task is adpated directly from Breese and Hillis
% 2003 <https://www.ncbi.nlm.nih.gov/pubmed/15010231> and is sensitive to
% auditory comprehension deficits. The expectations are that pictures
% activates phonological and semantic 'units' (neurons coding the object 
% features) and that in patients, impoverished representation lead to more
% error in matching pitures and incongruent words but still perform fine on 
% congruehnt items while an absence of semantic would create deficit even 
% in the congruent condition (i.e. chance level everywhere).

current = pwd;
if ~strcmp(current(end-3:end),'code')
    error('not running from code directory, CD to the right place then run the code')
end

% get to the bids root dir
cd ..
bids_dir = pwd;

%% data processing for healhty subjects

% get healthy subject derivative directory paths
subj = dir('sub-healthy*');
subjects = {};
for idx = 1:size(subj,1)
     subjects{idx} = [bids_dir filesep 'sourcedata' filesep subj(idx).name];
end

% read the data and compute score, accuracy, d prime and reaction time
N = length(subjects);
Scores = NaN(N,1);
Perf   = NaN(N,4);
RTs    = NaN(N,4);
dprime = NaN(N,3);

for s=1:N
    cd(subjects{s});                     % go the source directory
    name = dir('*.mat');                 % find the mat file
    load(name.name);                     % load the mat file
    conditions = data.conditions1;       % same, phono relatred, sem related, unrelated
    perf       =  data.perf1;            % 1 correct, 0 incorrect
    rt         = data.rt1';              % reaction times
    
    for i=1:4
        index = find(conditions == i);
        accuracy(:,i) = perf(index);     % performance for each condition 
        RT(:,i) = rt(index);             % reaction time for each condition
        if i>1                           % dprime relative to congruent condition
            [~,~,dprime(s,i-1)] = rst_sdt1(sum(conditions == 1),sum(conditions == i),nansum(accuracy(:,1)),length(index)-nansum(accuracy(:,i)),0);
        end
    end
    
    Scores(s) = sum(data.score1);        % total score for the test, already coded in the mat file     
    Perf(s,:) = nanmean(accuracy).*100;  % percentage correct
    RTs(s,:) = nanmedian(RT);            % median reaction time
    
    
    
    clear conditions perf rt accuracy RT
end

% make a figure to show the data (shared in the participants.tsv file)

figure('Name','behavioural results')
subplot(2,9,[1 2 3 4]); [~] = rst_data_plot(Perf,  'estimator','median','kernel','on','newfig','no'); title('Accuracy')
subplot(2,9,5);         [~] = rst_data_plot(Scores,'estimator','median','kernel','on','newfig','no'); title('Scores')
subplot(2,9,[6 7 8 9]); [~] = rst_data_plot(RTs,   'estimator','median','kernel','on','newfig','no'); title('Reaction times')

subplot(2,9,[10:12]);   
[dprimeM,dprimeCI] = rst_data_plot(dprime,'estimator','median','kernel','on','newfig','no'); title('d prime')
subplot(2,9,[13:15]);  
[PerfM,PerfCI]     = rst_data_plot(Perf(:,[2 3 4])-repmat(Perf(:,1),[1 3]),'estimator','median','kernel','on','newfig','no'); title('incongruency scores')
subplot(2,9,[16:18]);   
[RTsM,RTsCI]       = rst_data_plot(RTs(:,[2 3 4]) -repmat(RTs(:,1),[1 3]), 'estimator','median','kernel','on','newfig','no'); title('incongruency reaction times')

% The CI of dprime show that phonological and semantic foils are
% perceptually close to the target (CI do not differ from 0) while the
% unrelated foils are significantly different.

% Incongruency effects are observed for the phonological foils in accuracy 
% and RT, for the semantic foils in RT only, but not for the unrelated
% foils. This suggests that both phonologcal and semantic units of the
% target and foils compete (perceptually close) which leads to increase RT.
% However, decision is only affected by phonological foils (lower accuracy).

% correlation analysis dprime with accuracy/RT in phono and RT in semantic
[r,t,h,outid,hboot,CI] = skipped_correlation([dprime(:,1) dprime(:,[1 2])], [Perf(:,2) RTs(:,[2 3])],1);

% As the perceptual distance of phonolocal foils increase, the accuracy decreases
% As the perceptual distance of semantic foils decrease, the RTs decrease
% and conversely when the perceptual distance increases, the RTs increase
% (true for phonology too, but not significant)


% The drif diffusion model <https://github.com/CPernet/LanguageDecision/blob/master/language_decision/language_decision_model.ipynb>
% shows slower drift rates for all three incongruent foils (post prob 0.7 for unrelated, 0.99 for semantic, and 1 for phonology) 
% which indicates conflict in taking the decision - even with unrelated foils, which is not seens with RT or accuracy.
% The non-decision time ?? 




