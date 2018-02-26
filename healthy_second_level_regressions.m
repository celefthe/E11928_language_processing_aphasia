
%% fMRI second level regressions for healthy subjects only

% clear memoery, load SPM defaults and job manager
clear variables
spm_jobman('initcfg');
spmroot = fileparts(which('spm'));
spm('Defaults','fmri')

% the scripts assumes we run from inside the code directory
current = pwd;
if ~strcmp(current(end-3:end),'code')
    error('not running from code directory, CD to the right place then run the code')
end

% get to the bids root dir
cd ..
bids_dir = pwd;

% get d_primes csv
dprimes = readtable([bids_dir filesep 'derivatives' filesep 'd_primes' filesep 'healthy_dprimes.csv']);
drift_rates = readtable([bids_dir filesep 'derivatives' filesep 'ddm' filesep 'healthy_driftrates.csv']);

% get healthy subject derivative directory paths
subj = dir('sub-healthy*');

% for each subject, get the boosted incongruence contrasts, boosted betas, as well as d' and ddm
% parameters
subjects = {};
for idx = 1:size(subj,1)
     subjects{idx}.path = [bids_dir filesep 'derivatives' filesep subj(idx).name];
     
     subjects{idx}.con.phonology = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0007.nii,1'];
     subjects{idx}.con.semantic = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0008.nii,1'];
     subjects{idx}.con.unrelated = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0009.nii,1'];
     subjects{idx}.dprime.phonology = dprimes.SS_CP(idx);
     subjects{idx}.dprime.semantic = dprimes.SS_CS(idx);
     subjects{idx}.dprime.unrelated = dprimes.SS_US(idx);
     
     subjects{idx}.beta.ss = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0001.nii,1']
     subjects{idx}.beta.cp = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0005.nii,1']
     subjects{idx}.beta.cs = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0009.nii,1']
     subjects{idx}.beta.us = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0013.nii,1']
     subjects{idx}.subj_idx = drift_rates.subj_idx(idx);
     subjects{idx}.drift.ss = drift_rates.SS(idx);
     subjects{idx}.drift.cp = drift_rates.CP(idx);
     subjects{idx}.drift.cs = drift_rates.CS(idx);
     subjects{idx}.drift.us = drift_rates.US(idx);
end

%% regressor batch

% set up directory tree
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    {[bids_dir filesep 'derivatives' filesep 'group' filesep 'healthy']};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'regressions';

%% regression with d'
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''regressions''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'd_prime';

% multiple regression with d'
matlabbatch{3}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''d_prime''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% get the boosted contrasts for phono, semantic and unrelated
matlabbatch{3}.spm.stats.factorial_design.des.mreg.scans = {}
for idx=1:size(subjects,2)
    matlabbatch{3}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.con.phonology};
    matlabbatch{3}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.con.semantic};
    matlabbatch{3}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.con.unrelated};
end

% get the matching d'
% for each patient, these are [d';0;0] for phono, [0;d';0] for sem, and [0;0;d']
% for unrelated to match the boosted contrast image order
phono_dprimes = [];
for idx = 1:size(subjects,2)
    subj_dp = [subjects{idx}.dprime.phonology; 0; 0];
    phono_dprimes = [phono_dprimes ; subj_dp];
end
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(1).c = phono_dprimes;
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'phonology';
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 1;

sem_dprimes = [];
for idx = 1:size(subjects,2)
    subj_dp = [0; subjects{idx}.dprime.semantic; 0];
    sem_dprimes = [sem_dprimes ; subj_dp];
end
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(2).c = sem_dprimes;
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'semantic';
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 1;

unr_dprimes = [];
for idx = 1:size(subjects,2)
    subj_dp = [0; 0; subjects{idx}.dprime.unrelated];
    unr_dprimes = [unr_dprimes ; subj_dp];
end
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(3).c = unr_dprimes;
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(3).cname = 'unrelated';
matlabbatch{3}.spm.stats.factorial_design.des.mreg.mcov(3).iCC = 1;

matlabbatch{3}.spm.stats.factorial_design.des.mreg.incint = 1;

matlabbatch{3}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{3}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{3}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{3}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{3}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{3}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{3}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{4}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{4}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{4}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{5}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'phonology';
matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0];
matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{2}.tcon.name = 'semantic';
matlabbatch{5}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 0];
matlabbatch{5}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{3}.tcon.name = 'unrelated';
matlabbatch{5}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 1];
matlabbatch{5}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{4}.tcon.name = 'con_sum';
matlabbatch{5}.spm.stats.con.consess{4}.tcon.weights = [0 1 1 1];
matlabbatch{5}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.consess{5}.fcon.name = 'incongruence';
matlabbatch{5}.spm.stats.con.consess{5}.fcon.weights = [zeros(3,1) eye(3)];
matlabbatch{5}.spm.stats.con.consess{5}.fcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.delete = 1;


%% regression with drift rate per condition
matlabbatch{6}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''regressions''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{6}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'drift_rate_per_condition';


matlabbatch{7}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''drift_rate_per_condition''', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{7}.spm.stats.factorial_design.des.mreg.scans = {};
for idx=1:size(subjects,2)
    matlabbatch{7}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.ss};
    matlabbatch{7}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.cp};
    matlabbatch{7}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.cs};
    matlabbatch{7}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.us};
end

matlabbatch{7}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});

drift_ss = [];
for idx = 1:size(subjects,2)
    subj_v = [subjects{idx}.drift.ss; 0; 0; 0];
    drift_ss = [drift_ss ; subj_v];
end
matlabbatch{7}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{7}.spm.stats.factorial_design.cov(1).c = drift_ss;
matlabbatch{7}.spm.stats.factorial_design.cov(1).cname = 'ss';
matlabbatch{7}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{7}.spm.stats.factorial_design.cov(1).iCC = 1;

drift_cp = [];
for idx = 1:size(subjects,2)
    subj_v = [0; subjects{idx}.drift.cp; 0; 0];
    drift_cp = [drift_cp ; subj_v];
end
matlabbatch{7}.spm.stats.factorial_design.cov(2).c = drift_cp;
matlabbatch{7}.spm.stats.factorial_design.cov(2).cname = 'cp';
matlabbatch{7}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{7}.spm.stats.factorial_design.cov(2).iCC = 1;

drift_cs = [];
for idx = 1:size(subjects,2)
    subj_v = [0; 0; subjects{idx}.drift.cs; 0];
    drift_cs = [drift_cs ; subj_v];
end
matlabbatch{7}.spm.stats.factorial_design.cov(3).c = drift_cs;
matlabbatch{7}.spm.stats.factorial_design.cov(3).cname = 'cs';
matlabbatch{7}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{7}.spm.stats.factorial_design.cov(3).iCC = 1;

drift_us = [];
for idx = 1:size(subjects,2)
    subj_v = [0; 0; 0; subjects{idx}.drift.us];
    drift_us = [drift_us ; subj_v];
end
matlabbatch{7}.spm.stats.factorial_design.cov(4).c = drift_us;
matlabbatch{7}.spm.stats.factorial_design.cov(4).cname = 'us';
matlabbatch{7}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{7}.spm.stats.factorial_design.cov(4).iCC = 1;
matlabbatch{7}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{7}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{7}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{7}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{7}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{7}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{7}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{8}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{8}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{8}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{9}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{9}.spm.stats.con.consess{1}.tcon.name = 'ss';
matlabbatch{9}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0];
matlabbatch{9}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{9}.spm.stats.con.consess{2}.tcon.name = 'cp';
matlabbatch{9}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 0 0];
matlabbatch{9}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{9}.spm.stats.con.consess{3}.tcon.name = 'cs';
matlabbatch{9}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 1 0];
matlabbatch{9}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{9}.spm.stats.con.consess{4}.tcon.name = 'us';
matlabbatch{9}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 1];
matlabbatch{9}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{9}.spm.stats.con.consess{5}.tcon.name = 'beta_sum';
matlabbatch{9}.spm.stats.con.consess{5}.tcon.weights = [0 1 1 1 1];
matlabbatch{9}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{9}.spm.stats.con.delete = 1;
matlabbatch{9}.spm.stats.con.consess{6}.fcon.name = 'condition_activation';
matlabbatch{9}.spm.stats.con.consess{6}.fcon.weights = [zeros(4,1) eye(4)];
matlabbatch{9}.spm.stats.con.consess{6}.fcon.sessrep = 'none';

%% run batch
regressor_batch = spm_jobman('run',matlabbatch);