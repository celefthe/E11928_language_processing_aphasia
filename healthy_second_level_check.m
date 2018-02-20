%% fMRI second level data analysis -- check if model works 

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

% get healthy subject derivative directory paths
subj = dir('sub-healthy*');
subjects = {};
for idx = 1:size(subj,1)
     subjects{idx} = [bids_dir filesep 'derivatives' filesep subj(idx).name];
end

%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {[bids_dir filesep 'derivatives']};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'group';
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''group''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'healthy';
matlabbatch{3}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% Paired t-test for each hrf/hrf-boost pair of scans - look at task effect
for idx = 1:size(subj, 1)
     matlabbatch{3}.spm.stats.factorial_design.des.pt.pair(idx).scans = {
        [subjects{idx} filesep 'stats/con_0002.nii,1'],
        [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0002.nii,1']
     };
end

matlabbatch{3}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{3}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{3}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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
matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'boost effect';

% contrast hrf vs hrf-boost 
matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = [-1 1 zeros(1, size(subj,1))];
matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{5}.spm.stats.con.delete = 0;

% run batch
stat_batch = spm_jobman('run',matlabbatch);
