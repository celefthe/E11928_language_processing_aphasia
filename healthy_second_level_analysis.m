%% fMRI second level data analysis for healthy subjects only

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

% set up dir structure for group-level analyses
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {[bids_dir filesep 'derivatives']};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'group';
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''group''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'healthy';


%% One-sample t-test for ss condition

matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'ss';

% factorial design for one-sample t-test on beta for ss
matlabbatch{4}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''ss''', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% get boosted beta for ss
boosted_ss_beta =  {};
for idx = 1:size(subj, 1)
     boosted_ss_beta{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_beta_0001.nii,1'];
end

matlabbatch{4}.spm.stats.factorial_design.des.t1.scans = boosted_ss_beta';

matlabbatch{4}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{4}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{4}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{4}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{4}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{4}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{4}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{4}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{5}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{6}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{6}.spm.stats.con.consess{1}.tcon.name = 'ss';
matlabbatch{6}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{6}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{6}.spm.stats.con.delete = 0;

%% One-sample t-test for cp effect

matlabbatch{7}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{7}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'cp';

% factorial design for one-sample t-test on beta for cs
matlabbatch{8}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''cp''', ...
    substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% get boosted beta for cp
boosted_cp_beta =  {};
for idx = 1:size(subj, 1)
     boosted_cp_beta{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_beta_0005.nii,1'];
end

matlabbatch{8}.spm.stats.factorial_design.des.t1.scans = boosted_cp_beta';

matlabbatch{8}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{8}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{8}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{8}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{8}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{8}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{8}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{8}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{9}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{9}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{9}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{10}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{10}.spm.stats.con.consess{1}.tcon.name = 'semantic';
matlabbatch{10}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{10}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{10}.spm.stats.con.delete = 0;


%% One-sample t-test for cs effect

matlabbatch{11}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{11}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'cs';

% factorial design for one-sample t-test on beta for cs
matlabbatch{12}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''cs''', ...
    substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% get boosted beta for cs
boosted_cs_beta =  {};
for idx = 1:size(subj, 1)
     boosted_cs_beta{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_beta_0009.nii,1'];
end

matlabbatch{12}.spm.stats.factorial_design.des.t1.scans = boosted_cs_beta';

matlabbatch{12}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{12}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{12}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{12}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{12}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{12}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{12}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{12}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{13}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{13}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{13}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{14}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{14}.spm.stats.con.consess{1}.tcon.name = 'cs';
matlabbatch{14}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{14}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{14}.spm.stats.con.delete = 0;


%% One-sample t-test for us effect

matlabbatch{15}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{15}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'unrelated';

% factorial design for one-sample t-test on beta for us
matlabbatch{16}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''us''', ...
    substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% get boosted beta for semantic
boosted_us_beta =  {};
for idx = 1:size(subj, 1)
     boosted_us_beta{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_beta_0013.nii,1'];
end

matlabbatch{16}.spm.stats.factorial_design.des.t1.scans = boosted_us_beta';

matlabbatch{16}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{16}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{16}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{16}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{16}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{16}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{16}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{16}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{17}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{17}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{17}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{18}.spm.stats.con.spmmat(1) = ...\
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{17}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{18}.spm.stats.con.consess{1}.tcon.name = 'us';
matlabbatch{18}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{18}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{18}.spm.stats.con.delete = 0;



%% run batch
stat_batch = spm_jobman('run',matlabbatch);