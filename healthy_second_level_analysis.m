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
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Make Directory: Make Directory ''group''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'healthy';


%% One-sample t-test for  effect

matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Make Directory: Make Directory ''healthy''', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'incongruency';

% factorial design for one-sample t-test on delta constrast of incogruency
% effect
matlabbatch{4}.spm.stats.factorial_design.dir(1) = cfg_dep('Make Directory: Make Directory ''incongruency''', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));

% get boosted contrast for incongruency
boosted_inc_con =  {};
for idx = 1:size(subj, 1)
     boosted_inc_con{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0004.nii,1'];
end

matlabbatch{4}.spm.stats.factorial_design.des.t1.scans = boosted_inc_con';

matlabbatch{4}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{4}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{4}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{4}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{4}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{4}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{4}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{4}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{6}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.con.consess{1}.tcon.name = 'incongruency';
matlabbatch{6}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{6}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{6}.spm.stats.con.delete = 0;

%% One-sample t-test for semantic effect

matlabbatch{7}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Make Directory: Make Directory ''healthy''', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
matlabbatch{7}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'semantic';

% factorial design for one-sample t-test on delta constrast of semantic
% effect
matlabbatch{8}.spm.stats.factorial_design.dir(1) = cfg_dep('Make Directory: Make Directory ''semantic''', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));

% get boosted contrast for semantic
boosted_semantic =  {};
for idx = 1:size(subj, 1)
     boosted_semantic{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0008.nii,1'];
end

matlabbatch{8}.spm.stats.factorial_design.des.t1.scans = boosted_semantic';

matlabbatch{8}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{8}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{8}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{8}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{8}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{8}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{8}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{8}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{9}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{9}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{9}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{10}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{10}.spm.stats.con.consess{1}.tcon.name = 'semantic';
matlabbatch{10}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{10}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{10}.spm.stats.con.delete = 0;


%% One-sample t-test for phonology effect

matlabbatch{11}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Make Directory: Make Directory ''healthy''', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
matlabbatch{11}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'phonology';

% factorial design for one-sample t-test on delta constrast of phonology
% effect
matlabbatch{12}.spm.stats.factorial_design.dir(1) = cfg_dep('Make Directory: Make Directory ''phonology''', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));

% get boosted contrast for phono
boosted_phono =  {};
for idx = 1:size(subj, 1)
     boosted_phono{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0007.nii,1'];
end

matlabbatch{12}.spm.stats.factorial_design.des.t1.scans = boosted_phono';

matlabbatch{12}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{12}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{12}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{12}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{12}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{12}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{12}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{12}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{13}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{13}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{13}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{14}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{14}.spm.stats.con.consess{1}.tcon.name = 'phonology';
matlabbatch{14}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{14}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{14}.spm.stats.con.delete = 0;


%% One-sample t-test for unrelated effect

matlabbatch{15}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Make Directory: Make Directory ''healthy''', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
matlabbatch{15}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'unrelated';

% factorial design for one-sample t-test on delta constrast of unrelated
% effect
matlabbatch{16}.spm.stats.factorial_design.dir(1) = cfg_dep('Make Directory: Make Directory ''unrelated''', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));

% get boosted contrast for semantic
boosted_unrelated =  {};
for idx = 1:size(subj, 1)
     boosted_unrelated{idx} = [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0009.nii,1'];
end

matlabbatch{16}.spm.stats.factorial_design.des.t1.scans = boosted_semantic';

matlabbatch{16}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{16}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{16}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{16}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{16}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{16}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{16}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{16}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{17}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{17}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{17}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{18}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{17}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{18}.spm.stats.con.consess{1}.tcon.name = 'unrelated';
matlabbatch{18}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{18}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{18}.spm.stats.con.delete = 0;



%% run batch
stat_batch = spm_jobman('run',matlabbatch);