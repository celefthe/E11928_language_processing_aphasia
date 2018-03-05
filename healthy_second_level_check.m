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

%% QA Boost Effect

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {[bids_dir filesep 'derivatives']};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'group';


matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''group''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'healthy';


matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''healthy''', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'model_qa';


matlabbatch{4}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''model_qa''', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{4}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'boost_qa';


matlabbatch{5}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''boost_qa''', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

% Paired t-test for each hrf/hrf-boost pair of scans - look at task effect
for idx = 1:size(subj, 1)
     matlabbatch{5}.spm.stats.factorial_design.des.pt.pair(idx).scans = {
        [subjects{idx} filesep 'stats/con_0002.nii,1'],
        [subjects{idx} filesep 'boosted_stats/hrf_boost/sboost_con_0002.nii,1']
     };
end

matlabbatch{5}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{5}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{5}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{5}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{5}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{5}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{5}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{5}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{5}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{5}.spm.stats.factorial_design.globalm.glonorm = 1;


matlabbatch{6}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{6}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{6}.spm.stats.fmri_est.method.Classical = 1;


matlabbatch{7}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{7}.spm.stats.con.consess{1}.tcon.name = 'boost effect';

% contrast hrf vs hrf-boost 
matlabbatch{7}.spm.stats.con.consess{1}.tcon.weights = [-1 1 zeros(1, size(subj,1))];
matlabbatch{7}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{7}.spm.stats.con.delete = 1;

% run batch
boost_qa_batch = spm_jobman('run',matlabbatch);
clear matlabbatch;

%% Canonical HRF vs time derivative QA

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = boost_qa_batch{3}.dir;
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'hrf_vs_temporal';


matlabbatch{2}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''hrf_vs_temporal''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).name = 'conditions';
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).levels = 4;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).name = 'basis func';
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{2}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

% ss hrf
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(1).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(1).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0001.nii,1']};
end

% cp hrf
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0005.nii,1']};
end

% cs hrf§
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0009.nii,1']};
end

% us hrf
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).levels = [4
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0013.nii,1']};
end

% ss time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).levels = [1
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0002.nii,1']};
end

% cp time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).levels = [2
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0006.nii,1']};
end

% cs time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).levels = [3
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0010.nii,1']};
end

% us time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).levels = [4
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).scans(end+1,:) = ...
         {[subjects{idx} filesep 'stats/beta_0012.nii,1']};
end

matlabbatch{2}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{2}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{2}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.glonorm = 1;

% estimate model
matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

% generate contrasts - one for hrf only and an informed basis one incorporating
% the time derivative
matlabbatch{4}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{4}.spm.stats.con.consess{1}.fcon.name = 'hrf';
A = eye(8); A([2 4 6 8],[2 4 6 8])=0;
matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = A;
matlabbatch{4}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{2}.fcon.name = 'time deriv';
A = eye(8); A([1 3 5 7],[1 3 5 7])=0;
matlabbatch{4}.spm.stats.con.consess{2}.fcon.weights = A;
matlabbatch{4}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.delete = 0;


deriv_qa_batch = spm_jobman('run',matlabbatch);
