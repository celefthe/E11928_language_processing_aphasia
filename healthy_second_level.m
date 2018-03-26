%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fMRI Second Level Analysis for Healthy Participants %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up SPM environment
clear variables
spm_jobman('initcfg');
spmroot = fileparts(which('spm'));
spm('Defaults','fmri');

% if matlab is ran from the cmd, prevent graphics window from opening
%if ~usejava('desktop')
%    spm_get_defaults('cmdline',true);
%end

% set memory settings to improve performance 
spm_get_defaults('stats.maxmem', 2^31); 
spm_get_defaults('stats.resmem', true);

% the scripts assumes we run from inside the code directory
current = pwd;
if ~strcmp(current(end-3:end),'code')
    error('not running from code directory, CD to the right place then run the code');
end

% get to the bids root dir
cd ..
bids_dir = pwd;

% get participants.tsv file containing behavioural parameters
participants = readtable([bids_dir filesep 'participants.tsv'], 'FileType', 'text', 'Delimiter', '\t');

% get healthy subject derivative directory paths and assign behavioural modelling params to each subject
subj = dir('sub-healthy*');
subjects = {};
for idx = 1:size(subj,1)
     subjects{idx}.path = [bids_dir filesep 'derivatives' filesep subj(idx).name];

     subjects{idx}.beta.ss = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0001.nii,1'];
     subjects{idx}.beta.cp = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0005.nii,1'];
     subjects{idx}.beta.cs = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0009.nii,1'];
     subjects{idx}.beta.us = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_beta_0013.nii,1'];

     subjects{idx}.con.phonology = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0007.nii,1'];
     subjects{idx}.con.semantic = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0008.nii,1'];
     subjects{idx}.con.unrelated = [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0009.nii,1'];
     
     subjects{idx}.dprime.phonology = participants.dprimephono_sess_1(idx);
     subjects{idx}.dprime.semantic = participants.dprimesem_sess_1(idx);
     subjects{idx}.dprime.unrelated = participants.dprimeunrelated_sess_1(idx);
     
     subjects{idx}.beta.ss_time = [subjects{idx}.path filesep 'stats/beta_0002.nii,1'];
     subjects{idx}.beta.cp_time = [subjects{idx}.path filesep 'stats/beta_0006.nii,1'];
     subjects{idx}.beta.cs_time = [subjects{idx}.path filesep 'stats/beta_0010.nii,1'];
     subjects{idx}.beta.us_time = [subjects{idx}.path filesep 'stats/beta_0014.nii,1'];
     
     subjects{idx}.subj_idx = participants.id(idx);
     subjects{idx}.drift.ss = participants.drifratess_sess_1(idx);
     subjects{idx}.drift.cp = participants.drifratecp_sess_1(idx);
     subjects{idx}.drift.cs = participants.drifratecs_sess_1(idx);
     subjects{idx}.drift.us = participants.drifrateus_sess_1(idx);
     subjects{idx}.threshold = participants.ddmthreshold_sess_1(idx);
     subjects{idx}.bias = participants.ddmbias_sess_1(idx);
     subjects{idx}.nondec = participants.ddmnondec_sess_1(idx);

end

% set up dir structure for group-level analyses
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {[bids_dir filesep 'derivatives']};
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'group';
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''group''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'healthy';

setup_batch = spm_jobman('run', matlabbatch);
clear matlabbatch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Model QA                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% QA Task Effect
% Check task effect for standard and boosted hrf

matlabbatch{4}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''model_qa''', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{4}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'task_effect_qa';

matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''task_effect_qa''', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'hrf';


% One-sample t-test hrf for task effect
matlabbatch{6}.spm.stats.factorial_design.dir = ...
    cfg_dep('Make Directory: Make Directory ''hrf''', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{6}.spm.stats.factorial_design.des.t1.scans = {};
for idx = 1:size(subj,1)
    matlabbatch{6}.spm.stats.factorial_design.des.t1.scans(end+1,:) = {
        [subjects{idx}.path filesep 'stats/con_0001.nii,1']
    };
end
matlabbatch{6}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{6}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{6}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{6}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{6}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{6}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{6}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{6}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{7}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{7}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{7}.spm.stats.fmri_est.method.Classical = 1;


% One-sample t-test boosted hrf for task effect
matlabbatch{8}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''task_effect_qa''', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{8}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'boosted_hrf';

matlabbatch{9}.spm.stats.factorial_design.dir = ...
    cfg_dep('Make Directory: Make Directory ''boosted_hrf''', ...
    substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{9}.spm.stats.factorial_design.des.t1.scans = {};
for idx = 1:size(subj,1)
    matlabbatch{9}.spm.stats.factorial_design.des.t1.scans(end+1,:) = {
        [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0001.nii,1']
    };
end
matlabbatch{9}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{9}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{9}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{9}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{9}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{9}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{9}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{9}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{10}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{10}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{10}.spm.stats.fmri_est.method.Classical = 1;

% run batch
task_hrf_batch = spm_jobman('run',matlabbatch);
clear matlabbatch;


%% QA Boost Effect

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = task_hrf_batch{3}.dir;
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'boost_qa';

% Paired t-test for each hrf/hrf-boost pair of scans - look at task effect
matlabbatch{2}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''boost_qa''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{2}.spm.stats.factorial_design.des.pt.pair(idx).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.pt.pair(idx).scans = {
        [subjects{idx}.path filesep 'stats/con_0001.nii,1'],
        [subjects{idx}.path filesep 'boosted_stats/hrf_boost/sboost_con_0001.nii,1']
     };
end
matlabbatch{2}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{2}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{2}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{2}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.glonorm = 1;


matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;


% contrast hrf vs hrf-boost 
matlabbatch{4}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'boost effect';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [-1 1 zeros(1, size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.delete = 1;

matlabbatch{5}.spm.stats.results.spmmat = ...
    cfg_dep('Contrast Manager: SPM.mat File', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;  % first contrast is hrf vs hrf_boost
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.export{1}.binary.basename = 'binary';

% run batch
boost_qa_batch = spm_jobman('run',matlabbatch);
clear matlabbatch;


%% Canonical HRF vs time derivative ANOVA QA

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = task_hrf_batch{3}.dir;
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'hrf_deriv';

matlabbatch{2}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''hrf_deriv''', ...
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
         {[subjects{idx}.path filesep 'stats/beta_0001.nii,1']};
end

% cp hrf
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(2).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0005.nii,1']};
end

% cs hrf§
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(3).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0009.nii,1']};
end

% us hrf
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).levels = [4
                                                                    1];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(4).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0013.nii,1']};
end

% ss time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).levels = [1
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(5).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0002.nii,1']};
end

% cp time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).levels = [2
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(6).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0006.nii,1']};
end

% cs time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).levels = [3
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(7).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0010.nii,1']};
end

% us time deriv
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).levels = [4
                                                                    2];
matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).scans = {};
for idx = 1:size(subj, 1)
     matlabbatch{2}.spm.stats.factorial_design.des.fd.icell(8).scans(end+1,:) = ...
         {[subjects{idx}.path filesep 'stats/beta_0012.nii,1']};
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
matlabbatch{4}.spm.stats.con.delete = 1;

matlabbatch{5}.spm.stats.results.spmmat = ...
    cfg_dep('Contrast Manager: SPM.mat File', ...
    substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.contrast.contrasts = 2;
matlabbatch{5}.spm.stats.results.conspec.mask.contrast.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.mask.contrast.mtype = 0;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.ps = true;
matlabbatch{5}.spm.stats.results.export{2}.tspm.basename = 'overlay';
matlabbatch{5}.spm.stats.results.export{3}.binary.basename = 'binary';

deriv_qa_batch = spm_jobman('run',matlabbatch);
clear matlabbatch; 


%% Compute overlap between deriv_qa and boost_qa

deriv_qa_vol_headers = spm_vol(deriv_qa_batch{5}.filtered{1});
deriv_qa_vol_data = spm_read_vols(deriv_qa_vol_headers);

boost_qa_vol_headers = spm_vol(boost_qa_batch{5}.filtered{1});
boost_qa_vol_data = spm_read_vols(boost_qa_vol_headers);

overlap = sum(boost_qa_vol_data(:)>0 & deriv_qa_vol_data(:)>0) / ...
        (length((boost_qa_vol_data(:)>0)) + length((deriv_qa_vol_data(:)>0)));
fprintf('\nderiv_qa/boost_qa overlap: %3.4f %%\n\n', overlap * 100);

fileid = fopen([task_hrf_batch{3}.dir{1} filesep 'deriv-boost_overlap.txt'],'w');
fprintf(fileid,'%3.4f %%\n', overlap * 100);
fclose(fileid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Model Incongruency Effects %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Incongruence effect - boosted hrf

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = task_hrf_batch{2}.dir;
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'incongruence_effect';

matlabbatch{2}.spm.stats.factorial_design.dir = ...
    cfg_dep('Make Directory: Make Directory ''incongruence_effect''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

for idx = 1:size(subj,1)
    matlabbatch{2}.spm.stats.factorial_design.des.anovaw.fsubject(idx).scans = ...
        {[subjects{idx}.beta.ss]; [subjects{idx}.beta.cp]; [subjects{idx}.beta.cs]; [subjects{idx}.beta.us]};
    matlabbatch{2}.spm.stats.factorial_design.des.anovaw.fsubject(idx).conds = [1 2 3 4]; % SS, CP, CS, US
end
matlabbatch{2}.spm.stats.factorial_design.des.anovaw.dept = 1;
matlabbatch{2}.spm.stats.factorial_design.des.anovaw.variance = 1;
matlabbatch{2}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{2}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{2}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{2}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{4}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));

% contrast for overall incongruence
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'ss_all-';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 -1/3 -1/3 -1/3 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'ss_all+';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 1/3 1/3 1/3 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

% ss vs specific incongruence pair contrasts
matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'ss_cp-';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 -1 0 0 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'ss_cs-';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [1 0 -1 0 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{5}.tcon.name = 'ss_us-';
matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [1 0 0 -1 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{6}.tcon.name = 'ss_cp+';
matlabbatch{4}.spm.stats.con.consess{6}.tcon.weights = [-1 1 0 0 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{7}.tcon.name = 'ss_cs+';
matlabbatch{4}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1 0 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{8}.tcon.name = 'ss_us+';
matlabbatch{4}.spm.stats.con.consess{8}.tcon.weights = [-1 0 0 1 zeros(1,size(subj,1))];
matlabbatch{4}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

% subject effect contrast
%matlabbatch{4}.spm.stats.con.consess{5}.tcon.name = 'subj_effect';
%matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [zeros(1,4) ones(1,size(subj,1))*[1/size(subj,1)]];
%matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

% subject effect regressions with conditions
matlabbatch{4}.spm.stats.con.consess{9}.tcon.name = 'ss_regr+';
matlabbatch{4}.spm.stats.con.consess{9}.tcon.weights = [1 0 0 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{10}.tcon.name = 'cp_regr+';
matlabbatch{4}.spm.stats.con.consess{10}.tcon.weights =  [0 1 0 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{11}.tcon.name = 'cs_regr+';
matlabbatch{4}.spm.stats.con.consess{11}.tcon.weights = [0 0 1 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{12}.tcon.name = 'us_regr+';
matlabbatch{4}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 1 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

% negative subject effect regressions with conditions
matlabbatch{4}.spm.stats.con.consess{13}.tcon.name = 'ss_regr-';
matlabbatch{4}.spm.stats.con.consess{13}.tcon.weights = [-1 0 0 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{14}.tcon.name = 'cp_regr-';
matlabbatch{4}.spm.stats.con.consess{14}.tcon.weights =  [0 -1 0 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{15}.tcon.name = 'cs_regr-';
matlabbatch{4}.spm.stats.con.consess{15}.tcon.weights = [0 0 -1 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

matlabbatch{4}.spm.stats.con.consess{16}.tcon.name = 'us_regr-';
matlabbatch{4}.spm.stats.con.consess{16}.tcon.weights = [0 0 0 -1 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{4}.spm.stats.con.consess{16}.tcon.sessrep = 'none';


matlabbatch{4}.spm.stats.con.delete = 0;


%% Incongruence effect - time deriv
matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
    cfg_dep('Make Directory: Make Directory ''incongruence_effect''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'time_deriv';

matlabbatch{6}.spm.stats.factorial_design.dir = ...
    cfg_dep('Make Directory: Make Directory ''time_deriv''', ...
    substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));

for idx = 1:size(subj,1)
    matlabbatch{6}.spm.stats.factorial_design.des.anovaw.fsubject(idx).scans = ...
        {[subjects{idx}.beta.ss_time]; [subjects{idx}.beta.cp_time]; 
        [subjects{idx}.beta.cs_time]; [subjects{idx}.beta.us_time]};
    matlabbatch{6}.spm.stats.factorial_design.des.anovaw.fsubject(idx).conds = [1 2 3 4]; % SS, CP, CS, US
end
matlabbatch{6}.spm.stats.factorial_design.des.anovaw.dept = 1;
matlabbatch{6}.spm.stats.factorial_design.des.anovaw.variance = 1;
matlabbatch{6}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{6}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{6}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{6}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{6}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{6}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{6}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{6}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{6}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{6}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{7}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{7}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{7}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{8}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));

% contrast for overall incongruence
matlabbatch{8}.spm.stats.con.consess{1}.tcon.name = 'ss_all-';
matlabbatch{8}.spm.stats.con.consess{1}.tcon.weights = [1 -1/3 -1/3 -1/3 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{2}.tcon.name = 'ss_all+';
matlabbatch{8}.spm.stats.con.consess{2}.tcon.weights = [-1 1/3 1/3 1/3 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

% ss vs specific incongruence pair contrasts
matlabbatch{8}.spm.stats.con.consess{3}.tcon.name = 'ss_cp-';
matlabbatch{8}.spm.stats.con.consess{3}.tcon.weights = [1 -1 0 0 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{8}.spm.stats.con.consess{4}.tcon.name = 'ss_cs-';
matlabbatch{8}.spm.stats.con.consess{4}.tcon.weights = [1 0 -1 0 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{8}.spm.stats.con.consess{5}.tcon.name = 'ss_us-';
matlabbatch{8}.spm.stats.con.consess{5}.tcon.weights = [1 0 0 -1 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{6}.tcon.name = 'ss_cp+';
matlabbatch{8}.spm.stats.con.consess{6}.tcon.weights = [-1 1 0 0 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{8}.spm.stats.con.consess{7}.tcon.name = 'ss_cs+';
matlabbatch{8}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1 0 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{8}.spm.stats.con.consess{8}.tcon.name = 'ss_us+';
matlabbatch{8}.spm.stats.con.consess{8}.tcon.weights = [-1 0 0 1 zeros(1,size(subj,1))];
matlabbatch{8}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

% subject effect regressions with conditions
matlabbatch{8}.spm.stats.con.consess{9}.tcon.name = 'ss_regr+';
matlabbatch{8}.spm.stats.con.consess{9}.tcon.weights = [1 0 0 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{10}.tcon.name = 'cp_regr+';
matlabbatch{8}.spm.stats.con.consess{10}.tcon.weights =  [0 1 0 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{11}.tcon.name = 'cs_regr+';
matlabbatch{8}.spm.stats.con.consess{11}.tcon.weights = [0 0 1 0 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{12}.tcon.name = 'us_regr+';
matlabbatch{8}.spm.stats.con.consess{12}.tcon.weights = [0 0 0 1 ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

% negative subject effect regressions with conditions
matlabbatch{8}.spm.stats.con.consess{13}.tcon.name = 'ss_regr-';
matlabbatch{8}.spm.stats.con.consess{13}.tcon.weights = [-1 0 0 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{14}.tcon.name = 'cp_regr-';
matlabbatch{8}.spm.stats.con.consess{14}.tcon.weights =  [0 -1 0 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{15}.tcon.name = 'cs_regr-';
matlabbatch{8}.spm.stats.con.consess{15}.tcon.weights = [0 0 -1 0 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

matlabbatch{8}.spm.stats.con.consess{16}.tcon.name = 'us_regr-';
matlabbatch{8}.spm.stats.con.consess{16}.tcon.weights = [0 0 0 -1 -ones(1,size(subj,1))*[1/size(subj,1)]];
matlabbatch{8}.spm.stats.con.consess{16}.tcon.sessrep = 'none';


matlabbatch{8}.spm.stats.con.delete = 0;


incongruence_batch = spm_jobman('run', matlabbatch);
clear matlabbatch;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Regressions                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

matlabbatch{3}.spm.stats.factorial_design.des.mreg.inc = 1;

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


%% regression with drift rate + nondec 

matlabbatch{10}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = ...
    cfg_dep('Make Directory: Make Directory ''regressions''', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{10}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'ddm_with_nondec';


matlabbatch{11}.spm.stats.factorial_design.dir(1) = ...
    cfg_dep('Make Directory: Make Directory ''ddm_with_nondec''', ...
    substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','dir'));
matlabbatch{11}.spm.stats.factorial_design.des.mreg.scans = {};
for idx=1:size(subjects,2)
    matlabbatch{11}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.ss};
    matlabbatch{11}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.cp};
    matlabbatch{11}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.cs};
    matlabbatch{11}.spm.stats.factorial_design.des.mreg.scans(end+1,:) = {subjects{idx}.beta.us};
end

matlabbatch{11}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});

drift_ss = [];
for idx = 1:size(subjects,2)
    subj_v = [subjects{idx}.drift.ss; 0; 0; 0];
    drift_ss = [drift_ss ; subj_v];
end
matlabbatch{11}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(1).c = drift_ss;
matlabbatch{11}.spm.stats.factorial_design.cov(1).cname = 'v_{ss}';
matlabbatch{11}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(1).iCC = 1;

drift_cp = [];
for idx = 1:size(subjects,2)
    subj_v = [0; subjects{idx}.drift.cp; 0; 0];
    drift_cp = [drift_cp ; subj_v];
end
matlabbatch{11}.spm.stats.factorial_design.cov(2).c = drift_cp;
matlabbatch{11}.spm.stats.factorial_design.cov(2).cname = 'v_{cp}';
matlabbatch{11}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(2).iCC = 1;

drift_cs = [];
for idx = 1:size(subjects,2)
    subj_v = [0; 0; subjects{idx}.drift.cs; 0];
    drift_cs = [drift_cs ; subj_v];
end
matlabbatch{11}.spm.stats.factorial_design.cov(3).c = drift_cs;
matlabbatch{11}.spm.stats.factorial_design.cov(3).cname = 'v_{cs}';
matlabbatch{11}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(3).iCC = 1;

drift_us = [];
for idx = 1:size(subjects,2)
    subj_v = [0; 0; 0; subjects{idx}.drift.us];
    drift_us = [drift_us ; subj_v];
end
matlabbatch{11}.spm.stats.factorial_design.cov(4).c = drift_us;
matlabbatch{11}.spm.stats.factorial_design.cov(4).cname = 'v_{us}';
matlabbatch{11}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(4).iCC = 1;

nondec = [];
for idx = 1:size(subjects,2)
    subj_t = [
        subjects{idx}.nondec;
        subjects{idx}.nondec;
        subjects{idx}.nondec;
        subjects{idx}.nondec];
    nondec = [nondec ; subj_t];
end
matlabbatch{11}.spm.stats.factorial_design.cov(5).c = nondec;
matlabbatch{11}.spm.stats.factorial_design.cov(5).cname = 'nondec';
matlabbatch{11}.spm.stats.factorial_design.cov(5).iCFI = 1;
matlabbatch{11}.spm.stats.factorial_design.cov(5).iCC = 1;


matlabbatch{11}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{11}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{11}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{11}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{11}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{11}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{11}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{12}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('Factorial design specification: SPM.mat File', ...
    substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{12}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{12}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{13}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File', ...
    substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{13}.spm.stats.con.consess{1}.tcon.name = 'v_ss';
matlabbatch{13}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0 0];
matlabbatch{13}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{2}.tcon.name = 'v_cp';
matlabbatch{13}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 0 0 0];
matlabbatch{13}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{3}.tcon.name = 'v_cs';
matlabbatch{13}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 1 0 0];
matlabbatch{13}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{4}.tcon.name = 'v_us';
matlabbatch{13}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 1 0];
matlabbatch{13}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{5}.tcon.name = 'v_beta_sum';
matlabbatch{13}.spm.stats.con.consess{5}.tcon.weights = [0 1 1 1 1 0];
matlabbatch{13}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{6}.fcon.name = 'v_condition_activation';
matlabbatch{13}.spm.stats.con.consess{6}.fcon.weights = [zeros(4,1) eye(4)];
matlabbatch{13}.spm.stats.con.consess{6}.fcon.sessrep = 'none';
matlabbatch{13}.spm.stats.con.consess{7}.tcon.name = 'nondec';
matlabbatch{13}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 1];
matlabbatch{13}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch{13}.spm.stats.con.delete = 1;

% run batch
regressor_batch = spm_jobman('run',matlabbatch);

% go back to bids code dir
cd(current)
