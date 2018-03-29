%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate Average Images for Healthy Participants   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
spm_jobman('initcfg');
spmroot = fileparts(which('spm'));
spm('Defaults','fmri')

% if matlab is ran from the cmd, prevent graphics window from opening
if ~usejava('desktop')
    spm_get_defaults('cmdline',true);
end

% set memory settings to improve performance 
spm_get_defaults('stats.maxmem', 2^31); 
spm_get_defaults('stats.resmem', true);

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

for idx = 1:size(subj,1)
     subjects{idx}.path = [bids_dir filesep 'derivatives' filesep subj(idx).name];
     subjects{idx}.t1img = [subjects{idx}.path filesep 'anat' filesep 'normal_' subj(idx).name '_ses-01_T1w.nii'];

     if idx == 1
         t1_images = [subjects{idx}.t1img];
     else
         t1_images = [t1_images; subjects{idx}.t1img];
     end
end

%% Generate average T1 image from subject data
volume_headers = spm_vol(t1_images);
volume_data = spm_read_vols(volume_headers);
t1_average = mean(volume_data,4);
W = volume_headers(1);
W.fname = [bids_dir filesep 'derivatives' filesep 'group' filesep 'healthy_t1avg.nii'];
W.descrip = 'Average of 16 healthy normalized T1 controls';
W.private.descrip = W.descrip;
spm_write_vol(W,t1_average);


%% Segment average T1 image
matlabbatch{1}.spm.spatial.preproc.channel.vols(1) = {W.fname};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

segmentation_batch = spm_jobman('run',matlabbatch);

%% Create surface skullstripped image from normalised c1, c2 & c3 images

[c1 c2 c3 c4 c5] = segmentation_batch{1}.tiss.c

c1_header = spm_vol(c1);
c1_data = spm_read_vols(c1_header{1});
c2_header = spm_vol(c2);
c2_data = spm_read_vols(c2_header{1});
c3_header = spm_vol(c3);
c3_data = spm_read_vols(c3_header{1});

% Surface extraction using c1 & c2
spm_surf([c1{1}; c2{1}], 3);

% Surface extraction using c1/c2 & skullstripped T1
%TODO - revise below
t1_surf = (c1_data + c2_data) .* (t1_average - ((c1_data + c2_data + c3_data))>0);

W.fname = [bids_dir filesep 'derivatives' filesep 'group' filesep 'surf_healthy_t1avg.nii'];
W.descrip = 'Average of 16 healthy normalized T1 controls';
W.private.descrip = W.descrip;
spm_write_vol(W,t1_surf);

spm_surf(W.fname, 3);
