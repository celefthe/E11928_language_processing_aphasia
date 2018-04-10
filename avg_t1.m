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

     subjects{idx}.c1 = [subjects{idx}.path filesep 'anat' filesep 'c1' subj(idx).name '_ses-01_T1w.nii'];
     subjects{idx}.c2 = [subjects{idx}.path filesep 'anat' filesep 'c2' subj(idx).name '_ses-01_T1w.nii'];
     subjects{idx}.c3 = [subjects{idx}.path filesep 'anat' filesep 'c3' subj(idx).name '_ses-01_T1w.nii'];

     subjects{idx}.def = [subjects{idx}.path filesep 'anat' filesep 'y_' subj(idx).name '_ses-01_T1w.nii']; 
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



%% Normalise subject c1 & c2 images
for idx = 1:size(subj, 1)
    matlabbatch{1}.spm.spatial.normalise.write.subj(idx).def = {subjects{idx}.def};
    matlabbatch{1}.spm.spatial.normalise.write.subj(idx).resample = {
        subjects{idx}.c1;
        subjects{idx}.c2;
        subjects{idx}.c3
    };
end
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

seg_norm_batch = spm_jobman('run', matlabbatch);

%% Create surface skullstripped image from normalised c1, c2 & c3 images

seg_norm_files = {seg_norm_batch{1}.files}
for idx = 1:size(subj, 1)
    if idx == 1
        wc1_images = [seg_norm_files{idx}{1}];
        wc2_images = [seg_norm_files{idx}{2}];
        wc3_images = [seg_norm_files{idx}{3}];
    else 
        wc1_images = [wc1_images; seg_norm_files{idx}{1}];
        wc2_images = [wc2_images; seg_norm_files{idx}{2}];
        wc3_images = [wc3_images; seg_norm_files{idx}{3}];
    end
end

wc1_headers = spm_vol(wc1_images);
wc2_headers = spm_vol(wc2_images);
wc3_headers = spm_vol(wc3_images);

wc1_data = spm_read_vols(wc1_headers);
wc2_data = spm_read_vols(wc2_headers);
wc3_data = spm_read_vols(wc3_headers);

wc1_average = mean(wc1_data,4);
wc2_average = mean(wc2_data,4);
wc3_average = mean(wc3_data,4);

%spm_surf([wc1_images; wc2_images], 3);

t1_surf = (wc1_average + wc2_average) .* (t1_average - ((wc1_average + wc2_average + wc3_average))>0);
%c1_header = spm_vol(c1);
%c1_data = spm_read_vols(c1_header{1});
%c2_header = spm_vol(c2);
%c2_data = spm_read_vols(c2_header{1});
%c3_header = spm_vol(c3);
%c3_data = spm_read_vols(c3_header{1});

% Surface extraction using c1 & c2
%spm_surf([c1{1}; c2{1}], 3);

% Surface extraction using c1/c2 & skullstripped T1
%TODO - revise below
%t1_surf = (c1_data + c2_data) .* (t1_average - ((c1_data + c2_data + c3_data))>0);

W.fname = [bids_dir filesep 'derivatives' filesep 'group' filesep 'surfr_healthy_t1avg.nii'];
W.descrip = 'Average of 16 healthy normalized T1 controls';
W.private.descrip = W.descrip;
spm_write_vol(W,t1_surf);

spm_surf(W.fname, 3);
