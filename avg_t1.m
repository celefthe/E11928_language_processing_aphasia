%% Generate average T1 image from subject data

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
     t1_images = [t1_images; subjects{idx}.t1img];
end

volume_headers = spm_vol(t1_images);
volume_data = spm_read_vols(volume_headers(1));
t1_average = zeros(size(volume_data));

for idx=1:size(t1_images,1),
   volume_data = spm_read_vols(volume_headers(idx));
   t1_average = t1_average + volume_data/size(t1_images,1);
end

volume_headers(1).fname = ...
    [bids_dir filesep 'derivatives' filesep 'group' filesep 'healthy' filesep 'anat' filesep 'healthy_t1avg.nii']
spm_write_vol(volume_headers(1),t1_average);
