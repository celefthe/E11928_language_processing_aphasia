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
     if idx == 1
         t1_images = [subjects{idx}.t1img];
     else
         t1_images = [t1_images; subjects{idx}.t1img];
     end
end

volume_headers = spm_vol(t1_images);
volume_data = spm_read_vols(volume_headers);
t1_average = mean(volume_data,4);
W = volume_headers(1);
W.fname = [bids_dir filesep 'derivatives' filesep 'group' filesep 'healthy_t1avg.nii'];
W.descrip = 'Average of 16 healthy normalized T1 controls';
W.private.descrip = W.descrip;
spm_write_vol(W,t1_average);


