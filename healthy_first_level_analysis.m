%% fMRI data analysis for healthy subjects only

%% set up environment
% clear memoery, load SPM defaults and job manager
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

skip = 3; % how many initial volumes to skip

% the scripts assumes we run from inside the code directory
current = pwd;
if ~strcmp(current(end-3:end),'code')
    error('not running from code directory, CD to the right place then run the code')
end

cd ..
bids_dir = pwd;
subj = dir('sub*');
index = 1;

%% first level analysis
for subject = 1:size(subj,1)
    if strcmp(subj(subject).name(5:11),'healthy')
        % work in the derivative directory (create if not exist)
        mkdir([bids_dir filesep 'derivatives' filesep subj(subject).name])
        cd(subj(subject).name)
        
        % unzip anatomical and functional images
        cd('anat'); file = dir('*.nii.gz');
        T1name=gunzip(file.name, [bids_dir filesep 'derivatives' filesep subj(subject).name filesep 'anat']);
        cd ..
        cd('func'); 
        events = dir('*.tsv'); fmri_events = tdfread(events.name);
        file = dir('*.nii.gz'); 
        for f=1:size(file,1)
            if ~isempty(strfind(file(1).name,'audi'))
                bold_file = file(1).name;
            end
        end
        fMRIname=gunzip(bold_file, [bids_dir filesep 'derivatives' filesep subj(subject).name filesep 'func']);
        cd([bids_dir filesep 'derivatives' filesep subj(subject).name filesep 'func'])
        ThreeD_names = spm_file_split(spm_vol(cell2mat(fMRIname))); clear fMRI
        for f=1+skip:size(ThreeD_names,1)
            fMRI{f-skip} = ThreeD_names(f).fname;
        end
        
        % set the (0,0,0) coordinate
        fMRI = fMRI';
        P = fMRI;
        P{size(fMRI,1)+1}=cell2mat(T1name);
        spmup_auto_reorient(P,size(P,1));
        clear file bold_file events fMRIname P 
        
        %% preprocessing batch
        for s=(1+skip):size(ThreeD_names,1)
            matlabbatch{1}.spm.temporal.st.scans{1}((s-skip),:) = {ThreeD_names(s).fname};
        end        
        matlabbatch{1}.spm.temporal.st.nslices = 30;
        matlabbatch{1}.spm.temporal.st.tr = 2.5;
        matlabbatch{1}.spm.temporal.st.ta = 1.91666666666667;
        matlabbatch{1}.spm.temporal.st.so = ...
            [0 1.209 0.0806 1.2896 0.1612 1.3702 0.2418 1.4508 0.3224 1.5314 0.403 1.612 0.4836 1.6926 0.5642 1.7732 ...
            0.6448 1.8538 0.7254 1.9344 0.806 2.015 0.8866 2.0956 0.9672 2.1762 1.0478 2.2568 1.1284 2.3374];
        matlabbatch{1}.spm.temporal.st.refslice = matlabbatch{1}.spm.temporal.st.so(2);
        matlabbatch{1}.spm.temporal.st.prefix = 'sl-timing_';
        
        matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = ...
            cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val',...
            '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 1;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'rp_';
        
        matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = ...
            cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','rmean'));
        matlabbatch{3}.spm.spatial.coreg.estimate.source = T1name;
        matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = ...
            [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        
        matlabbatch{4}.spm.util.checkreg.data(1) = ...
            cfg_dep('Coregister: Estimate: Coregistered Images', ...
            substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','cfiles'));
        matlabbatch{5}.spm.util.print.fname = 'coregistered';
        matlabbatch{5}.spm.util.print.fig.figname = 'Graphics';
        matlabbatch{5}.spm.util.print.opts = 'pdf';
        
        matlabbatch{6}.spm.spatial.preproc.channel.vols(1) = ...
            cfg_dep('Coregister: Estimate: Coregistered Images', ...
            substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','cfiles'));
        matlabbatch{6}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{6}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{6}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{6}.spm.spatial.preproc.tissue(1).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,1']};
        matlabbatch{6}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{6}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(2).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,2']};
        matlabbatch{6}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{6}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(3).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,3']};
        matlabbatch{6}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{6}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(4).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,4']};
        matlabbatch{6}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{6}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(5).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,5']};
        matlabbatch{6}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{6}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(6).tpm = {[spmroot filesep 'tpm' filesep 'TPM.nii,6']};
        matlabbatch{6}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{6}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{6}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{6}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{6}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{6}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{6}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{6}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{6}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{6}.spm.spatial.preproc.warp.write = [0 1];
        
        matlabbatch{7}.spm.util.checkreg.data(1) = ...
            cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
        matlabbatch{7}.spm.util.checkreg.data(2) = cfg_dep('Segment: c2 Images', ...
            substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
        matlabbatch{8}.spm.util.print.fname = 'segmented';
        matlabbatch{8}.spm.util.print.fig.figname = 'Graphics';
        matlabbatch{8}.spm.util.print.opts = 'pdf';
        
        matlabbatch{9}.spm.spatial.normalise.write.subj.def(1) = ...
            cfg_dep('Segment: Forward Deformations', ...
            substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{9}.spm.spatial.normalise.write.subj.resample(1) = ...
            cfg_dep('Coregister: Estimate: Coregistered Images', ...
            substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','cfiles'));
        matlabbatch{9}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{9}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{9}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{9}.spm.spatial.normalise.write.woptions.prefix = 'normal_';
        
        matlabbatch{10}.spm.tools.oldnorm.estwrite.subj.source(1) = ...
            cfg_dep('Realign: Estimate & Reslice: Mean Image', ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','rmean'));
        matlabbatch{10}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
        matlabbatch{10}.spm.tools.oldnorm.estwrite.subj.resample(1) = ...
            cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','sess', '()',{1}, '.','cfiles'));
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.template = ...
            {[spmroot filesep 'toolbox' filesep 'OldNorm' filesep 'EPI.nii,1']};
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.interp = 4;
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{10}.spm.tools.oldnorm.estwrite.roptions.prefix = 'oldnorm_';
    
        matlabbatch{11}.spm.spatial.smooth.data(1) = ...
            cfg_dep('Old Normalise: Estimate & Write: Normalised Images (Subj 1)', ...
            substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','files'));
        matlabbatch{11}.spm.spatial.smooth.fwhm = [6 6 6];
        matlabbatch{11}.spm.spatial.smooth.dtype = 0;
        matlabbatch{11}.spm.spatial.smooth.im = 0;
        matlabbatch{11}.spm.spatial.smooth.prefix = 'smooth_';
        
        % low smooth to compensate for hrf boost
        matlabbatch{12}.spm.spatial.smooth.data(1) = ...
            cfg_dep('Old Normalise: Estimate & Write: Normalised Images (Subj 1)', ...
            substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('()',{1}, '.','files'));
        matlabbatch{12}.spm.spatial.smooth.fwhm = [2 2 2];
        matlabbatch{12}.spm.spatial.smooth.dtype = 0;
        matlabbatch{12}.spm.spatial.smooth.im = 0;
        matlabbatch{12}.spm.spatial.smooth.prefix = 'low_smooth_';
        
        try
            preprocess = spm_jobman('run',matlabbatch);
        catch exc
           warning('Could not preprocess %s - %s', subj(subject).name, exc.message);
           cd(bids_dir);
           clear matlabbatch;
           continue;
        end
        
        clear matlabbatch;
        
        
        %% stats batch
        flags = struct('motion_parameters','on','globals','on','volume_distance','off','movie','off');
        new_file = spmup_realign_qa(preprocess{11}.files,flags);
        
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
            {[bids_dir filesep 'derivatives' filesep subj(subject).name]};
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'stats';
        
        matlabbatch{2}.spm.stats.fmri_spec.dir(1) = ...
            cfg_dep('Make Directory: Make Directory ''stats''', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','dir'));
        matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 2.5;
        matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = 30;
        matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 15;
        for f=1:size(preprocess{11}.files,1)
            matlabbatch{2}.spm.stats.fmri_spec.sess.scans(f,:) = preprocess{11}.files(f);
        end
        
        % fmri_events.Condition = 1,2,3,4 = SS, CP,CS, US
        matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name = 'SS';
        matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name = 'CP';
        matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name = 'CS';
        matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).name = 'US';

        for cond = 1:4
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).onset = fmri_events.Onset(fmri_events.Condition==cond);
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).duration = ...
                mean(fmri_events.Duration(fmri_events.Condition==cond));
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).tmod = 0;
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).pmod.name = 'RT';
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).pmod.param = ...
                fmri_events.Duration(fmri_events.Condition==cond) - ...
                mean(fmri_events.Duration(fmri_events.Condition==cond));
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).pmod.poly = 1;
            matlabbatch{2}.spm.stats.fmri_spec.sess.cond(cond).orth = 1;
        end
        
         if isfield(fmri_events,'Error_onset')
            charn = size(fmri_events.Error_onset); % because of 'nil' sees entries as chars            
            errors = sum(fmri_events.Error_onset == ['nil' repmat(' ',[1 charn(2)-3])],2) ~= charn(2);
            if sum(errors)~=0
                matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).name = 'Errors';
                matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).onset = str2num(fmri_events.Error_onset(errors,:));
                matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).duration = str2num(fmri_events.Error_duration(errors,:));
                matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
            end
        end

        matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg(1) = {new_file};
        matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';

        matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = ...
            cfg_dep('fMRI model specification: SPM.mat File', ...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

        matlabbatch{4}.spm.stats.con.spmmat(1) = ...
            cfg_dep('Model estimation: SPM.mat File', ...
            substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','spmmat'));
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'task';
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
        matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'task-hrf';
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'task trial variations';
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0];
        matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'incongruency';
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = ...
            [-1 0 -1 0 ...
            0.333333333333333 0 0.333333333333333 0 ...
            0.333333333333333 0 0.333333333333333 0 ...
            0.333333333333333 0 0.333333333333333 0];
        matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{5}.tcon.name = 'incongruency hrf';
        matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = ...
            [-1 0 0 0 ...
            0.333333333333333 0 0 0 ...
            0.333333333333333 0 0 0 ...
            0.333333333333333 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{6}.tcon.name = 'incongruency trial variations';
        matlabbatch{4}.spm.stats.con.consess{6}.tcon.weights = ...
            [0 0 -1 0 ...
            0 0 0.333333333333333 0 ...
            0 0 0.333333333333333 0 ...
            0 0 0.333333333333333 0];
        matlabbatch{4}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{7}.tcon.name = 'phonology-hrf';
        matlabbatch{4}.spm.stats.con.consess{7}.tcon.weights = [-1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{8}.tcon.name = 'semantic-hrf';
        matlabbatch{4}.spm.stats.con.consess{8}.tcon.weights = [-1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.consess{9}.tcon.name = 'unrelated-hrf';
        matlabbatch{4}.spm.stats.con.consess{9}.tcon.weights = [-1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
        matlabbatch{4}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
        matlabbatch{4}.spm.stats.con.delete = 0;
        
        try
            stat_batch = spm_jobman('run',matlabbatch);
        catch exc
           warning('Could not calculate stats for %s - %s', subj(subject).name, exc.message);
           cd(bids_dir);
           clear matlabbatch;
           continue;
        end
        

        %% boosted stats batch
        % re-run the stats batch but on un-smoothed data and smooth the
        % boosted beta and con

        % boost stats
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = ...
            {[bids_dir filesep 'derivatives' filesep subj(subject).name]};
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'boosted_stats';

        matlabbatch{2}.spm.stats.fmri_spec.dir(1) = ...
            cfg_dep('Make Directory: Make Directory ''boosted_stats''', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
            substruct('.','dir'));

        for f=1:size(preprocess{12}.files,1)
            matlabbatch{2}.spm.stats.fmri_spec.sess.scans(f,:) = preprocess{12}.files(f);
        end

        try
            boost_batch = spm_jobman('run',matlabbatch);
        catch exc
           warning('Could not boost stats for %s - %s', subj(subject).name, exc.message);
           cd(bids_dir);
           clear matlabbatch;
           continue;
        end

        boosted_files = spmup_hrf_boost([bids_dir filesep 'derivatives' ...
            filesep subj(subject).name filesep 'boosted_stats' filesep 'SPM.mat']);
        spmup_smooth_boostedfiles(boosted_files{1},[8 8 8]); % smooth boosted beta files
        spmup_smooth_boostedfiles(boosted_files{2},[8 8 8]); % smooth boosted con files
        clear matlabbatch;
        cd(bids_dir);
    end
end
