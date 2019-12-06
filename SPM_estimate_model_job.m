%-----------------------------------------------------------------------
% Job saved on 18-Nov-2019 14:42:19 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% load the results_dir variable - to ensure the path we use in this script 
% matches the output folder specified in SPM_batch.m
load results_dir

matlabbatch{1}.spm.stats.fmri_est.spmmat = {[results_dir 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
