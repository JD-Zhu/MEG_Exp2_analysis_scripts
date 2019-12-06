%-----------------------------------------------------------------------
% Job saved on 06-Dec-2019 18:40:17 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% load the results_dir variable - to ensure the path we use in this script 
% matches the output folder specified in SPM_batch.m
load results_dir

matlabbatch{1}.spm.stats.results.spmmat = {[results_dir 'SPM.mat']};

matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Intercept';
matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec(1).thresh = 0.1;
matlabbatch{1}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{1}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec(1).mask.none = 1;

matlabbatch{1}.spm.stats.results.conspec(2).titlestr = 'MainEffect_Context';
matlabbatch{1}.spm.stats.results.conspec(2).contrasts = 2;
matlabbatch{1}.spm.stats.results.conspec(2).threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec(2).thresh = 0.1;
matlabbatch{1}.spm.stats.results.conspec(2).extent = 0;
matlabbatch{1}.spm.stats.results.conspec(2).conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec(2).mask.none = 1;

matlabbatch{1}.spm.stats.results.conspec(3).titlestr = 'MainEffect_Ttype';
matlabbatch{1}.spm.stats.results.conspec(3).contrasts = 3;
matlabbatch{1}.spm.stats.results.conspec(3).threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec(3).thresh = 0.1;
matlabbatch{1}.spm.stats.results.conspec(3).extent = 0;
matlabbatch{1}.spm.stats.results.conspec(3).conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec(3).mask.none = 1;

matlabbatch{1}.spm.stats.results.conspec(4).titlestr = 'Interaction';
matlabbatch{1}.spm.stats.results.conspec(4).contrasts = 4;
matlabbatch{1}.spm.stats.results.conspec(4).threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec(4).thresh = 0.1;
matlabbatch{1}.spm.stats.results.conspec(4).extent = 0;
matlabbatch{1}.spm.stats.results.conspec(4).conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec(4).mask.none = 1;

matlabbatch{1}.spm.stats.results.units = 2;
matlabbatch{1}.spm.stats.results.export{1}.ps = true;
matlabbatch{1}.spm.stats.results.export{2}.png = true;
matlabbatch{1}.spm.stats.results.export{3}.xls = true;
