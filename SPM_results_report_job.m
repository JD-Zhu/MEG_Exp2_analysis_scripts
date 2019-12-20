%-----------------------------------------------------------------------
% Job saved on 06-Dec-2019 18:40:17 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% load the results_dir variable - to ensure the path we use in this script 
% matches the output folder specified in SPM_batch.m
load results_dir
% load the which_contrasts settings, so that we know which contrasts to
% report
load which_contrasts


matlabbatch{1}.spm.stats.results.spmmat = {[results_dir 'SPM.mat']};


if strcmp(which_contrasts, 'bi_only')
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Switch_cost';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 1;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(1).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(1).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(1).mask.none = 1;

    matlabbatch{1}.spm.stats.results.conspec(2).titlestr = 'Mixing_cost';
    matlabbatch{1}.spm.stats.results.conspec(2).contrasts = 2;
    matlabbatch{1}.spm.stats.results.conspec(2).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(2).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(2).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(2).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(2).mask.none = 1;
    
elseif (strcmp(which_contrasts, '3x3') || strcmp(which_contrasts, '3x2')) 
    
    % for '3x3', there are 3 contrasts to report (below)
    % for '3x2', we can use the same spec to report the sw$ first
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'MainEffect_Context';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 1;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(1).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(1).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(1).mask.none = 1;

    matlabbatch{1}.spm.stats.results.conspec(2).titlestr = 'MainEffect_Ttype';
    matlabbatch{1}.spm.stats.results.conspec(2).contrasts = 2;
    matlabbatch{1}.spm.stats.results.conspec(2).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(2).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(2).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(2).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(2).mask.none = 1;

    matlabbatch{1}.spm.stats.results.conspec(3).titlestr = 'Interaction';
    matlabbatch{1}.spm.stats.results.conspec(3).contrasts = 3;
    matlabbatch{1}.spm.stats.results.conspec(3).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(3).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(3).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(3).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(3).mask.none = 1;

    % Overall effect - we no longer tested this coz we are not interested
    %{
    matlabbatch{1}.spm.stats.results.conspec(4).titlestr = 'Intercept';
    matlabbatch{1}.spm.stats.results.conspec(4).contrasts = 4;
    matlabbatch{1}.spm.stats.results.conspec(4).threshdesc = 'FWE';
    matlabbatch{1}.spm.stats.results.conspec(4).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
    matlabbatch{1}.spm.stats.results.conspec(4).extent = 0;
    matlabbatch{1}.spm.stats.results.conspec(4).conjunction = 1;
    matlabbatch{1}.spm.stats.results.conspec(4).mask.none = 1;
    %}


    % for '3x2', we now add the reporting of mix$
    if strcmp(which_contrasts, '3x2') 
        matlabbatch{1}.spm.stats.results.conspec(4).titlestr = 'MIX - MainEffect_Context';
        matlabbatch{1}.spm.stats.results.conspec(4).contrasts = 4;
        matlabbatch{1}.spm.stats.results.conspec(4).threshdesc = 'FWE';
        matlabbatch{1}.spm.stats.results.conspec(4).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
        matlabbatch{1}.spm.stats.results.conspec(4).extent = 0;
        matlabbatch{1}.spm.stats.results.conspec(4).conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec(4).mask.none = 1;

        matlabbatch{1}.spm.stats.results.conspec(5).titlestr = 'MIX - MainEffect_Ttype';
        matlabbatch{1}.spm.stats.results.conspec(5).contrasts = 5;
        matlabbatch{1}.spm.stats.results.conspec(5).threshdesc = 'FWE';
        matlabbatch{1}.spm.stats.results.conspec(5).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
        matlabbatch{1}.spm.stats.results.conspec(5).extent = 0;
        matlabbatch{1}.spm.stats.results.conspec(5).conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec(5).mask.none = 1;

        matlabbatch{1}.spm.stats.results.conspec(6).titlestr = 'MIX - Interaction';
        matlabbatch{1}.spm.stats.results.conspec(6).contrasts = 6;
        matlabbatch{1}.spm.stats.results.conspec(6).threshdesc = 'FWE';
        matlabbatch{1}.spm.stats.results.conspec(6).thresh = 0.1; % set alpha threshold here (p<0.05 or p<0.1)
        matlabbatch{1}.spm.stats.results.conspec(6).extent = 0;
        matlabbatch{1}.spm.stats.results.conspec(6).conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec(6).mask.none = 1;
    end
end
    

matlabbatch{1}.spm.stats.results.units = 2;
matlabbatch{1}.spm.stats.results.export{1}.ps = true;
matlabbatch{1}.spm.stats.results.export{2}.png = true;
matlabbatch{1}.spm.stats.results.export{3}.xls = true;
