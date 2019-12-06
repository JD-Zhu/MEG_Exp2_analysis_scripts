%-----------------------------------------------------------------------
% Job saved on 22-Nov-2019 10:35:09 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% Select folder - this should match the output folder selected in SPM_batch.m
results_dir = 'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3_freeori\SPM_results\';

matlabbatch{1}.spm.stats.factorial_design.dir = {results_dir};

% Specify factors in the model
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Context';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Ttype';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;


% Specify all the images (in this order: NatSingle, NatStay, NatSw,
% ArtSingle, ArtStay, ArtSw, BiSingle, BiStay, BiSw,
% then repeat the above for subj2, subj3 ...)
load([results_dir 'conds.mat']);
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans = conds;
    % Example of the required format:
    %{
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_ArtSingle_subj10_evoked\condition_Undefined.nii,1'
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_ArtSingle_subj11_evoked\condition_Undefined.nii,1'
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_ArtSwitch_subj8_evoked\condition_Undefined.nii,1'
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_ArtSwitch_subj9_evoked\condition_Undefined.nii,1'
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_BiSingle_subj10_evoked\condition_Undefined.nii,1'
    'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\LACC_BiSingle_subj11_evoked\condition_Undefined.nii,1'
    %}
      
% Specify the design matrix - see ref below:
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;ed07c72c.1602
% https://www.jiscmail.ac.uk/cgi-bin/wa.exe?A2=SPM;f35bc681.1203
load([results_dir 'design.mat']);
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = design_complete;
    % BU of simplified version (ie. not including the Subject column):
    %{
    [1 1
    1 2
    1 3
    2 1
    2 2
    2 3
    3 1
    3 2
    3 3];
    %}

    
% Other things for the model (default settings)                                                               
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 3;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{3}.inter.fnums = [2
                                                                                  3];
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
