%-----------------------------------------------------------------------
% Job saved on 18-Nov-2019 11:18:16 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% Select output folder (output files will be saved here)
results_dir = 'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_results\';

matlabbatch{1}.spm.stats.factorial_design.dir = {results_dir};

% Specify factors in the model
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Context';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 3;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0; % independent factor
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0; % assuming equal variance across all conds (ie. satisfies the test of sphericity), so no need to perform correction
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Ttype';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 3;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;


% Specify the images for each cell in the model
load([results_dir 'conds.mat']);

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = cond_1_1';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = cond_1_2';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [1 3];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = cond_1_3';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = cond_2_1';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).levels = [2 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).scans = cond_2_2';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).levels = [2 3];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).scans = cond_2_3';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).levels = [3 1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).scans = cond_3_1';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).levels = [3 2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).scans = cond_3_2';

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(9).levels = [3 3];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(9).scans = cond_3_3';
 

% Other things for the model (default settings)                                                               
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%{
function images = genImageList(image_dir, ROI_name, cond_name)
    images = [];
    NumSubjects = 24;
    
    for i = 1:NumSubjects
        image_file = [image_dir '\' ROI_name '_' cond_name '_subj' int2str(i) '_evoked\condition_Undefined.nii,1'];
        images{i} = image_file;
    end
end
%}