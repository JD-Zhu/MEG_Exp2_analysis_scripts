%-----------------------------------------------------------------------
% Job saved on 06-Dec-2019 18:00:32 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% load the results_dir variable - to ensure the path we use in this script 
% matches the output folder specified in SPM_batch.m
load results_dir

matlabbatch{1}.spm.stats.con.spmmat = {[results_dir 'SPM.mat']};


% To construct the appropriate contrast matrices, simply use:
% Con = spm_make_contrasts([3 3]); % for a 3x3 design

matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Intercept (overall effect)';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 1 1 -1 -1 -1 0 0 0
                                                        0 0 0 1 1 1 -1 -1 -1]; % tmply replaced using Main effect of A
%matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 1 1 1 1 1 1 1 1];  % not working for some reason
%matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 1 1 1 1 1 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333 0.333333333333333];
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'Main effect of Factor A (Context)';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.weights = [1 1 1 -1 -1 -1 0 0 0
                                                        0 0 0 1 1 1 -1 -1 -1];  % we don't need to put in zeros(1,24) - SPM will autoly fill the 24 trailing cells with 0
%matlabbatch{1}.spm.stats.con.consess{2}.fcon.weights = [1 -1 0 0 0 0 0.333333333333333 0.333333333333333 0.333333333333333 -0.333333333333333 -0.333333333333333 -0.333333333333333 0 0 0
%                                                        0 1 -1 0 0 0 0 0 0 0.333333333333333 0.333333333333333 0.333333333333333 -0.333333333333333 -0.333333333333333 -0.333333333333333];
matlabbatch{1}.spm.stats.con.consess{2}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = 'Main effect of Factor B (Ttype)';
matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = [1 -1 0 1 -1 0 1 -1 0
                                                        0 1 -1 0 1 -1 0 1 -1];
%matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = [0 0 0 1 -1 0 0.333333333333333 -0.333333333333333 0 0.333333333333333 -0.333333333333333 0 0.333333333333333 -0.333333333333333 0
%                                                        0 0 0 0 1 -1 0 0.333333333333333 -0.333333333333333 0 0.333333333333333 -0.333333333333333 0 0.333333333333333 -0.333333333333333];
matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{4}.fcon.name = 'Interaction A*B (proper 3x3)';
matlabbatch{1}.spm.stats.con.consess{4}.fcon.weights = [1 -1 0 -1 1 0 0 0 0
                                                        0 1 -1 0 -1 1 0 0 0
                                                        0 0 0 1 -1 0 -1 1 0
                                                        0 0 0 0 1 -1 0 -1 1];
%matlabbatch{1}.spm.stats.con.consess{4}.fcon.weights = [0 0 0 0 0 0 1 -1 0 -1 1 0 0 0 0
%                                                        0 0 0 0 0 0 0 1 -1 0 -1 1 0 0 0
%                                                        0 0 0 0 0 0 0 0 0 1 -1 0 -1 1 0
%                                                        0 0 0 0 0 0 0 0 0 0 1 -1 0 -1 1];
matlabbatch{1}.spm.stats.con.consess{4}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 0;
