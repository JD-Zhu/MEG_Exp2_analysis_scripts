% load epoched data for one condition
load('E:\Judy\Exp2\6_MEG-data\Results_ERF\TSPCA10000_3\A01-XC-3489_erf.mat')
a = trials_clean.NatStay;
%save a a

% prepare S
S.D = spm_eeg_ft2spm(a, 'temp');
S.mode = 'scalp x time';

% convert 2 SPM images
[images, outroot] = spm_eeg_convert2images(S);


%{
% interpolate missing channels
load('neighbours.mat');
load('all_labels.mat');
a_interp = repair_bad_channels(a, neighbours, all_labels);
%}