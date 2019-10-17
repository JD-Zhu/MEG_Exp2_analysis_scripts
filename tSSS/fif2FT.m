%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fif2FT.m: loads fif files into FieldTrip
%
% This function should be used after running mq_con2fif_maxfilter.py.
% It reads in the fif files (output from tSSS) and converts them into a
% fieldtrip-compatible data format.
%
% NOTE: this script relies on the default file names generated by tSSS. 
% Therefore, please do not rename the fif files!
%
% Author: Judy Zhu (di.zhu@mq.edu.au)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
%   dir_name      =  the folder containing the con file & fif files
%   file_name     =  name of the con file
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
%   alldata_tsss  =  data from the fif files, in a fieldtrip-compatible format
%
%%%%%%%%%%%%%%%%%
% Example usage:
%%%%%%%%%%%%%%%%%
%
%   alldata_tsss = fif2FT('E:\Judy\MEG-data\', '3591_KW_ME180_2019_09_28.con')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alldata_tsss = fif2FT(dir_name, file_name)
    %% Checking input

    % make sure the folder location is valid
    if isempty(dir(dir_name))
        fprintf('\nError: Directory does not exist. Please check the dir_name you specified.\n');
    end
    
    % remove the '.con' at the end of file_name
    if strcmp(file_name(end-3:end), '.con')
        file_name = file_name(1:end-4);
    end
    
    fullpath = fullfile(dir_name, file_name); % 'E:\Judy\MEG-data\3591_KW_ME180_2019_09_28';
    
    fprintf('\nLoading fif files and converting to FieldTrip format...\n\n');
    

    %% Read in the original data (i.e. con file), to create the "alldata" structure

    rawfile = [fullpath '.con'];            

    cfg                      = [];
    cfg.trialfun             = 'ft_trialfun_general';
    cfg.headerfile           = rawfile;
    cfg.datafile             = rawfile;
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials     = 1; % read in all data as a single segment
    cfg = ft_definetrial(cfg);

    cfg.continuous = 'yes';
    alldata = ft_preprocessing(cfg);


    %% Use the data after tSSS to replace the original data
    % We can read in the fif files using ft_read_data. This gives us a matrix 
    % which contains the actual MEG measurements, so we just need to 
    % insert this matrix into the "alldata" struct created earlier

    % If we have a large dataset (i.e. > 2087 seconds of MEG recording), 
    % the tSSS script should output more than one .fif files, 
    % so we need to join these together

    % check how many fif files there are
    fif_files = dir([fullpath '_raw_tsss*.fif']);
    fif_files = {fif_files.name}; % extract the filenames

    % read in the first fif file
    tsss = ft_read_data([fullpath '_raw_tsss.fif']);
    tmp = tsss(1:160, :); % only take the MEG channels (1-160)

    % Read in each subsequent fif file
    for i = 2:length(fif_files)
        tsss = ft_read_data([fullpath '_raw_tsss-' mat2str(i-1) '.fif']);
        tmp = [tmp tsss(1:160, :)]; % append to the previous matrix
    end

    % use the tSSS data to replace the original data
    alldata_tsss = alldata;
    alldata_tsss.trial{1,1}(1:160, :) = tmp;


    %% Plot
    cfg           = [];
    cfg.viewmode  = 'vertical';
    cfg.continous = 'yes';
    cfg.blocksize = 60; % display 60-sec segments
    cfg.ylim      = [ -4e-13   4e-13 ];
    cfg.channel = 'MEG';

    ft_databrowser(cfg, alldata);       % before tSSS
    ft_databrowser(cfg, alldata_tsss);  % after tSSS
    
end