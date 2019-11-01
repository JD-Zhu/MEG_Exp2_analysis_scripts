% Script to compute spatial filter that evolves over time
% using the EMSf toolbox:
% https://bitbucket.org/emsf/emsf_matlab/src/default/
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

%%%%%% USER INPUT REQUIRED %%%%%%

% Specify the run_name here:
run_name = 'TSPCA10000_3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% run the #define section
global DataFolder; global ResultsFolder; global filename_suffix; 
global eventnames_real; global collapse_across_langs;
global colours;
common();

% where to read in the processed data for each subject
output_name = ['output\\' run_name '\\']; % intermediate output files inside each SubjectFolder
S2_output_filename = ['S2_after_visual_rejection' filename_suffix '.mat']; % Stage 2 output (stored inside each Subject folder)

% where to save the output from EMSf
ResultsFolder_thisrun = [ResultsFolder run_name '_EMSf' '\\']; % ERF results for all subjects
mkdir(ResultsFolder_thisrun);
S3_output_filename = ['_erf' filename_suffix '.mat'];

% add the EMSf toolbox to Matlab search path
% DO NOT use genpath(), only add the root dir. Otherwise it causes name conflict with the 'nearest' fn in FT
addpath('C:\Users\43606024\Documents\MATLAB\emsf-emsf_matlab-75c285db6472'); 
        

%% Start

% find all subject folders containing raw MEG recording
SubjectIDs = [dir([DataFolder 'A*']); dir([DataFolder 'B*'])];
%SubjectIDs([14]) = []; % remove certain subjects from the list
SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array

% loop thru all subjects
for i = 1:length(SubjectIDs)

    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    output_path = [SubjectFolder output_name];
    S3_output_file = [ResultsFolder_thisrun SubjectID S3_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S3_output_file, 'file') ~= 2)    
        
        % make sure we have a clean start (i.e. no leftover var contents from last subject)
        clear erf; clear erf_clean;
        clear trials; clear trials_clean;

        % load the S2_output_file
        load([output_path S2_output_filename]);


        % SELECT which 2 conds to compare
        %TEMP - just looking at stay vs sw in bivalent context for now:
        stay_trials = [events_allBlocks.BiStayC; events_allBlocks.BiStayE];
        switch_trials = [events_allBlocks.BiSwitchC; events_allBlocks.BiSwitchE];

        % one column for each cond, using logical indices to indicate which
        % trials belong to this cond
        conds = zeros(length(all_blocks_clean.trial), 2); % 2 columns = 2 conditions
        conds(stay_trials, 1) = 1;   % column 1: stay trials
        conds(switch_trials, 2) = 1; % column 2: switch trials

        % convert cell array of trials {channel x time} => 3d matrix (channel x time x trial)
        data = convert_FT_to_eeglab(all_blocks_clean.trial, 'EMSf'); 

        % call EMS-filtering
        [trl, msf] = ems_2cond_diff(data, conds, []);


        %% Average across trials to compute ERF
        trials_clean = [];

        cfg = [];
        cfg.trials = stay_trials;
        trials_clean.Stay = ft_redefinetrial(cfg, all_blocks_clean);
        cfg = [];
        cfg.trials = switch_trials;
        trials_clean.Switch = ft_redefinetrial(cfg, all_blocks_clean);
        %{
        for j = 1:length(eventnames_real)
            cfg = [];
            cfg.trials = events_allBlocks.(eventnames_real{j});
            trials_clean.(eventnames_real{j}) = ft_redefinetrial(cfg, all_blocks_clean);
        end
        %}

        % remove the EMS-filtering toolbox from search path,
        % otherwise the 'nearest' fn in this toolbox will cause an issue for
        % compute_ERF (which is meant to use the FT 'nearest' fn)
        %rmpath('C:\Users\43606024\Documents\MATLAB\emsf-emsf_matlab-75c285db6472');

        [erf_clean, erf_allconds] = compute_ERF(trials_clean);   

        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', 'run_name', ...
            'erf_clean', 'erf_allconds', '-v7.3');
    else
        load(S3_output_file);
    end

    
    % === Plot ERF & GFP (can use this to regen all plots from saved erf results) ===
    load('lay.mat');
    plot_ERF_general([], erf_clean, erf_allconds, lay, false, true, false);

    % save the plot as an image
    png_output_file = [ResultsFolder_thisrun 'Figures_GFP\\' SubjectID '_GFP.png'];
    if (exist(png_output_file, 'file') ~= 2) 
        saveas(gcf, png_output_file);
    end


    %% TODO %%
    % 1. Run stats across all subjects
    
    
    % 2. we can also plot the spatial filter using ems_topo_explore.
    % More things to do? see:
    % https://bitbucket.org/emsf/emsf_matlab/src/default/readme.txt
        

end
