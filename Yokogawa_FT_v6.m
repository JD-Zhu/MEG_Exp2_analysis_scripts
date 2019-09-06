%{
Version log:

v3: latest stable working version with all the features
v4: most EVAL statements are replaced by using dynamic field names rather
than dynamic variable names (script should run faster)
v5: added baseline correction after calc'ing erf -> later removed, in order to preserve cov matrix;
    added error-trial rejection, based on "errorsheet" from beh data checking;
    added downsampling (for saving);
    added computation of cov matrix in timelockanalysis (to enable creation of spatial filters in sourceanalysis)
v6: split whole script into 3 main stages of processing & each stage now loops thru all subjects;
    abstracted certain sections into its own function;
    created "common.m" - equivalent to a .h file which can be included/executed at the top of each script
%}

%%
%close all
%clear all % disable this line if u want breakpoints to work

% run the #define section
global DataFolder; global ResultsFolder; global filename_suffix; 
global eventcodes; global eventnames; global eventnames_real; global collapse_across_langs;
global colours;
common();


% enable access to 'SubjectID' from inside "trig_fun_160_...", so that 
% correct code_delay value can be set for each subject (30ms for first 5 subjects, 100ms for the others)
global SubjectID; 

% find all subject folders containing raw MEG recording
SubjectIDs = [dir([DataFolder 'A*']); dir([DataFolder 'B*'])];
SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array
%SubjectIDs = {'B02-YW-3523'}; % or manually select which subjects to process


% === Settings ===

% Please adjust as required:

% > create a name for this run (this will be the name of the output folder)
%run_name = 'hpfilt0.5';
%run_name = 'raw';
run_name = 'CALM6001';
%run_name = 'TSPCA10000_3';
%run_name = 'detrend1000';
output_name = ['output\\' run_name '\\']; % location to save intermediate output files inside each SubjectFolder
ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % ERF results for all subjects
mkdir(ResultsFolder_thisrun);

% > name of confile in the SubjectFolder (can use wildcards)
confile_name = ['*' run_name '.con'];

% > which steps to run?
DO_HPF = false;
DO_DETREND = false;
DO_ICA = true; % if we used human subjects (i.e. not dry run), set this to true
RUN_ICA_ON_1HZ_FILTERED_DATA = false; % for MEG, we prob don't need to apply 1Hz HPF before running ICA
DO_BEH_CHECK = false; % if subjects produced beh responses, set this to true
DO_PCA = false; % if subjects produced vocal responses, set this to true

% when running many subjects in one batch, process all auto steps until the first manual step
RUN_UP_TO_DETRENDING = false;
RUN_UP_TO_MANUAL_ARTEFACT = false;  % 1st manual processing
RUN_UP_TO_ICA = false;               % 2nd manual processing

% > other options:
CHANNEL_REPAIR = false; % repair bad/rejected channels?
CALC_UNCLEANED_ERF = false; % calculate uncleaned erf? (for quality check of response-component rejection)

REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS = false; % remove trigger-leak artefact? (spike around 55ms before cue onset & target onset)
REMOVE_TRIGGER_ARTEFACT_ON_AVG_ERF = false;
REMOVE_TRIGGER_LEAK_CHANNELS = false; % remove all channels affected by trigger leak (81-89)?

% =================


% check the Settings, and modify stuff accordingly
%if (CHANNEL_REPAIR)
%    ResultsFolder = [ResultsFolder(1:end-2) '_channelrepair\\']; % modify ResultsFolder location
%end    

% set filenames for saving the output from each stage (so that we don't have to rerun the whole thing from beginning every time)
S1_output_filename = 'S1_preprocessed_data.mat'; % Stage 1 output (stored inside each Subject folder)
S2_output_filename = ['S2_after_visual_rejection' filename_suffix '.mat']; % Stage 2 output (stored inside each Subject folder)
S3_output_filename = ['_erf' filename_suffix '.mat']; % ERF output (stored in ResultsFolder for all subjects)

% load nececessary files
load('lay.mat');
load('neighbours.mat');
%load('all_labels.mat');

        

%% Stage 1: preprocessing & downsampling

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    % create output folder for all the save files
    output_path = [SubjectFolder output_name];
    mkdir(output_path);


    S1_output_file = [output_path S1_output_filename];
    
    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S1_output_file, 'file') ~= 2)          
        
        %{
        Preprocessing steps for Yokogawa MEG data:
        - filtering
        - manually mark bad sections & bad channels
        - ICA artefact rejection
        - trigger-based trial definition (i.e. epoching)
        %}
        
        
        % This is the combined con file (combined in MEG160),
        % so we only need to process one file now
        files = dir([SubjectFolder confile_name]);
        rawfile = [SubjectFolder files(1).name];

        
        % >>>
        % Step 1: filtering
        output_file = [output_path 'filtered.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            [alldata] = filtering(rawfile, DO_HPF);
            
            % save now, because the 0.1Hz HPF (blackman-windowed sinc FIR)
            % TAKES FOREVER to run (17 minutes for 1G con file)
            if (DO_HPF) % to save disk space, only save if we did HPF
                save(output_file, 'alldata', '-v7.3');
            end
        else
            load(output_file);
        end

        
        % >>>
        % Step 2: detrend (optional)
        if (DO_DETREND)
           output_file = [output_path 'after_detrend.mat'];

            % if haven't already processed this before, do it now & save a copy
            if (exist(output_file, 'file') ~= 2) 
                load('chanlocs.mat');
                [alldata] = detrend_on_continuous_data(alldata, chanlocs);
                
                save(output_file, 'alldata', '-v7.3');
            else
                load(output_file);
            end
        end
        
        % If running in batch, skip to next subject now
        if (RUN_UP_TO_DETRENDING)
            continue;
        end
        
        
        % Print out SubjectID so we know which subject we are working on
        fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

        
        % >>>
        % Step 3: manually mark artefacts
        output_file = [output_path 'arft.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            [arft] = mark_artefact(alldata);
            save(output_file, 'arft', '-v7.3');
        else
            load(output_file);
        end

        % reject the manually marked artefact
        arft.artfctdef.reject = 'nan'; % This fills those sections with NaNs.
                                      % You will then need to write code to
                                      % remove any trials containing NaNs
                                      % before computing ERF.
        alldata = ft_rejectartifact(arft, alldata);

        
        % >>>
        % Step 4: manually mark bad channels & reject them
        output_file = [output_path 'selChLabel.mat'];
        
        % Remove the trigger channels b4 identifying bad channels
        cfg         = [];
        cfg.channel = alldata.label(1:160);
        alldata     = ft_selectdata(cfg, alldata);
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            [selChLabel] = reject_bad_channels(alldata);
            save(output_file, 'selChLabel', '-v7.3');
        else
            load(output_file);
            % remove the bad channels
            cfg                         = [];
            cfg.channel                 = selChLabel;
            alldata                     = ft_selectdata(cfg, alldata);
        end

        % If running in batch, skip to next subject now
        if (RUN_UP_TO_MANUAL_ARTEFACT)
            continue;
        end        
        
        
        % >>>
        % Step 5: ICA
        if (DO_ICA)
            output_file = [output_path 'afterICA.mat'];
        
            % if haven't already processed this before, do it now & save a copy
            if (exist(output_file, 'file') ~= 2)  

                % Run ICA to identify components (if haven't done this yet)
                output_file_ICA = [output_path 'ICA_comps.mat'];        
                if (exist(output_file_ICA, 'file') ~= 2)    
                    
                    if (RUN_ICA_ON_1HZ_FILTERED_DATA) % apply 1Hz HPF before running ICA
                        [comp] = ICA_run(true, rawfile, arft, selChLabel);
                    else % directly run ICA without applying 1Hz HPF
                        [comp] = ICA_run(false, alldata);
                    end

                    % Immediately save a copy, as ICA takes a long time to run
                    save(output_file_ICA, 'comp', '-v7.3');
                else
                    load(output_file_ICA);
                end

                % If running in batch, skip to next subject now
                if (RUN_UP_TO_ICA)
                    continue;
                end


                % Select which ICA components to reject
                % set filename for diary (to record the components selected)
                set(0,'DiaryFile', [output_path 'ICA_log.txt']);

                [alldata_afterICA] = ICA_reject_comps(alldata, comp, lay);

                save(output_file, 'alldata_afterICA', '-v7.3');
            else
                load(output_file);
            end      
        else % if not doing ICA, just update the variable name
            alldata_afterICA = alldata;
        end
        
        
        % >>>
        % Step 6: Epoching
        output_file = [output_path 'epoched.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)  
            % Define trials using custom trialfun
            cfg                   = [];
            cfg.dataset           = rawfile;
            cfg.continuous        = 'yes';
            cfg.trialfun          = 'trig_fun_160_basic_v2';
            cfg.trialdef.prestim  = 0.8;      % pre-stimulus interval
            cfg.trialdef.poststim = 1;        % post-stimulus interval
            trialinfo_b = ft_definetrial(cfg);

            alldata_afterICA = ft_redefinetrial(trialinfo_b, alldata_afterICA);
            
            save(output_file, 'alldata_afterICA', '-v7.3');
        else
            load(output_file);
        end
        
        cfg         = [];
        cfg.demean  = 'yes'; % subtracts the mean of the time window from all samples (i.e. centres the waveform on 0)
        all_blocks = ft_preprocessing(cfg, alldata_afterICA);

        
        % >>>
        % Step 7: downsample the data for saving
        %all_blocks.time(1:end) = all_blocks.time(1); % this avoids numeric round off issues in the time axes upon resampling
        cfg            = [];
        cfg.resamplefs = 200; % sampling freq was 1000Hz, best to use a divisor of it (200Hz is commonly used)
        cfg.detrend    = 'no';
        all_blocks     = ft_resampledata(cfg, all_blocks);

        % SAVE preprocessed data - takes a while!!
        save(S1_output_file, 'all_blocks', 'trialinfo_b', '-v7.3');
    end
end


%%  Stage 2: trial exclusions

for k = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(k));
    SubjectFolder = [DataFolder SubjectID '\\'];

    output_path = [SubjectFolder output_name];
    S2_output_file = [output_path S2_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S2_output_file, 'file') ~= 2)    

        load([output_path S1_output_filename]);
   
        
        % divide up the master event list, to create 1 list for each cond
        events_allBlocks = identify_event_types(SubjectID, trialinfo_b);
        
        % in each list, remove the indices corresponding to error trials
        if (DO_BEH_CHECK)
            events_allBlocks = exclude_beh_errors(SubjectID, events_allBlocks);
        end

    
        % === PCA artefact removal ===
        % remove mouth-movement artefact by extracting main components from the "response" epochs
        % and projecting these out of all trials

        if (DO_PCA)
            [all_blocks_clean, response_comp] = remove_artefact_PCA(all_blocks, events_allBlocks, lay, 'response');
            %TODO: Based on B11-pilot, it seems that removing top 5 comps is too
            % much. The pairwise comparison plots show that this removes the N1 peak etc
            % maybe change to 1:3 comps or sth, be more conservative.
        
        else % if not doing ICA, just update the variable name
            all_blocks_clean = all_blocks;
        end
        
        
        % remove trigger-leak artefact (if needed)
        %if (REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS)
        %    load('lay.mat');
        %    [all_blocks_clean, trigger_comp] = remove_artefact_PCA(all_blocks_clean, events_allBlocks, lay, 'trigger');
        %end        
       
        
        % === Reject Outlier Trials ===

        % Print out SubjectID so we know which subject we are working on
        fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

        % Display visual trial summary to reject outlier trials
        cfg              = [];
        cfg.feedback     = 'no'; % suppress console output (so that it's easy to find the SubjectID we printed out above)
        cfg.method       = 'summary';
        cfg.metric       = 'zvalue'; % default is 'var'
        cfg.keepchannel  = 'no';
        cfg.keeptrial    = 'nan'; % we keep the rejected trials as 'NaN' here,
            % because if we remove them, that will change the indices of all subsequent trials,
            % which will no longer match the indices we are using in events_allBlocks
        all_blocks_clean = ft_rejectvisual(cfg, all_blocks_clean);
        
    
        save([output_path S2_output_filename], 'all_blocks_clean', 'events_allBlocks', '-v7.3'); 
        % 'all_blocks' was not changed in Stage 2, so don't need to save again
    end
end


%% Stage 3: time-domain analysis (i.e. compute erf)

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

        load([output_path S2_output_filename]);
        
        % remove trigger-leak artefact (if needed)
        if (REMOVE_TRIGGER_ARTEFACT_ON_INDI_EPOCHS)
            load('lay.mat');
            [all_blocks_clean, trigger_comp] = remove_artefact_PCA(all_blocks_clean, events_allBlocks, lay, 'trigger');
        end

        if (REMOVE_TRIGGER_LEAK_CHANNELS)
            cfg         = [];
            cfg.channel = {'all', '-AG083', '-AG087', '-AG088', '-AG082', '-AG084', '-AG086', '-AG081', '-AG085', '-AG089'}; % {'MEG'};
            all_blocks_clean = ft_selectdata(cfg, all_blocks_clean);
        end

        % perform channel repair if needed
        if (CHANNEL_REPAIR)
            load('neighbours.mat');
            all_labels = all_blocks_clean.cfg.channel; % full list of 160 labels
            all_blocks_clean = repair_bad_channels(all_blocks_clean, neighbours, all_labels);
        end
        
        
        % === ft_redefine all event types (i.e. 8 real conditions + 'response' event) ===

        % in uncleaned data
        if (CALC_UNCLEANED_ERF)
            % load 'all_blocks'
            load([output_path S1_output_filename]);

            for j = 1:length(eventnames)
                cfg = [];
                cfg.trials = events_allBlocks.(eventnames{j});
                trials.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks);
            end
        end
        
        % in cleaned data
        for j = 1:length(eventnames)
            cfg = [];
            cfg.trials = events_allBlocks.(eventnames{j});
            trials_clean.(eventnames{j}) = ft_redefinetrial(cfg, all_blocks_clean);
        end

                
        % === Remove trials containing NaN ===

        % Remove any trials containing NaN sections, otherwise we would end
        % up with partial trials.
        % Note 1: we remove these trials from both the DATA (all_blocks) 
        % and the TRIAL INDEX LIST (trial_info_b.event),
        % to keep them consistent with each other
        % Note 2: this step must be done after exclude_beh_errors,
        % otherwise the trial list here won't match up with the errorsheet
        
        %[all_blocks_clean, clean_event] = remove_nan_trials(all_blocks_clean, trialinfo_b.event);
        %trialinfo_b.event = clean_event;        
        
        % in uncleaned data
        if (CALC_UNCLEANED_ERF)
            for j = 1:length(eventnames)
                [trials.(eventnames{j}), ~] = remove_nan_trials(trials.(eventnames{j}), events_allBlocks.(eventnames{j}));
            end
        end
        
        % call this fn for each cond separately
        for j = 1:length(eventnames)
            [trials_clean.(eventnames{j}), events_allBlocks.(eventnames{j})] = remove_nan_trials(trials_clean.(eventnames{j}), events_allBlocks.(eventnames{j}));
        end
        
        
        % === collapse across langs (if needed) ===
        if collapse_across_langs
            % combine trials_clean(1 & 2), trials_clean(3 & 4), etc.
            for j = 1:length(eventnames_real)
                %fprintf(['j = ' num2str(j) ';  eventnames_real{j} = ' eventnames_real{j} '\n']);
                % copy over the general struct
                trials_clean_collapsed.(eventnames_real{j}) = trials_clean.(eventnames{j*2-1});
                % append these two trial lists together: j*2-1, j*2
                trials_clean_collapsed.(eventnames_real{j}).trial = [trials_clean.(eventnames{j*2-1}).trial, trials_clean.(eventnames{j*2}).trial];
                trials_clean_collapsed.(eventnames_real{j}).time = [trials_clean.(eventnames{j*2-1}).time, trials_clean.(eventnames{j*2}).time];
                trials_clean_collapsed.(eventnames_real{j}).sampleinfo = vertcat(trials_clean.(eventnames{j*2-1}).sampleinfo, trials_clean.(eventnames{j*2}).sampleinfo);
            end
            
            trials_clean = trials_clean_collapsed;
        end
        
        
        
        % delete vars that are no longer needed (to clear up some memory)
        %clear all_blocks;
        %clear all_blocks_clean;
        %clear trialinfo_b;

        %save([output_path 'before_computing_erf.mat'], 'trials', 'trials_clean', 'response_comp', '-v7.3');


        % === compute ERFs ===

        % in uncleaned data (just for quality check of PCA component rejection)
        if (CALC_UNCLEANED_ERF)
            for j = 1:length(eventnames_real)
                cfg         = [];
                %cfg.nanmean = 'yes';
                %trials.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
                erf.(eventnames_real{j}) = ft_timelockanalysis(cfg, trials.(eventnames_real{j})); % Do this to create average field and timelock struct
            end
        end
        
        % in cleaned data (compute erfs & cov matrices)
        [erf_clean, erf_allconds] = compute_ERF(trials_clean);

        % Do not perform baseline correction here, because it will remove the cov matrix from the timelock output:
        % "baseline correction invalidates previous covariance estimate, removing cov"
        %{
        % baseline correction
        for j = 1:length(eventnames)
            cfg = [];
            cfg.baseline = [-0.1 0];
            erf.(eventnames{j}) = ft_timelockbaseline(cfg, erf.(eventnames{j})); % uncleaned data
            erf_clean.(eventnames{j}) = ft_timelockbaseline(cfg, erf_clean.(eventnames{j})); % cleaned data
        end
        %}
        
        % remove trigger artefact on averaged ERF
        if REMOVE_TRIGGER_ARTEFACT_ON_AVG_ERF
            % here we use remove_artefact_PCA just to compute the trigger_comp
            % (we ignore the cleaned data it returns),
            % then manually call ft_rejectcomponent on the avg ERF for each cond
            load('lay.mat');
            [~, trigger_comp] = remove_artefact_PCA(all_blocks_clean, events_allBlocks, lay, 'trigger');
            
            % ft_rejectcomponent automatically converts timelock data to raw (and then
            % converts back) if you pass in one condition at a time.
            % the "cov" and "var" fields get removed, but we only need "cov" from 
            % erf_cue_combined & erf_target_combined anyway, not from the indi conditions.

            for j = 1:length(eventnames)
                cfg              = [];
                cfg.demean       = 'no';
                cfg.component    = 1:1; % project out the 1st principal component
                erf_clean.(eventnames{j}) = ft_rejectcomponent(cfg, trigger_comp, erf_clean.(eventnames{j}));
            end

            % Note: the covmatrix computed above (in compute_ERF.m) 
            % should still be valid, because it is not based on any data 
            % containing the trigger artefact (we modified cfg.covariancewindow
            % to 0~650ms following cue onset)
        end
        
        %### TODO: insert TF analysis here ###
        
        
        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', 'run_name', ...
            'trials_clean', 'erf_clean', 'erf_allconds', '-v7.3'); %'erf',       
    else
        load(S3_output_file);
    end

    
    % === Plot ERF & GFP (can use this to regen all plots from saved erf results) ===
    
    if (CALC_UNCLEANED_ERF)
        %plot_ERF(erf, erf_clean, [], lay, true, true, false);
    else % clean erf only
        plot_ERF([], erf_clean, erf_allconds, lay, false, true, true);
    end    
    
end
