%{
Performs preprocessing for Yokogawa MEG data:
- filtering
- manually mark bad sections & bad channels
- ICA artefact rejection
- trigger-based trial definition (i.e. epoching)
%}

function [all_blocks, trialinfo_b] = preprocessing(SubjectFolder, lay)
    
    % find all .con files
    %{
    files = dir([SubjectFolder '*.con']);
    % skip "_test.con" file
    for i = 1:length(files)
        if (strcmp(files(i).name(end-8:end), '_test.con') == 1)
            files(i) = [];
        end
    end
    %}
    
    % This is the combined con file (combined in MEG160)
    % so we only need to process one file now
    files = dir([SubjectFolder '*_B1-concat.con']);
    
    
    % create output folder for all the save files
    output_path = [SubjectFolder 'output\\hpfilt0.5\\'];
    mkdir(output_path);

    % each cycle processes one '.con' file.
    % we do it in multiple steps (multiple calls to ft_prepreocessing)
    % to ensure things happen in the desired order
    for i = 1:length(files)
        rawfile = [SubjectFolder files(i).name];
        hdr = ft_read_header(rawfile, 'dataformat','yokogawa_con'); % read header file

        output_file = [output_path 'filtered_B' num2str(i) '.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)            
            % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
            cfg                      = [];
            cfg.trialfun             = 'ft_trialfun_general';
            cfg.headerfile           = rawfile;
            cfg.datafile             = rawfile;
            cfg.trialdef.triallength = Inf;
            cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
            cfg = ft_definetrial(cfg);

            cfg.continuous = 'yes';
            alldata = ft_preprocessing(cfg);
            
            % Select the gradiometers (i.e. actual MEG channels)
            % plus the "all triggers" channel (193)
            cfg         = [];
            cfg.channel = [hdr.grad.label; '191'; '193']; % alldata.label(1:160);
            alldata = ft_selectdata(cfg, alldata);
            
            % ft_preprocessing: reads in MEG data
            %{
            cfg.bpfilter   = 'yes';
            cfg.bpfreq     = [0.2 30]; % bandpass filter [0.5 30], successfully filtered out the low-freq drift!!
            cfg.bpfiltwintype = 'blackman'; % may help to get rid of the ringing effect?
            cfg.bpfiltord = 64; % default is 4
            alldata = ft_preprocessing(cfg);
            %}        

            % highpass filter
% TEMPORARILY disabling this - we applied online filter 0.1Hz, so no need to do this again offline
%
            cfg            = [];
            cfg.hpfilter   = 'yes';
            cfg.hpfilttype = 'firws';
            cfg.hpfreq     = 0.5; % 0.5 +- 0.1Hz
            cfg.hpfiltdf   = 0.2; % transition window width (for firws; this param overrides order)
                                  % hpfreq - (hpfiltdf / 2) must be >= 0
            cfg.hpfiltwintype = 'blackman';
            cfg.hpfiltdir  = 'onepass-zerophase';
            alldata = ft_preprocessing(cfg, alldata);
%
            
            % lowpass filter
            cfg         = [];
            cfg.lpfilter   = 'yes';
            cfg.lpfilttype = 'firws';
            cfg.lpfreq     = 40; % 40 +- 10Hz
            cfg.lpfiltdf   = 20; % wider transition window means it will run much faster
            cfg.lpfiltwintype = 'blackman';
            cfg.lpfiltdir  = 'onepass-zerophase';
            alldata = ft_preprocessing(cfg, alldata);
                                       
            % deal with 50Hz line noise (necessary even after bandpass filter, coz the 50Hz noise is huge)
            cfg          = [];
            cfg.bsfilter = 'yes';
            cfg.bsfreq   = [49.5 50.5];
            % alternatively, can use: (but you need to pad the data to 5~10 seconds)
            % http://www.fieldtriptoolbox.org/faq/what_kind_of_filters_can_i_apply_to_my_data/
            %cfg.dftfilter = 'yes';
            %cfg.dftfreq   = [50 100 150];
            alldata = ft_preprocessing(cfg, alldata);
                        
            % save now, because the 0.1Hz HPF (blackman-windowed sinc FIR)
            % TAKES FOREVER to run (17 minutes for 1G con file)
            save(output_file, 'cfg', 'alldata', '-v7.3');
        else
            load(output_file);
        end;

        
        % browse raw data to remove noisy segments (e.g. clenched jaw)
        % http://www.fieldtriptoolbox.org/walkthrough/#visual-data-inspection
        %
        % Note: doing 'partial' rejection here because there is only 1 trial,
        % but this gives weird partial trials later on, and the final avg timecourse
        % (ERFs) are also partial - which is bad. 
        % Alternative is to do this step after epoching, however in that case
        % the data we inspect will be after filtering, so it's harder to
        % identify the artefacts (coz muscle artefacts are usually high freq),
        % plus the filtering itself will create some weird stuff which 
        % may be mistakenly identified as artefact
        %
        % Maybe best to use ICA later instead of doing manual selection here
        % see remove_artefact_ICA.m (called by Yokogawa_FT_v6.m)
        %
        output_file = [output_path 'arft_B' num2str(i) '.mat'];
        if (exist(output_file, 'file') ~= 2)    
            cfg           = [];
            cfg.viewmode  = 'vertical';
            cfg.continous = 'yes';
            cfg.blocksize = 120; % display 2-min segments
            cfg.ylim      = [ -4e-13   4e-13 ];
            % first use this scaling to see trigger channels (to mark breaks)
            cfg.blocksize = 300; % display 5-min segments
            cfg.ylim      = [ -0.25  0.25 ];
            cfg.channel   = alldata.label(end-1:end); % display channels [191 193];

            % randomly sample a few channels, to check if the bandpass filter 
            % can successfully filter out the low-freq drift
            %{
            %cfg.channel = [1 8 15 23 30 37 44 51 58 65 72 78 84 91 98 104 111 119 125 132 139 146 153]; 
            cfg.channel = [1 11 35 49 60 72 85 98 111 125 139 153]; 
            cfg.ylim = [ -1.0941e-12  1.0941e-12 ];
            %}
            arft = ft_databrowser(cfg, alldata);
            save(output_file, 'arft');
        else
            load(output_file);
        end;

        % reject the manually marked artefact
        arft.artfctdef.reject = 'nan'; % This fills those sections with NaNs.
                                      % You will then need to write code to
                                      % remove any trials containing NaNs
                                      % before computing ERF.
        alldata = ft_rejectartifact(arft, alldata);
        
        
        % Remove the trigger channels b4 identifying bad channels
        cfg         = [];
        cfg.channel = alldata.label(1:160);
        alldata     = ft_selectdata(cfg, alldata);
        
        % Mark bad channels manually 
        output_file = [output_path 'selChLabel_B' num2str(i) '.mat'];
        if (exist(output_file, 'file') ~= 2)            
            cfg                         = [];
            cfg.method                  = 'channel';
            %cfg.method                  = 'summary';
            cfg.alim                    = 1e-10;
            cfg.keepchannel             = 'no';
            cfg.keeptrial               = 'nan';
            alldata = ft_rejectvisual (cfg, alldata); % this removes the bad channels

            selChLabel                  = alldata.label;
            save(output_file, 'selChLabel');
        else
            load(output_file);
            % remove the bad channels
            cfg                         = [];
            cfg.channel                 = selChLabel;
            alldata                     = ft_selectdata(cfg, alldata);
        end;
 
 
        % === ICA artefact removal ===
        % remove eye blinks & other large muscle artefacts (e.g. jaw clenching, hand movements)
        
        % remove the 3 "ref channels"
        %{
        cfg         = [];
        cfg.channel = all_blocks.label([1:91 93:94 97:160]);
        all_blocks = ft_selectdata(cfg, all_blocks);
        %}
        output_file = [output_path 'data_afterICA_B' num2str(i) '.mat'];
        if (exist(output_file, 'file') ~= 2)  
            % need to run ICA on 1Hz-filtered data, so we prep that data now
            cfg                         = [];
            cfg.trialfun                = 'ft_trialfun_general';  
            cfg.channel                 = hdr.grad.label; 
            cfg.continuous              = 'yes';
            cfg.hpfilter                = 'yes';
            cfg.hpfilttype              = 'firws';
            cfg.hpfreq                  = 1;
            cfg.hpfiltdf                = 1.5;
            cfg.hpfiltwintype           = 'blackman';
            cfg.hpfiltdir               = 'onepass-zerophase';
            cfg.dftfreq                 = 50; % removal line noise
            cfg.headerfile              = rawfile;
            cfg.datafile                = rawfile;
            data4ICA                    = ft_preprocessing(cfg);

            %lowpass filter for ICA
            cfg                         = [];
            cfg.lpfilter                = 'yes';
            cfg.lpfilttype              = 'firws';
            cfg.lpfreq                  = 40;
            cfg.lpfiltdf                = 10;
            cfg.lpfiltwintype           = 'blackman';
            cfg.lpfiltdir               = 'onepass-zerophase';
            data4ICA                    = ft_preprocessing(cfg, data4ICA);

            %reject the artifacts and channels in data4ICA that have been marked in script of Step 4 
            arft.artfctdef.reject       = 'nan';
            data4ICA = ft_rejectartifact(arft, data4ICA);

            cfg                         = [];
            cfg.channel                 = selChLabel;
            data4ICA                    = ft_selectdata(cfg, data4ICA);
    
            
            % Set filename for diary (to record the components selected for rejection in ICA)
            set(0,'DiaryFile', [output_path 'ICA_log_B' num2str(i) '.txt']);
            
            % Now call fn to run ICA
            [alldata_afterICA, comp] = remove_artefact_ICA(alldata, lay, data4ICA);
            
            
            save(output_file, 'alldata_afterICA','-v7.3');
            save([output_path 'ICA_comps_B' num2str(i) '.mat'], 'comp','-v7.3');
        else 
            load(output_file);
        end
        
        
        % === Epoching ===
        % Define trials using custom trialfun
        cfg                   = [];
        cfg.dataset           = rawfile;
        cfg.continuous        = 'yes';
        cfg.trialfun          = 'trig_fun_160_basic_v2';
        cfg.trialdef.prestim  = 0.8;      % pre-stimulus interval
        cfg.trialdef.poststim = 1;        % post-stimulus interval
        trialinfo_b(i) = ft_definetrial(cfg);
        
        alldata_afterICA = ft_redefinetrial(trialinfo_b(i), alldata_afterICA);

        
        cfg         = [];
        cfg.demean  = 'yes'; % subtracts the mean of the time window from all samples (i.e. centres the waveform on 0)
        %cfg.detrend = 'yes'; % removes low-frequency drift
        
        %cfg.polyremoval = 'yes'; % == remove higher order polynomial (default: 2nd-order)
        %cfg.polyorder = 40; % customise the order (higher order => affects higher freqs?)
        %
        block(i) = ft_preprocessing(cfg, alldata_afterICA);
        
        % <2> locdetrend
        %{
        % first, convert to EEGLAB format
        data = block(i);
        output_file = [output_path 'epoched_B' num2str(i) '.mat'];
        save(output_file, 'data', '-v7.3');
        [EEG] = fieldtrip2eeglab(output_file);
        % fill in the chanlocs var
        temp = load([ResultsFolder 'chanlocs.mat']);
        EEG.chanlocs = temp.chanlocs;
        
        %EEG = robust_locdetrend(EEG, 'no', 300, 150, 'yes'); % ~0.83Hz HPF; works perfectly (epoches are completely flat after lcodetrend), but is this removing too much?
        %EEG = robust_locdetrend(EEG, 'no', 600, 300, 'yes'); % ~0.42Hz HPF; a little bit of change to the waveform, still some drifts remaining
        EEG = robust_locdetrend(EEG, 'no', 900, 450, 'yes'); % ~0.42Hz HPF; a little bit of change to the waveform, still some drifts remaining
        %EEG = robust_locdetrend(EEG, 'no', 1000, 500, 'yes'); % ~0.25Hz HPF; almost no change to the waveform
        
        block(i) = ft_preprocessing(cfg, EEG);      
        %}


        % Create layout file for later (plotting)
        %{
        cfg      = [];
        cfg.grad = alldata.grad; % struct containing gradiometer definition
        lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations
        %save([ResultsFolder 'lay.mat'], 'lay');
        %}
        
        % Prepare neighbours & save for use in stat analysis
        %{
        cfg_neighb        = [];
        cfg_neighb.method = 'triangulation'; % 'distance' may have some issues
                   % 'triangulation' seems to give each channel more neighbours 
                   % (about twice as many as 'distance' method gives)
        neighbours        = ft_prepare_neighbours(cfg_neighb, all_blocks);
        save([ResultsFolder 'neighbours.mat'], 'neighbours');
        %}

    end

    % combine data from all the blocks into one dataset
    all_blocks = block(1);
    for i = 2:length(block) % each cycle appends one block
        all_blocks = ft_appenddata([], all_blocks, block(i));
    end
    
end