% Run ICA to identify eye blinks & other large muscle artefacts 
% (e.g. jaw clenching, hand movements)
%
% @param filter_first - run ICA on 1Hz-filtered data (usually applicable to EEG data,
%                       which has less channels, e.g. 32-channel. If you don't filter at 1Hz first,
%                       slow waves will take up too many comps & you might not get a comp for eye blinks)
% if yes: then need to supply these arguments: rawfile, arft, selChLabel
%       e.g. [comp] = ICA_run(true, rawfile, arft, selChLabel);
% if no:  then need to supply this argument: alldata
%       e.g. [comp] = ICA_run(false, alldata);

function [comp] = ICA_run(filter_first, varargin)
    %% prep the data for ICA (apply 1Hz HPF if need to) 
    if (filter_first)
        rawfile    = varargin{1};
        artf       = varargin{2};
        SelChLabel = varargin{3};
        
        hdr = ft_read_header(rawfile, 'dataformat','yokogawa_con'); % read header file
    
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

        % low-pass filter for ICA
        cfg                         = [];
        cfg.lpfilter                = 'yes';
        cfg.lpfilttype              = 'firws';
        cfg.lpfreq                  = 40;
        cfg.lpfiltdf                = 20;
        cfg.lpfiltwintype           = 'blackman';
        cfg.lpfiltdir               = 'onepass-zerophase';
        data4ICA                    = ft_preprocessing(cfg, data4ICA);

        % reject the artifacts and channels in data4ICA that have been marked Steps 3 & 4 
        arft.artfctdef.reject       = 'nan';
        data4ICA                    = ft_rejectartifact(arft, data4ICA);

        cfg                         = [];
        cfg.channel                 = selChLabel;
        data4ICA                    = ft_selectdata(cfg, data4ICA);

    else % if not applying the 1Hz filter, just read in the data
        
        data4ICA = varargin{1};
    
    end        

    
    %% Run ICA (don't use 'fastica' method, use 'runica' method)
    disp('About to run ICA using the Runica method')
    cfg            = [];
    cfg.method     = 'runica';
    cfg.channel    = 'all'; 
    comp           = ft_componentanalysis(cfg, data4ICA);

end