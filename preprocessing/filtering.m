% do_HPF: if yes, apply 0.5Hz high-pass filter

function [alldata] = preprocessing(rawfile, do_HPF)
    
    hdr = ft_read_header(rawfile, 'dataformat','yokogawa_con'); % read header file
         
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

    % highpass filter (optional)
    if (do_HPF)
        cfg            = [];
        cfg.hpfilter   = 'yes';
        cfg.hpfilttype = 'firws';
        cfg.hpfreq     = 0.5; % 0.5 +- 0.1Hz
        cfg.hpfiltdf   = 0.2; % transition window width (for firws; this param overrides order)
                              % hpfreq - (hpfiltdf / 2) must be >= 0
        cfg.hpfiltwintype = 'blackman';
        cfg.hpfiltdir  = 'onepass-zerophase';
        alldata = ft_preprocessing(cfg, alldata);
    end

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
    
end