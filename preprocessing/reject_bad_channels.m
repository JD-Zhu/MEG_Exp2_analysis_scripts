% manually mark bad channels & reject them

function [selChLabel] = reject_bad_channels(alldata)
    cfg                         = [];
    cfg.method                  = 'channel';
    %cfg.method                  = 'summary';
    cfg.alim                    = 1e-10;
    cfg.keepchannel             = 'no';
    cfg.keeptrial               = 'nan';
    alldata = ft_rejectvisual (cfg, alldata); % this removes the bad channels

    selChLabel                  = alldata.label;
end