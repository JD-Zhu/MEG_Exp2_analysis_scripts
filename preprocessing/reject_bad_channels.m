% manually identify bad channels
% Note: this fn does not change the data (ie. it does not remove the bad channels), 
% it only produces a list containing the channels that we select to keep

function [selChLabel] = reject_bad_channels(alldata)
    cfg                         = [];
    cfg.method                  = 'channel';
    %cfg.method                  = 'summary';
    cfg.alim                    = 1e-10;
    cfg.keepchannel             = 'no';
    cfg.keeptrial               = 'nan';
    ft_rejectvisual (cfg, alldata); % visual tool to help identify bad channels

    selChLabel                  = alldata.label;
end