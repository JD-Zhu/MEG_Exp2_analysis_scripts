% Computes ERFs & covariance matrices for by-condition data, combined cue & combined target

function [erf_clean, erf_allconds] = compute_ERF (trials_clean)
    % run the #define section
    %global conds_cue; global conds_target; 
    global eventnames_real;
    common();
    
    % We kept bad trials as NaN, so need to exclude them now
    for j = 1:length(eventnames_real)
        %length(trials_clean.(eventnames_real{j}).trial)
        good_trials_idx = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),trials_clean.(eventnames_real{j}).trial,'uni',0)))); % Just need to evaluate the first element as all samples in bad trial are NaN
        cfg             = [];
        cfg.trials      = good_trials_idx;
        trials_clean.(eventnames_real{j}) = ft_redefinetrial(cfg, trials_clean.(eventnames_real{j}));
    end

    % An alternative method to remove NaN trials - doesn't seem to work!
    %{
    for j = 1:length(eventnames_real)
        cfg         = [];
        cfg.nanmean = 'yes';
        trials_clean.(eventnames_real{j}) = ft_selectdata(cfg, trials_clean.(eventnames_real{j})); % Do this because we kept bad trials as NaN
    end
    %}


    % Compute erf & cov matrix on the combined data (all conds appended together)
    cellarray = struct2cell(trials_clean); % convert struct to cell array, 
    trials_clean_allconds = ft_appenddata([], cellarray{:}); % then you can feed it in as 'varargin'

    cfg                  = [];
    cfg.covariance       = 'yes';
    cfg.covariancewindow = [0 0.5]; % do not include any period after vocal response onset
    erf_allconds  = ft_timelockanalysis(cfg, trials_clean_allconds);

    % Alternative method (not sure if this is 100% correct)
    %{
    % average ERF across all conds
    erf_allconds = erf_clean.(eventnames_real{1});
    for j = 2:length(eventnames_real) % add together the ERF from all conds
        erf_allconds.avg = erf_allconds.avg + erf_clean.(eventnames_real{j}).avg;    
    end
    erf_allconds.avg = erf_allconds.avg / length(eventnames_real); % find the mean
    %}
    
    % Compute erf & cov matrix for each condition 
    for j = 1:length(eventnames_real)
        erf_clean.(eventnames_real{j}) = ft_timelockanalysis(cfg, trials_clean.(eventnames_real{j}));
    end

end