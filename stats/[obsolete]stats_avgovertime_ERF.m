%% First, run the top 2 sections in stats_ERF.m




%% Then, do the following MANUAL STEPS (1-3)

% (1) CHOOSE which stat output file to load
load([ResultsFolder_thisrun '\stats_Interactions_minnbchan2.mat'])

% (2) SPECIFY which contrast to look at:
stat = MixCost_interaction;
stat_name = 'MixCost_interaction'; % to use in output file name
%conds_to_plot = [1 3 4 6 7 9];

% read out the sig channels & time points
[row,col] = find(stat.mask);
sig_chans = unique(row)';
sig_samples = unique(col)';

% display the time points in the cluster
fprintf('\nTime points in the cluster:'); 
stat.time(sig_samples)

% (3) check the .time field in "allSubjects_erf", to find the sample indices corresponding to this time interval
start_sample = 192; %193; % temporal span of cluster for minnbchan=3
end_sample = 201;  %197;


%% automatic processing from here onwards

% extract a single value for each subject in each condition
% (avg over sig channels & time points)

%load([ResultsFolder_thisrun 'allSubjects_erf.mat']); 

T = table(); % create an empty table (to store the output)

%GFP_all = []; % for "diagnostic plot"

cond_names = fieldnames(allSubjects_erf);
for j = 1:length(cond_names) % each cycle handles one cond (e.g. NatStay)
    cond_name = cond_names{j};
    allsubjects = allSubjects_erf.(cond_name); % all subjects' time courses for this cond
    
    for i = 1:length(allsubjects) % each cycle handles one subject
        
        % Paul: if using axial ERFs (contains +ve & -ve values), 
        % it's not advisable to take a plain avg over the sig channels
        % (due to sign-flipping issues). You need to use the GFP or RMS.
        % If you planar transform the data (makes everything +ve),
        % and run all sensor-space analysis after that,
        % then you can just take plain avg across channels here.
        
        % In the end, we decided to take plain avg for axial ERFs anyway,
        % because GFP is not usable - it shows discrepancy due to 
        % diff processing order (see notes & diagnostic plot below).
        
        if 1 %(PLANAR_TRANSFORM) % take plain avg of sig channels
            cfg         = [];
            cfg.channel = stat.label(sig_chans); % Always retrieve channels using the "label" field (which uses actual channel names)
                  % rather than directly indexing into the matrix. This way we can ensure that we get the correct channels.
            cfg.avgoverchan = 'yes';
            
            avg_across_chans = ft_selectdata(cfg, allsubjects{i}); % take plain avg across all channels in the cluster -> this results in a single timecourse
            avg = mean(avg_across_chans.avg(start_sample:end_sample)); % avg over all time points in the cluster -> this results in a single value
            
            % old method - prone to selecting wrong channels (see new method above)
            %{
            subject_timecourse = allsubjects{i}.avg; % chan * time
            subject_timecourse_part = subject_timecourse(sig_chans, start_sample:end_sample); % extract the sub-matrix which contains the sig channels & time points
            avg_across_chans = mean(subject_timecourse_part); % avg across channels -> this results in a single timecourse
            avg = mean(avg_across_chans); % avg across time points -> this results in a single value
            %}
        else % calculate the GFP of all sig channels
            cfg         = [];
            cfg.method  = 'power';
            cfg.channel = stat.label(sig_chans);        
            GFP = ft_globalmeanfield(cfg, allsubjects{i});
            
            % for "diagnostic plot":
            % store all subjects's GFP in this cond, then we can do GA plot below
            %GFP_all.(cond_name){i} = GFP;
            
            % avg over all time points in the cluster
            avg = mean(GFP.avg(start_sample:end_sample));
        end
        
        T.(cond_name)(i,1) = avg; % store into the appropriate cell in the table
    end   
end

% Diagnostic plot:
%{
cfg = [];
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no'; % average across subjects
for j = 1:length(eventnames_real)
    GA_erf_GFP.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, GFP_all.(eventnames_real{j}){:});  
end

figure('Name','GFP_all_subjects'); hold on
for j = 1:length(eventnames_real)
    if colours_and_lineTypes % use a combination of colours and line types to distinguish conds
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours{j}, 'LineStyle', lineTypes{j});
    else % just use colours
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours(j,:));
    end
    xlim(PLOT_XLIM);
end
legend(eventnames_real);
%}
% This showed that we did not make a mistake. The discrepancy btwn the avgovertime results here
% & the overall GFP plot is prob due to diff processing order (GA then GFP vs. GFP on each subject then GA).
            

%% save to csv file
start_time = allSubjects_erf.NatStay{1}.time(start_sample) * 1000;
end_time   = allSubjects_erf.NatStay{1}.time(end_sample) * 1000;
output_file = [ResultsFolder_thisrun 'avgovertime\\' stat_name '_' int2str(start_time) '-' int2str(end_time) 'ms.csv'];
writetable(T, output_file);
