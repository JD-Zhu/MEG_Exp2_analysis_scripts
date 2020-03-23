%% First, run the top 3 sections in stats_ERF.m
% (to load the data)



%% Then, run the following code to extract a single value for each 
% subject in each condition (avg over sig channels & time points)


% SPECIFY which contrast to look at:
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
% MANUAL STEP - check the .time field in "allSubjects_erf", to find the sample indices corresponding to this time interval
start_sample = 192; %193; % temporal span of cluster for minnbchan=3
end_sample = 201;  %197;


T = table(); % create an empty table (to store the output)

cond_names = fieldnames(allSubjects_erf);
for j = 1:length(cond_names) % each cycle handles one cond (e.g. NatStay)
    cond_name = cond_names{j};
    allsubjects = allSubjects_erf.(cond_name); % all subjects' time courses for this cond
    
    for i = 1:length(allsubjects) % each cycle handles one subject
        
        % Paul: It's not advisable to take a plain avg over the sig channels here
        % (due to sign-flipping issues). You need to use the GFP or RMS.
        % Or you can planar transform the data (makes everything +ve),
        % then run all sensor-space analysis after that
        
        % Opt 1: calculate the GFP of all sig channels
        cfg         = [];
        cfg.method  = 'power';
        cfg.channel = stat.label(sig_chans);        
        GFP = ft_globalmeanfield(cfg, allsubjects{i});
        % avg over all time points in the cluster
        avg = mean(GFP.avg(start_sample:end_sample));
        
        % Opt 2: take a plain avg over all sig channels <- not advisable
        %{
        subject_timecourse = allsubjects{i}.avg; % chan * time
        subject_timecourse_part = subject_timecourse(sig_chans, start_sample:end_sample); % extract the sub-matrix which contains the sig channels & time points
        avg_across_chans = mean(subject_timecourse_part); % avg across channels -> this results in a single timecourse
        avg = mean(avg_across_chans); % avg across time points -> this results in a single value
        %}
        
        T.(cond_name)(i,1) = avg; % store into the appropriate cell in the table
    end   
end

%% save to csv file
start_time = allSubjects_erf.NatStay{1}.time(start_sample) * 1000;
end_time   = allSubjects_erf.NatStay{1}.time(end_sample) * 1000;
output_file = [ResultsFolder_thisrun 'avgovertime\\' stat_name '_' int2str(start_time) '-' int2str(end_time) 'ms.csv'];
writetable(T, output_file);
