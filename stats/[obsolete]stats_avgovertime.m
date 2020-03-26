% First, run the top 2 sections in stats_ROI.m
% (to load the data)



% Then, run the following code to extract a single value for each 
% subject in each condition (avg over time)

% [obsolete] For the "test bivalent first" logic (find sw$ & mix$ in bivalent,
% then avgovertime in those clusters to test for interaction)
%{
% Based on the finding - choose one of these: 
% (1) Bi_sw in LIFG, p = 0.077, between 425~450 ms
ROI_name = 'LIFG';
start_sample = 246; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 251;

% (2) Bi_sw in RIFG, p = 0.004, between 210~260 ms
ROI_name = 'RIFG';
start_sample = 203; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 213;
%}

% Now we use this to unpack main effects.
% Choose one of the following to unpack:
% (1) Main effect of Context in RACC, 125-150ms
ROI_name = 'RACC';
start_sample = 186; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 191;
% (2) Main effect of Ttype in LSMA, 230~245 ms
ROI_name = 'LSMA';
start_sample = 207; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 210;
% (3) Main effect of Ttype in LSMA, 395~425 ms
ROI_name = 'LSMA';
start_sample = 240; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 246;


data = allSubjects_ROIs.(ROI_name); % data for the current ROI
T = table(); % create an empty table (to store the output)
    
cond_names = fieldnames(data);
for j = 1:length(cond_names) % each cycle handles one cond (e.g. NatStay)
    cond_name = cond_names{j};
    allsubjects = data.(cond_name); % all subjects' time courses for this cond
    
    for i = 1:length(allsubjects) % each cycle handles one subject
        subject_timecourse = allsubjects{i}.avg;
        avg = mean(subject_timecourse(start_sample:end_sample)); % take the average over the pre-defined time interval
        T.(cond_name)(i,1) = avg; % store into the appropriate cell in the table
    end   
end

% save to csv file
start_time = data.NatStay{1,1}.time(start_sample) * 1000;
end_time   = data.NatStay{1,1}.time(end_sample) * 1000;
output_file = [ResultsFolder_ROI_thisrun 'avgovertime\\' ROI_name '_' int2str(start_time) '-' int2str(end_time) 'ms.csv'];
writetable(T, output_file);
