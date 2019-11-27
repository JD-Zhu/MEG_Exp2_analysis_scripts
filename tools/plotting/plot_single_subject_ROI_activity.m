% You can use this script to plot the source activity timecourse in a single subject.
% "Source activity" can be for an ROI (i.e. virtual sensor), or for a single vertex,
% just load the desired source timecourse into the VE variable

% Step 1: 
% Load the ROI output file (e.g. A01-XC-3489_ROI.mat)
% by dragging it into the workspace

% Step 2:
% SELECT which ROI to plot (or you can loop thru all ROIs)
VE = ROI_activity.LIFG;
eventnames_real = fieldnames(VE);

% SELECT which conds to plot:
%conds = 1:9; % plot all conds
conds = [1 3 4 6]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; 
for j = conds
   plot(VE.(eventnames_real{j}).time, VE.(eventnames_real{j}).avg);
   xlim([-0.8 1]);
end
legend(eventnames_real(conds));
