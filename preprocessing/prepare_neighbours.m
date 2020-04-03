% prepare neighbours - only need to do this once & save for later

%{
Default distance is 4.
For MQ system:
  default distance gives on average 4.6 neighbours per channel;       (prob too few)
  distance of 5 gives on average 8.0 neighbours per channel;
  distance of 6 gives on average 11.4 neighbours per channel;
  distance of 7 gives on average 15.0 neighbours per channel;         
  'triangulation' method gives on average 7.6 neighbours per channel. (I used this for MEG Exp 1)
%}

cfg        = [];
cfg.method = 'distance'; % or 'triangulation'
cfg.neighbourdist = 7; % distance threshold (only applicable to 'distance' method)
%cfg.feedback = 'yes'; % this will produce a plot showing the neighbour relationships
% http://www.fieldtriptoolbox.org/faq/how_can_i_define_neighbouring_sensors/

neighbours = ft_prepare_neighbours(cfg, all_blocks); % pass in any 'S1_proprocessed_data.mat'
save('neighbours.mat', 'neighbours');