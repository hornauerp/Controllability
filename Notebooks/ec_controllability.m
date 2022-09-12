%% Import external packages
cd('/home/phornauer/Git/Scheid_EC_network_tools')
addpath('/home/phornauer/Git/npy-matlab/npy-matlab')
addpath(genpath('/home/phornauer/Git/Scheid_EC_network_tools'))
addpath(genpath('/home/phornauer/Git/GenLouvain'))
addpath('/home/phornauer/Git/NetworkCommunityToolbox')
%% Load data
fr_data_path = '/home/phornauer/Git/Controllability/Data/smoothed.npy';
fr_data = readNPY(fr_data_path);

%%
spk_t = double(readNPY('/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/manuelsc/Christian/Chip2602/baseline/spike_times.npy'));
spk_c = readNPY('/net/bs-filesvr02/export/group/hierlemann/intermediate_data/Mea1k/manuelsc/Christian/Chip2602/baseline/spike_templates.npy');
%% Cut data
t_min = 3.5; %[s]
t_min = t_min*1000;
t_max = 5.5;%[s]
t_max = t_max*1000; %conversion to ms
fr_input = fr_data(:,t_min+1:t_max);

id_selection = [1:50];
fr_input = fr_input(id_selection,:);
figure;plot(sum(fr_input,1))
%%
windowSize = 100; % [ms]
Networks = getNets('./', fr_input, windowSize, 0.5, [0, 0.1, .5]);
%%
figure;
subplot(2,1,1)
imagesc(Networks.wSim_0)
subplot(2,1,2)
plot(sum(fr_input,1))
%% Partition Network similarity matrix into desired number of communities
netw = Networks.wSim_0;
[~,nan_idx] = max(sum(isnan(netw)));
netw(nan_idx,:) = []; netw(:,nan_idx) = [];
init_gamma = (0.8:0.05:1.2); % Set of initial spatial resolution parameters 
nTarget = 3;               % Target number of communities 
Qiter = 100;               % iterations of Louvain algorithm
Partitions = findSimComms(netw, init_gamma, Qiter, nTarget);

%%
gamma = Partitions.gamma(Partitions.nTargInd);
k = full(sum(netw));
B = full(netw - gamma*(k'*k)/sum(k));

[S,Q]=genlouvain(B, 10000, 0);
%%
figure;
subplot(3,1,1)
plot(sum(fr_input,1))
subplot(3,1,2)
imagesc(B)
hold on
xline(find(diff(M))+0.5)
subplot(3,1,3)
plot(M)
%%
figure;imagesc(Networks.pcm(:,:,1));colorbar
%%
corrcoef(Networks.pcm(:,:,6),Networks.pcm(:,:,8))
%%
t = 200; %time in [s]
spk_max = t*20000;
sel_times = spk_t(spk_t<spk_max)/20000;
sel_channel = spk_c(spk_t<spk_max);
binned = histcounts(sel_times,0:0.01:double(max(sel_times)));
figure;
subplot(2,1,1)
plot(sel_times,sel_channel,'k.','MarkerSize',0.1)

subplot(2,1,2)
plot(binned)