addpath(genpath('/home/phornauer/Git/Controllability'))
addpath('/home/phornauer/Git/npy-matlab/npy-matlab')

%% Load data
fr_data_path = '/home/phornauer/Git/Controllability/Data/smoothed.npy';
fr_data = readNPY(fr_data_path);

%% DDC params
thres = 0;
TR = 1;
win = 100;

%% Cut data
t_min = 3.5; %[s]
t_min = t_min*1000;
t_max = 5.5;%[s]
t_max = t_max*1000; %conversion to ms
fr_input = fr_data(:,t_min+1:t_max);

id_selection = [1:272];
fr_input = fr_input(id_selection,:);
% figure;plot(sum(fr_input,1))
%%
act_idx = find(sum(fr_input,2)<10000);
fr_input(act_idx,:) = []; 
%% Calculate DDC across time windows
win_idx = 1:win:length(fr_input);
input_array = cell(1,length(win_idx));
for w = 1:length(win_idx)
    input_array{w} = fr_input(:,win_idx(w):win_idx(w)+win-1);
end
DDC = cell(1,length(win_idx));

parfor w = 1:length(input_array)
    DDC{w} = calc_DDC(input_array{w}',thres,TR);
end
%%
figure;
tiledlayout('flow')
for i = 1:length(DDC)
    nexttile
    imagesc(DDC{i});colorbar
    caxis([-2e20 2e20])
end
%%
corr_mat = zeros(length(DDC));
for i = 1:length(DDC)-1
    for j = i+1:length(DDC)
        cc = corrcoef(DDC{i},DDC{j});
        corr_mat(i,j) = cc(1,2);
        corr_mat(j,i) = cc(1,2);
    end
end
%%
figure;
subplot(2,1,1)
imagesc(corr_mat);colorbar
subplot(2,1,2)
plot(sum(fr_input,1))
%%
[M,Q] = community_louvain(corr_mat,0.5,[],'negative_asym');
[X,Y,indsort] = grid_communities(M);
figure;imagesc(corr_mat(indsort,indsort))