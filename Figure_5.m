% Figure_5.m - Tuned cells, shuffle test examples, tuned cell agreement,
% image plots of tuning agreement
%
% Code to produce figure 5 of the 'Pitfalls of deconvolution' paper 2020 
% biorxiv: https://www.biorxiv.org/content/10.1101/871137v1
%
% M. Evans January 2020

%% Set figure parameters
clear all
run figure_properties_deconv.m

% Load session/task data
load Data/an197522_2013_03_07.mat

%% file names and other housekeeping
data_ID = {'ca';'ev';'Y_e';'S2P_k6';'ML_e2';'LZ_k';'S2P_t6';'ML_t';'LZ_t2'};

meth_names_paper = {'Calcium';'Peron';'Yaksi';'Suite2P_{kernel}';'MLSpike_{kernel}';'LZero_{kernel}';'Suite2P_{events}';'MLSpike_{events}';'LZero_{events}'};

[ncells,nt] = size(dat.timeSeriesArrayHash.value{1,2}.valueMatrix);

x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

meths = 1:9; 
nmeths = numel(meths);
meth_c = 1:6; % Continuous methods.
meth_s = 7:9; % Spike inference methods. 

%% Individual cell shuffle examples
ca = dat.timeSeriesArrayHash.value{1,2}.valueMatrix;

% Replace NaNs with zeros
ca(isnan(ca)) = 0;

% data psth
psth_ca = nan(100,50);
for c = 1:100
    s_ca = [];
    for i = 1:numel(trials)
        s_ca(i,:) = [ca(c,find(x == trials(i))),nan(1,100-numel(find(x == trials(i))))];
    end
    psth_ca(c,:) = nanmean(s_ca(:,1:50),1);
end

% Load example shuffle tests for two cells
load Data/shuff_ca_8_9_10000.mat;

% cell  9 - tuned
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.dists]);
nfth = []; % 95th percentile
subplot(1,2,1)
plot(squeeze(shuff_ca(2,:,:)),'color',[.5,.5,.5]);
[mx,ix] = max(squeeze(shuff_ca(2,:,:)));
nfth(1) = prctile(mx,95);
hold all
plot([0,50],nfth(1)*ones(2,1),'r','linewidth',widths.plot)
plot(psth_ca(9,:),'k','linewidth',widths.plot)
axis square
xlabel('Time (frames)')
ylabel('dF/F')

% cell 8 - not tuned
subplot(1,2,2)
plot(squeeze(shuff_ca(1,:,:)),'color',[.5,.5,.5]);
[mx,ix] = max(squeeze(shuff_ca(1,:,:)));
nfth(2) = prctile(mx,95);
hold all
plot([0,50],nfth(2)*ones(2,1),'r','linewidth',widths.plot)
plot(psth_ca(8,:),'k','linewidth',widths.plot)
axis square

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig5a_10K_tuning_shufftest_example_cells_9_8'],'-r600','-dsvg','-painters')

%% Tuned cells
load('Data/tuning_shuffle_results.mat')

num_tg = num_tuned(meths);
Tuned_g = Tuned{meths};
clear CIs JIs

% Of only the 'best' methods
JIs(1,:) = betainv(0.025,num_tg+0.5, 1552 - num_tg +0.5);
JIs(2,:) = betainv(0.975,num_tg+0.5, 1552 - num_tg +0.5);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]);

norm_tuned = num_tg/1552;
bar(norm_tuned,'k')
set(gca, 'XTick', [1:nmeths]);
set(gca, 'XTickLabels', meth_names_paper(meths),'XTickLabelrotation',45); 
hold all
errorbar(1:numel(norm_tuned),norm_tuned,norm_tuned-JIs(1,:),JIs(2,:)-norm_tuned,'color',[.5,.5,.5],'linewidth',widths.plot,'Linestyle','none')
ylim([0,1]);
xlim([0.5,nmeths+.5])
ylabel('P tuned cells')
axis square
title('Proportion of tuned cells')

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig5b_P_tuned_all'],'-dsvg')

%% Degree of agreement between methods on tuned status of a given cell
% Which cells are 'tuned' ? Are the the same across methods?
tuned_array = zeros(nmeths,ncells);
for j = 1:nmeths
    tuned_cells = counted_tuned{meths(j)};
    tuned_array(j,tuned_cells) = 1;
end

% Degree of agreement between methods
agreed = sum(tuned_array);
[agreement,sort_agreed] = sort(agreed);

% Histogram of agreement stats
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]);

h = hist(agreement,nmeths+1);
bar(linspace(0,nmeths,nmeths+1),h,'k')
title('Degree of agreement between methods')
xlabel('N methods that agree a given cell is tuned')
ylabel('N cells')
axis square

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig5c_tuned_cell_all_agreement'],'-dsvg')

%% Separate kernel vs event methods
% continuous
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]);

% Degree of agreement between continuous methods
agreed = sum(tuned_array(meth_c,:));
[agreement,sort_agreed] = sort(agreed);

h = hist(agreement,7);
bar(linspace(0,6,7),h,'k')
xlim([-1,7])
ylim([0,ncells])
axis square

title('Continuous methods')
xlabel('N methods that agree a given cell is tuned')
ylabel('N cells')

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig5d_tuned_cells_agreement_continuous'],'-dsvg')

% spike inference
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]);

% Degree of agreement between spike inference methods
agreed = sum(tuned_array(meth_s,:));
[agreement,sort_agreed] = sort(agreed);

h = hist(agreement,4);
bar(linspace(0,3,4),h,'k')
xlim([-1,4])
ylim([0,ncells])
axis square

title('Spike inference methods')
xlabel('N methods that agree a given cell is tuned')
ylabel('N cells')

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig5d_tuned_cells_agreement_spike'],'-dsvg')

%% Image plot of tuned cells sorted by peak time, with peak overlaid
% (Doing this separately for ~50 cells classed as tuned by 1, 3, 6 and 9 methods, for all vs separate methods)
load(['Data/data_peaks_and_psth.mat'])

% Best methods
agreed = sum(tuned_array);
[agreement,sort_agreed] = sort(agreed);
Nagreed_list = [1,3,6,9];
for i = 1:4
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.retina_square]);
    
    Nagreed = Nagreed_list(i);
    to_plot = find(agreed == Nagreed,50);
    ca_psth = data_psth{1};
    ca_peaks = data_peaks{1};
    ca_ix = data_ix{1};
    [sorted_peak,sorted_ix] = sort(ca_ix(to_plot,1));
    
    subplot(3,3,1);
    % imagesc(ca_psth(one(sorted_ix),:))
    imagesc(zscore(ca_psth(to_plot(sorted_ix),:)')')
    hold all
    % plot(sorted_peak,1:numel(one),'w.')
    title(meth_names_paper{1})
    axis square
    
    for j = 2:nmeths
        
        d_psth = data_psth{meths(j)};
        d_peaks = data_peaks{meths(j)};
        d_ix = data_ix{meths(j)};
        
        subplot(3,3,j);
        imagesc(d_psth(to_plot(sorted_ix),:))
        imagesc(zscore(d_psth(to_plot(sorted_ix),:)')')
        hold all
        %     plot(d_ix(one(sorted_ix)),1:numel(one),'w.')
        title(meth_names_paper{meths(j)})
        axis square
    end
    suptitle(['N agreed = ',num2str(Nagreed)])
    colormap cubehelix
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
    print([exportpath,'Fig5_agreement_image_',num2str(Nagreed),'_fontsize_8'],'-r600','-dpdf')
end