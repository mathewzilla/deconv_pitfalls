% Figure_1.m - Ground truth data example, firing rate error, correlation
% coefficient/Error Rate comparison over a parameter sweep
%
% Code to produce figure 1 of the 'Pitfalls of deconvolution' paper 2020 
% biorxiv: https://www.biorxiv.org/content/10.1101/871137v1
%
% M. Evans January 2020

%% Set figure parameters
clear all
run figure_properties_deconv.m

% Load data
load Data/CAI2.mat

%% Ground truth example data
load('Data/data_20120417_cell4_001.mat')
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.wide]);

plot(obj.timeSeriesArrayHash.value{4}.time,obj.timeSeriesArrayHash.value{4}.valueMatrix,'linewidth',widths.plot,'color',[.5,.5,.5])
hold all
plot(obj.timeSeriesArrayHash.value{5}.time(find(obj.timeSeriesArrayHash.value{5}.valueMatrix)),-ones(numel(find(obj.timeSeriesArrayHash.value{5}.valueMatrix)),1),'k*','linewidth',0.5)
plot(obj.timeSeriesArrayHash.value{1}.time,2.5+zscore(obj.timeSeriesArrayHash.value{1}.valueMatrix),'linewidth',widths.plot,'color','k')

xlim([168,218])
ylim([-2,6])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'Fig1a_ground_truth_example'],'-r600','-dpdf'); 

%% Normalized firing rate error (60Hz and 7Hz)
load Data/FR_data_210318

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.dists]);


subplot(2,1,1);
plot([0.25,1.75],[0,0],'k--','linewidth',widths.plot); hold all;
subplot(2,1,2);
plot([0.25,1.75],[0,0],'k--','linewidth',widths.plot); hold all; 

for i = 1:3
    % Normalized firing rate error
    subplot(2,1,1);
    fr_est = FR_data(:,1,i)';
    plot((0.5*i)+[-0.125,-0.025],mean((fr_est-fr)./fr)*ones(2,1),'k','linewidth',widths.plot)
    plot(0.5*i*ones(21,1) - 0.1 + 0.05*rand(21,1),(fr_est-fr)./fr,'.','color',cmap_gt(7,:),'markersize',5);
    
    fr_est = FR_data(:,2,i)';
    plot((0.5*i)+[0.025,0.125],mean((fr_est-fr)./fr)*ones(2,1),'k','linewidth',widths.plot)
    plot(0.5*i*ones(21,1) + 0.05 + 0.05*rand(21,1),(fr_est-fr)./fr,'.','color',cmap_gt(2,:),'markersize',5);
    
    subplot(2,1,2)
    hold all
    fr_est = FR_data(:,3,i)';
    plot((0.5*i)+[-0.125,-0.025],mean((fr_est-fr)./fr)*ones(2,1),'k','linewidth',widths.plot)
    plot(0.5*i*ones(21,1) - 0.1 + 0.05*rand(21,1),(fr_est-fr)./fr,'.','color',cmap_gt(7,:),'markersize',5);
    fr_est = FR_data(:,4,i)';
    plot((0.5*i)+[0.025,0.125],mean((fr_est-fr)./fr)*ones(2,1),'k','linewidth',widths.plot)
    plot(0.5*i*ones(21,1) + 0.05 + 0.05*rand(21,1),(fr_est-fr)./fr,'.','color',cmap_gt(2,:),'markersize',5);

end

subplot(2,1,1);
ylim([-1,2])
xlim([0.25,1.75])
title('60 Hz')
set(gca,'XTick',[0.5,1,1.5],'XTickLabel',{'Suite2P','MLSpike','LZero'}');
ylabel('Normalized firing rate error')

subplot(2,1,2);
ylim([-1,2])
xlim([0.25,1.75])
title('7 Hz')
set(gca,'XTick',[0.5,1,1.5],'XTickLabel',{'Suite2P','MLSpike','LZero'}');

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'Fig1bc_gt_stripe_plot'],'-dpdf'); 


%% MLSpike parameter sweep figure (PCC/ER vs firing rate)
load Data/MLSpike.mat
load Data/MLSpike_FRs.mat
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 6]);

par_names = {'a';'tau';'sigma'}; par_lets = {'A';'B';'C'}; meths = {'PCC';'ER'};
for i = 1:3
    for j = 1:2
        for k = 1:21
        ax(j*3+i-3) = subplot(2,3,j*3+i-3); 
        plot(eval([par_lets{i},'_FR(k,:)']),eval([par_lets{i},'_',meths{j},'(k,:)']),'color',cmap_gt(k,:),'linewidth',widths.plot);
        hold all
        if i == 2; xlabel('Estimated Firing Rate'); end 
        if i == 1; ylabel(meths(j)); end
        title(par_names{i})
        axis square
        ylim([0,1])
        end
    end
end

linkaxes(ax,'x')
xlim([0,5])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath 'Fig1de_MLSpike_FRvsER'],'-dpdf');

%% Suite 2P parameter sweep figure
load Data/CAI2.mat
load Data/S2P_results_130218.mat

num_spks_ds = S2P_results{1};
results_ds = S2P_results{2};
ER_ds = S2P_results{3};
PCC_ds = S2P_results{4};
num_spks_rt = S2P_results{5};
results_ds = S2P_results{6};
ER_rt = S2P_results{7};
PCC_rt = S2P_results{8};

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.dists]);

for k = 1:21
    subplot(1,2,1)
    plot(num_spks_rt(k,:)'/(numel(CAI{k}.f)/60),PCC_rt(k,:)','color',cmap_gt(k,:),'linewidth',widths.plot)
    axis square; hold all
    xlabel('Inferred event rate')
    ylabel('PCC')
    subplot(1,2,2)
    plot(num_spks_rt(k,:)'/(numel(CAI{k}.f)/60),ER_rt(k,:)','color',cmap_gt(k,:),'linewidth',widths.plot)
    axis square; hold all
    xlabel('Inferred event rate')
    ylabel('ER')
end

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'Fig1f_S2P_FRvsER_Hz'],'-dpdf');