% Figure_6.m - Touch tuning. One cell example and a histogram
%
% Code to produce figure 6 of the 'Pitfalls of deconvolution' paper 2020 
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

%% Find the first touch on every trial
first_touch = [];
for i = 1:numel(trials)
    tt = [];
    x1 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{1} == i,1,'first');
    tt1 = dat.eventSeriesArrayHash.value{2}.eventTimes{1}(x1);
    x2 = find(dat.eventSeriesArrayHash.value{2}.eventTrials{2} == i,1,'first');
    tt2 = dat.eventSeriesArrayHash.value{2}.eventTimes{2}(x2);
    tt = [tt1,tt2];
    if find(tt)
        first_touch = [first_touch, min(tt(find(tt)))];
    end
    
end

% Closest frame to touch time (as touch is measured in milliseconds)
ftf = [];
for i = 1:numel(first_touch)
    [mn,im] = min(abs(dat.timeSeriesArrayHash.value{2}.time - first_touch(i)));
    ftf(i) = im; 
end

%% Actual touch triggered distribution
TTA = {};
for j = 1:nmeths
    load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    tth = zeros(ncells,15);

    for c = 1:ncells
        this_c = data(c,:);
        tt_sample = zeros(numel(ftf),15);
        for i = 1:numel(ftf)
            tt_sample(i,:) = this_c(ftf(i)-7 : ftf(i)+7);
        end
        tth(c,:) = nanmean(tt_sample,1);
    end
    TTA{j} = tth;
end

%% Plot cell 495 as it is clearly touch tuned
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 5]);

for i = 1:nmeths 
    plot(zscore(TTA{i}(495,:)),'color',cmap_ca(i,:),'linewidth',widths.plot); hold all;
end
xlabel('Time (frames)')
ylabel('dF/F (z-scored)')
title('Touch triggered average (Cell 495)')
plot([7.5,7.5],[-2,4],'k--')
axis square
legend({meth_names_paper{meths},'Touch'},'location','bestoutside')

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig6a_touch_triggered_average_cell_495'],'-dsvg')

%% Shuffle test for touch tuning - uncomment to run. It takes a while!
% % 10K repeats of matched sample distribution vs touch triggered amplitude peak distribution
% % then take empirical p value
% 
% TT_shuff = {};
% 
% for j = 1:nmeths
%     load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
%     data = eval(data_ID{meths(j)});
% 
%     N_tt = zeros(ncells,1);    % Number of shuffles greater than empirical max
%     P_tt = zeros(ncells,1);    % empirical p value of touch tuning
%     
%     for c = 1 :ncells
%         display(['Meth ',num2str(j),'. Cell ',num2str(c)])
%         this_c = data(c,:);
%         tt_sample = zeros(numel(ftf),15);
% 
%         for i = 1:numel(ftf)
%             tt_sample(i,:) = this_c(ftf(i)-7 : ftf(i)+7);
%         end
%         
%         tth = nanmean(tt_sample(:,1:15),1);
%         tt_max = max(tth);
%         
%         % Do the same calculation 10K times on random chunks of data
%         
%         parfor n = 1:10000
%             % Equivalently sized array of random segments
%             tt_shuff = zeros(numel(ftf),15);
%             rsamp = randsample(1:nt-15,numel(ftf));
%             
%             % Random sample
%             for i = 1:numel(ftf)
%                 tt_shuff(i,:) = this_c(rsamp(i):rsamp(i)+14);
%             end
%             
%             shuff_max(n) = max(nanmean(tt_shuff(:,1:15),1));
%             
%         end
%         
%         % Count number of shuffled peaks are higher than empirical peak
%         N_tt(c) = numel(find(shuff_max>tt_max));
%         
%         % Empirical p value of touch tuning
%         P_tt(c) = N_tt(c)./10000;
%     end
%     TT_shuff{j,1} = N_tt;
%     TT_shuff{j,2} = P_tt;
% 
% end 
%
% save Data/TT_shuff.mat TT_shuff

%% Recompute touch tuning based on Benjamimi-Hochberg correction of p values
% Load shuffle results
load('Data/TT_shuff.mat')
for j = 1:nmeths
    p_values = TT_shuff{j,1}./10000; 
    
    [H,T] = benjaminihochberg(p_values',0.05); % Mark Humprhies' code
    TT_shuff{j,3} = H;
end

% 95% confidence interval (Jeffreys Interval) on number of tuned cells
for j = 1:nmeths
    N_tuned(j) = numel(find(TT_shuff{j,3}));
end

JIs(1,:) = betainv(0.025,N_tuned+0.5, 1552 - N_tuned +0.5);
JIs(2,:) = betainv(0.975,N_tuned+0.5, 1552 - N_tuned +0.5);

%% Bar graph of N touch tuned
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); % retina_square for OG fig
bar(N_tuned,'k')
set(gca, 'XTickLabels', meth_names_paper(meths),'XTickLabelRotation',45)
ylabel('N touch tuned')
hold all
errorbar(1:numel(N_tuned),N_tuned,N_tuned-ncells*JIs(1,:),ncells*JIs(2,:)-N_tuned,'color',[.5,.5,.5],'linewidth',widths.plot,'Linestyle','none')

axis square

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

print([exportpath,'Fig6b_N_touch_tuned'],'-dsvg')



