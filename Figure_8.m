% Figure_8.m - Dimensionality estimates
%
% Code to produce figure 8 of the 'Pitfalls of deconvolution' paper 2020 
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
methods = {'ca';'ev';'Y';'S2P_k6';'ML_e2';'LZ_k';'S2P_t6';'ML_t';'LZ_t2' };
[ncells,nt] = size(dat.timeSeriesArrayHash.value{1,2}.valueMatrix);

meths = 1:9; 
nmeths = numel(meths);

%% Dimensionality - Eigenvalues and variance explained
clear egs explained

for j = 1:nmeths 
    meth = methods{meths(j)}
    % Load data
    load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    data(find(isnan(data))) = 0;

    [V,D] = eig(cov(data'));
        
    egs(j,:) = sort(diag(D),1,'descend');
    explained(j,:) = 100*cumsum(egs(j,:))/sum(egs(j,:));
end

% For each method, work out nDims for >80 explained variance
for i = 1:nmeths
    th(i) = find(explained(i,:) >= 80,1,'first');
end

%% Plot variance explained and stripe plot side by side
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 12 12]);

for i = 1:nmeths
    subplot(2,3,[1,2,4,5])
    plot(explained(i,:),'color',cmap_ca(i,:),'linewidth',2); hold all
    ylabel('Variance explained (%)')
    xlabel('N eigenvectors')
    axis square
    
    subplot(2,3,6)
    plot(1,100*th(i)/1552,'.','color',cmap_ca(i,:),'markersize',20); hold all;
    title('Dimensionality (%)')
    ylim([0,100])
    axis square
end

subplot(2,3,[1,2,4,5])
plot([0,ncells],[80,80],'k--')
subplot(2,3,6)
set(gca, 'XTick', 1,'XTickLabel','');

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig8_dimensionality_combined'],'-r600','-dsvg')