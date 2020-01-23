% Figure_7.m - Pairwise correlation distributions, stability of
% correlations, differences in relative correlation structure across
% methods
%
% Code to produce figure 7 of the 'Pitfalls of deconvolution' paper 2020 
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

x = dat.timeSeriesArrayHash.value{1,2}.trial;
trials = unique(x);

meths = 1:9; 
nmeths = numel(meths);
meth_c = 1:6; % Continuous methods.
meth_s = 7:9; % Spike inference methods. 

%% Image plot of pairwise correlations
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.retina_square]);
    
CXY = {};
for j = 1:nmeths
    meth = methods{meths(j)};
    % Load data
    load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    
    pcc_data = corrcoef(data');
    pcc_data(find(eye(ncells))) = 0;
    pcc_data(find(isnan(pcc_data))) = 0;
    CXY{j} = pcc_data;
    ax(j) = subplot(3,3,j);
    imagesc(pcc_data); hold all
    title(meth_names_paper{meths(j)})
    caxis([-0.1,0.6]);
    
    % draw squares on top - shifted to make svg align
    plot([407,409.25,409.25,407,407],[409.25,409.25,407,407,409.25],'g','linewidth',widths.plot)
    plot([411,413.25,413.25,411,411],[413.25,413.25,411,411,413.25],'y','linewidth',widths.plot)
    axis square
    axis off
end

linkaxes(ax)
xlim([405,425])
ylim([405,425])

colormap cubehelix
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7c_PCC_image_all_caxis'],'-dsvg')

%% Put upper triangle of correlation matrices in a single array for easier plotting
X = ones(ncells);
Y = triu(X);
Ti = find(Y);

CXY_all = [];
for n = 1:nmeths
    meth_names_paper{meths(n)};
    T = CXY{n}(Ti);
    CXY_all(n,:) = T;
end

%% Scatter plot
siz = M-2;
alph = 0.1;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 8]);
hold all

for n = 1:nmeths
    
    npairs = numel(T);
    scatter(0.1*randn(npairs,1)+n+1*ones(npairs,1),CXY_all(n,:),siz,cmap_ca(n,:),'filled');
    plot(n+[0.75,1.25;0.75,1.25;0.75,1.25]',[prctile(CXY_all(n,:),[5,50,95]);prctile(CXY_all(n,:),[5,50,95])],'k','linewidth',widths.plot)
    
end

ylim([-0.4,1]);
xlim([1.25,10.75])
plot([0,11],[0,0],'k--','linewidth',widths.plot)
set(gca,'XTick',2:n+1,'XTickLabel',meth_names_paper(meths),'XTickLabelrotation',45)
pbaspect([2 1 1])
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7a_PCC_distribution_scatter_tinydots'],'-r600','-dsvg','-painters')

%% Data subset correlation analysis (uncomment to run)
% for d = 1:nmeths
%    pcc_dist = zeros(100,100); 
%    data = eval(data_ID{meths(d)});
%    pcc_all = corrcoef(data');
%    pcc_all(find(eye(ncells))) = 0;
%    pcc_all(find(isnan(pcc_all))) = 0;
% 
%    parfor i = 1:100
%        for j = 1:100
%            all_times = randperm(nt);
%            sample_times = all_times(1:round(nt*0.01*i));
%            data_subset = data(:,sample_times);
%            pcc_subset = corrcoef(data_subset');
%            pcc_subset(find(isnan(pcc_subset))) = 0;
%            pcc_subset(find(eye(ncells))) = 0;
%            pcc_compare = corrcoef(pcc_all,pcc_subset);
%            pcc_dist(i,j) = pcc_compare(1,2);
%        end
%    end
%    cxy_pcc{d} = pcc_dist;
% end
% save cxy_pcc_19 cxy_pcc

%% Load and plot
load Data/cxy_pcc_19.mat
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 6]);
plot([50,50],[0,1],'k--','linewidth',widths.plot); hold all

for i = 1:nmeths
    dotc = (cmap_ca(i,:)+0.2); dotc(find(dotc>1)) = 1;
    myeb(1:100,mean(cxy_pcc{i}'),std(cxy_pcc{i}'),cmap_ca(i,:),dotc)
end

pbaspect([2 1 1])
xlim([0,100])
xlabel('Subset duration (% of total)')
ylabel('PCC(Cxy_{subset},Cxy_{total})')
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7b_cxy_stability_mean_wide'],'-r600','-dsvg','-painters')

%% CXY of CXYs, all methods, then separate for spike inference and continuous methods
% Using Spearman's rank correlation
sxycxy = zeros(nmeths,nmeths);
for j = 1:nmeths
    for k = 1:nmeths
        if j == k
        else
            sxy = corr(CXY_all(j,:)',CXY_all(k,:)','Type','Spearman');
            sxycxy(j,k) = sxy;
        end
    end
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 8]);
imagesc(sxycxy)
axis square
set(gca,'XTick',1:nmeths,'XTickLabel',meth_names_paper(meths),'XTickLabelrotation',45)
set(gca,'YTick',1:nmeths,'YTickLabel',meth_names_paper(meths),'YTickLabelrotation',45)
colorbar; colormap cubehelix;

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7d_sxycxy_all_rot45'],'-r600','-dsvg','-painters')

%% Continuous methods only
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 6]);
imagesc(sxycxy(meth_c,meth_c))
axis square
set(gca,'XTick',1:numel(meth_c),'XTickLabel',meth_names_paper(meths(meth_c)),'XTickLabelrotation',45)
set(gca,'YTick',1:numel(meth_c),'YTickLabel',meth_names_paper(meths(meth_c)),'YTickLabelrotation',45)
colorbar; colormap cubehelix;

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7e_sxycxy_continuous_rot'],'-r600','-dsvg')

%% Spike inference methods only
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 6 6]);
imagesc(sxycxy(meth_s,meth_s))
axis square
set(gca,'XTick',1:numel(meth_s),'XTickLabel',meth_names_paper(meths(meth_s)),'XTickLabelrotation',45)
set(gca,'YTick',1:numel(meth_s),'YTickLabel',meth_names_paper(meths(meth_s)),'YTickLabelrotation',45)
colorbar; colormap cubehelix;
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig7e_sxycxy_spikeinf_rot'],'-r600','-dsvg')