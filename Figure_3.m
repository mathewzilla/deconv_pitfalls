% Figure_3.m - Example data from df/f calcium and 8 deconvolution methods
%
% Code to produce figure 3 of the 'Pitfalls of deconvolution' paper 2020 
% biorxiv: https://www.biorxiv.org/content/10.1101/871137v1
%
% M. Evans January 2020

%% Set figure parameters
clear all
run figure_properties_deconv.m

% Load session/task data
load Data/an197522_2013_03_07.mat

%% file names and other task data
data_ID = {'ca';'ev';'Y_e';'S2P_k6';'ML_e2';'LZ_k';'S2P_t6';'ML_t';'LZ_t2'};
meth_names_paper = {'Calcium';'Peron';'Yaksi';'Suite2P_{kernel}';'MLSpike_{kernel}';'LZero_{kernel}';'Suite2P_{events}';'MLSpike_{events}';'LZero_{events}'};

% Extract trial and session info
x = dat.timeSeriesArrayHash.value{1,2}.trial;

meths = 1:9; 
nmeths = numel(meths);

%% example data figure
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.dists]);

t_array = dat.timeSeriesArrayHash.value{2}.time/1000; 
start_t = 5200; 
end_t = 6500;
for j = 1:nmeths
    load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    norm_data = data(27,:)/max(data(27,:));
    plot(t_array-t_array(start_t),-j+0.9*norm_data,'color',cmap_ca(j,:),'linewidth',widths.plot);
    hold all
    
end
set(gca,'ytick',linspace(-9,-1,9),'yticklabel',meth_names_paper(fliplr(meths)))
xlabel('Time (s)')
xlim([0,t_array(end_t) - t_array(start_t)])

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig3b_example_data'],'-r600','-dpdf');
