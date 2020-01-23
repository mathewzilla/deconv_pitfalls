% Figure_2.m - Variation in best-fit spike deconvolution parameters across
% ground-truth recordings, delta Error Rate with change in parameter values
%
% Code to produce figure 2 of the 'Pitfalls of deconvolution' paper 2020 
% biorxiv: https://www.biorxiv.org/content/10.1101/871137v1
%
% M. Evans January 2020

%% Set figure parameters
clear all
run figure_properties_deconv.m

% Load data
load Data/CAI2.mat

%% Plot the best Suite2P parameters
load Data/S2P_metrics_210318.mat
Thresholds = logspace(-1,2,13);
M = 10;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 12 8]); 

% switch for grey lines between PCC and ER results (not in paper)
greylines = 0;

% switch for connecting same cells across recordings
connect_cells = 1;

% Which recordings correspond to which cells?
% Plotting lines connecting the same cells
% 9 cells, 21 recordings
cell_N = [1,1,1,2,2,3,3,4,4,4,4,5,6,6,6,7,7,7,8,9,9];

%% Figure code

% Density plot
plt(1) = subplot(6,6,1);
[k,xi] = ksdensity(best_params_ER,'bandwidth',0.5);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC,'bandwidth',0.5);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,14]); pbaspect([2,1,1]);
set(gca,'xtick',[1,5,13],'xticklabel',[10^-1,10^0,10^2])

plt(2) = subplot(6,6,[7,13]);
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(i),best_params_PCC(i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all
    plot(best_params_PCC(i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    
end
axis square
xlim([1,14]);
set(gca,'xtick',[1,5,13],'xticklabel',[10^-1,10^0,10^2])

% Downsampled version
load Data/S2P_metrics_ds_210318.mat

% Density plot
plt(3) = subplot(6,6,19);
[k,xi] = ksdensity(best_params_ER,'bandwidth',0.5);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC,'bandwidth',0.5);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,14]); pbaspect([2,1,1]);
set(gca,'xtick',[1,5,13],'xticklabel',[10^-1,10^0,10^2])

plt(4) = subplot(6,6,[25,31]);
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end


for i = 1:21
    if greylines
        plot([best_params_ER(i),best_params_PCC(i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all
    plot(best_params_PCC(i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
end
axis square
xlim([1,14]);
set(gca,'xtick',[1,5,13],'xticklabel',[10^-1,10^0,10^2])

% Plot the best MLSpike parameters
load Data/MLSpike_best_params
a_list = logspace(-2,0,21); tau_list = logspace(-2,0.6990,21); sigma_list = logspace(-2,0,21);

% Density plots
plt(5) = subplot(6,6,2);
[k,xi] = ksdensity(best_params_ER(1,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(1,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,10^-1,10^0])
xlim([1,21]);pbaspect([2,1,1]);

plt(6) = subplot(6,6,3);
[k,xi] = ksdensity(best_params_ER(2,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(2,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,18,21],'xticklabel',[10^-2,2,5])
xlim([1,21]);pbaspect([2,1,1]);

plt(7) = subplot(6,6,4);
[k,xi] = ksdensity(best_params_ER(3,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(3,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,0.1,1])
xlim([1,21]);pbaspect([2,1,1]);


plt(8) = subplot(6,6,[8,14]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(1,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(1,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(1,i),best_params_PCC(1,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(1,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(1,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]);
    set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,10^-1,10^0])
end
axis square

plt(9) = subplot(6,6,[9,15]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(2,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(2,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(2,i),best_params_PCC(2,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(2,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(2,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]); 
    set(gca,'xtick',[1,18,21],'xticklabel',[10^-2,2,5])
end
axis square

plt(10) = subplot(6,6,[10,16]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(3,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(3,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(3,i),best_params_PCC(3,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(3,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all; 
    plot(best_params_PCC(3,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]); 
    set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,0.1,1])
end
axis square

% Downsampled version
load Data/MLSpike_best_params_ds


% Density plots
plt(11) = subplot(6,6,20);
[k,xi] = ksdensity(best_params_ER(1,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(1,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,10^-1,10^0])
xlim([1,21]); pbaspect([2,1,1]);

plt(12) = subplot(6,6,21);
[k,xi] = ksdensity(best_params_ER(2,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(2,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,18,21],'xticklabel',[10^-2,2,5])
xlim([1,21]); pbaspect([2,1,1]);

plt(13) = subplot(6,6,22);
[k,xi] = ksdensity(best_params_ER(3,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(3,:),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,0.1,1])
xlim([1,21]); pbaspect([2,1,1]);



plt(14) = subplot(6,6,[26,32]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(1,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(1,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(1,i),best_params_PCC(1,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(1,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(1,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]); 
    set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,10^-1,10^0])
end
axis square

plt(15) = subplot(6,6,[27,33]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(2,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(2,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(2,i),best_params_PCC(2,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(2,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all; 
    plot(best_params_PCC(2,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]);
    set(gca,'xtick',[1,18,21],'xticklabel',[10^-2,2,5])
end
axis square

plt(16) = subplot(6,6,[28,34]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(3,this_c),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(3,this_c),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(3,i),best_params_PCC(3,i)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(3,i),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all; 
    plot(best_params_PCC(3,i),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21]);
    set(gca,'xtick',[1,11,21],'xticklabel',[10^-2,0.1,1])
end
axis square

% LZero
load Data/LZero_metrics_210318.mat
lambdas = [logspace(-1,1,13),linspace(11,20,10)];
scales = [logspace(-1,1,13),linspace(11,20,10)];

% Density plots
plt(17) = subplot(6,6,5);
[k,xi] = ksdensity(best_params_ER(:,1),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(:,1),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,21]); pbaspect([2,1,1]);
set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])

plt(18) = subplot(6,6,6);
[k,xi] = ksdensity(best_params_ER(:,2),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(:,2),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,21]); pbaspect([2,1,1]);
set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])


plt(19) = subplot(6,6,[11,17]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c,1),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c,1),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(i,1),best_params_PCC(i,1)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i,1),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;

    plot(best_params_PCC(i,1),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21])
    set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])
end
axis square

plt(20) = subplot(6,6,[12,18]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c,2),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c,2),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(i,2),best_params_PCC(i,2)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i,2),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(i,2),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21])
    set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])
end
axis square

% Downsampled version
load Data/LZero_metrics_ds_210318.mat

% Density plots
plt(21) = subplot(6,6,23);
[k,xi] = ksdensity(best_params_ER(:,1),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(:,1),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,21]); pbaspect([2,1,1]);
set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])

plt(22) = subplot(6,6,24);
[k,xi] = ksdensity(best_params_ER(:,2),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(2,:)); hold all
[k,xi] = ksdensity(best_params_PCC(:,2),linspace(1,21,100),'bandwidth',1);
plot(xi,k,'linewidth',widths.plot,'color',cmap_gt(7,:))
xlim([1,21]); pbaspect([2,1,1]);
set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])


plt(23) = subplot(6,6,[29,35]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c,1),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c,1),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(i,1),best_params_PCC(i,1)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i,1),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(i,1),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21])
    set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])
end
axis square

plt(24) = subplot(6,6,[30,36]); axis square
if connect_cells
    for i = 1:9
        this_c = find(cell_N == i);
        plot(best_params_ER(this_c,2),this_c,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
        plot(best_params_PCC(this_c,2),this_c+0.25,'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
end
for i = 1:21
    if greylines
        plot([best_params_ER(i,2),best_params_PCC(i,2)],[i+0.125,i+0.125],'linewidth',widths.plot,'color',[.5,.5,.5]); hold all
    end
    plot(best_params_ER(i,2),i,'.','color',cmap_gt(2,:),'markersize',M)
    hold all;
    plot(best_params_PCC(i,2),i+0.25,'.','color',cmap_gt(7,:),'markersize',M)
    xlim([1,21])
    set(gca,'xtick',[1,7,20],'xticklabel',[10^-1,1,20])
end
axis square

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

%% align sub figures 
pos1 = get(plt(1), 'Position');
pos2 = get(plt(2), 'Position');

pos3 = get(plt(3), 'Position');
pos4 = get(plt(4), 'Position');

pos5 = get(plt(5), 'Position'); pos6 = get(plt(6), 'Position'); pos7 = get(plt(7), 'Position');
pos8 = get(plt(8), 'Position'); pos9 = get(plt(9), 'Position'); pos10 = get(plt(10), 'Position');

pos11 = get(plt(11), 'Position'); pos12 = get(plt(12), 'Position'); pos13 = get(plt(13), 'Position');
pos14 = get(plt(14), 'Position'); pos15 = get(plt(15), 'Position'); pos16 = get(plt(16), 'Position');

pos17 = get(plt(17), 'Position'); pos18 = get(plt(18), 'Position');
pos19 = get(plt(19), 'Position'); pos20 = get(plt(20), 'Position');

pos21 = get(plt(21), 'Position'); pos22 = get(plt(22), 'Position');
pos23 = get(plt(23), 'Position'); pos24 = get(plt(24), 'Position');

% % Set and x alignment and width of first axes equal to second
pos1(1) = pos2(1); pos1(3) = pos2(3);
pos3(1) = pos4(1); pos3(3) = pos4(3);
pos5(1) = pos8(1); pos5(3) = pos8(3);
pos6(1) = pos9(1); pos6(3) = pos9(3);
pos7(1) = pos10(1); pos7(3) = pos10(3);
pos11(1) = pos14(1); pos11(3) = pos14(3);
pos12(1) = pos15(1); pos12(3) = pos15(3);
pos13(1) = pos16(1); pos13(3) = pos16(3);
pos17(1) = pos19(1); pos17(3) = pos19(3);
pos18(1) = pos20(1); pos18(3) = pos20(3);
pos21(1) = pos23(1); pos21(3) = pos23(3);
pos22(1) = pos24(1); pos22(3) = pos24(3);

set(plt(1),'Position',pos1)
set(plt(3),'Position',pos3)
set(plt(5),'Position',pos5)
set(plt(6),'Position',pos6)
set(plt(7),'Position',pos7)
set(plt(11),'Position',pos11)
set(plt(12),'Position',pos12)
set(plt(13),'Position',pos13)
set(plt(17),'Position',pos17)
set(plt(18),'Position',pos18)
set(plt(21),'Position',pos21)
set(plt(22),'Position',pos22)

print([exportpath,'Fig2ab_params_scatter_density_aligned_cell_lines'],'-r600','-dsvg','-painters');

%% Change in error rate as a function of change away from best params
load Data/MLSpike_ds

%% Find best parameters - get results in array format
clear all_ER_A all_ER_B all_ER_C all_PCC_A all_PCC_B all_PCC_C
for i = 1:21
    for j = 1:21
        all_ER_A(i,j) = results_A{i,j}.ER;
        all_ER_B(i,j) = results_B{i,j}.ER;
        all_ER_C(i,j) = results_C{i,j}.ER;
        
        all_PCC_A(i,j) = results_A{i,j}.PCC_FR;
        all_PCC_B(i,j) = results_B{i,j}.PCC_FR;
        all_PCC_C(i,j) = results_C{i,j}.PCC_FR;
    end
end

%% Plot ER shifted with respect to best params
clear bestERA bestERB bestERC minERA minERB minERC

for c = 1:21
    [bestERA(c),minERA(c)] = min(all_ER_A(c,:));
    [bestERB(c),minERB(c)] = min(all_ER_B(c,:));
    [bestERC(c),minERC(c)] = min(all_ER_C(c,:));
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 12 4]);

subplot(1,4,1)
for c = 1:21    
    plot((a_list-a_list(minERA(c))),all_ER_A(c,:)-all_ER_A(c,minERA(c)),'linewidth',widths.plot,'color',cmap_gt(c,:));
    hold all
end
plot([0,0],[0,1],'k--')
ylabel('\Delta ER')
axis square
xlabel('\Delta A')
title('A')

subplot(1,4,2)
for c = 1:21    
    plot((tau_list-tau_list(minERB(c))),all_ER_B(c,:)-all_ER_B(c,minERB(c)),'linewidth',widths.plot,'color',cmap_gt(c,:));
    hold all
end
plot([0,0],[0,1],'k--')
axis square
xlabel('\Delta tau')
title('tau')

subplot(1,4,3)
for c = 1:21    
    plot((sigma_list-sigma_list(minERC(c))),all_ER_C(c,:)-all_ER_C(c,minERC(c)),'linewidth',widths.plot,'color',cmap_gt(c,:));
    hold all
end
plot([0,0],[0,1],'k--')
axis square
xlabel('\Delta sigma')
title('sigma')

% Same for Suite2P
load Data/S2P_results_130218.mat
Thresholds = logspace(-1,2,13); 

ERs = S2P_results{3};

for c = 1:21
    [bestERS2P(c),minERS2P(c)] = min(ERs(c,:));
end

subplot(1,4,4)
for c = 1:21    
    plot((Thresholds-Thresholds(minERS2P(c))),ERs(c,:)-ERs(c,minERS2P(c)),'linewidth',widths.plot,'color',cmap_gt(c,:));
    hold all
end
plot([0,0],[0,1],'k--')
axis square
xlabel('\Delta Threshold')
title('Suite2P')

FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);
print([exportpath,'Fig2cd_delta_ER_MLSpike_S2P'],'-dpdf');












