% Figure properties for 'pitfalls of deconvolution' paper figures

format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
color = 'rgb';
dpi = 600;
fontsize = 8;
fontname = 'Arial';
M = 3; % marker size for univariate scatter plots
sym = 'o';  % markers for scatters and strip plots

Units = 'centimeters';

% line widths
widths.plot = 0.75;
widths.error = 0.5;
widths.axis = 0.5;

% panel sizes
figsize.square = [5 5];
figsize.linkedscatter = [3 5];
figsize.dists = [8 4];
figsize.wide = [8 2];
figsize.retina_square = [10 10];
figsize.retina_dists = [16 8];
figsize.retina_wide = [16 4];

% colours for ground truth plots
cmap_gt = varycolor(21);

% colours for deconvolution of real data
cmap = brewermap(9,'Set1'); 
cmap(6,:) = [0.8,0.8,0.2];  % changing yellow as too light to see
cmap_ca = circshift(cmap,[1,0]);

% exportpath
if ~exist('Fig_panels', 'dir')
    mkdir('Fig_panels');
end
relpath = pwd;
exportpath = [relpath,'/Fig_panels/'];

