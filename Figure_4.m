% Figure_4.m - Event rate distributions and silent cells
%
% Code to produce figure 4 of the 'Pitfalls of deconvolution' paper 2020 
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

[ncells,nt] = size(dat.timeSeriesArrayHash.value{1,2}.valueMatrix);

meths = 1:9; 
nmeths = numel(meths);
meth_c = 1:6; % Continuous methods.
meth_s = 7:9; % Spike inference methods. 

%% Event rate calculation
clear ER_g
for j = 1:nmeths
    load(['Data/deconv_nine_examples/',data_ID{meths(j)},'.mat'])
    data = eval(data_ID{meths(j)});
    
    
    if ismember(meths(j),meth_c) % deconvolved calcium methods
        for c = 1:ncells
            c;
            this_c = data(c,:);
            
            % Threshold based on std of residual noise (remove smooth version of data
            % first)
            smooth_c = conv(this_c,ones(4,1),'same')/4;
            resid = this_c - smooth_c;
            
            % Clear zeros
            resid = resid(find(resid));
            
            sig = std(resid);
            level = mean(this_c) + 3*sig;
            
            found_events = find(this_c>=level);
            
            % Event rate (Hz)
            ER_g(j,c) = 7 * numel(found_events)/nt;
        end
        
    elseif ismember(meths(j),[7,9]) % Spike inference methods - count events
        
        for c = 1:ncells
            this_c = data(c,:);
            
            ER_g(j,c) = 7 * (numel(find(this_c))/nt);
        end
    end
    
   if meths(j) == 8 % MLSpike as it returns spike counts
        
        for c = 1:ncells
            this_c = data(c,:);
            
            ER_g(j,c) = 7 * sum(this_c)/nt;
        end
    end
end


%% Fix buggy LZero result where all elements are 1 instead of 0 for silent cells
silly_cells = find(ER_g(9,:)>3.5); 
ER_g(9,silly_cells) = 0;           

%% Plot event/spike rate results

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 8 16]);

mx = max(ER_g(:));
mn = min(ER_g(:));

bins = linspace(mn,mx,100);
clear ax
for j = 1:nmeths
    ax(j) = subplot(nmeths,1,j);
    [F] = hist(ER_g(j,:),bins);
    
    bar(bins,F,'facecolor',cmap_ca(j,:),'edgecolor',cmap_ca(j,:))
    ylabel(meth_names_paper{meths(j)})
    set(get(gca,'ylabel'),'rotation',0, 'HorizontalAlignment','right')
    pbaspect([5 1 1])
%     axis off
end
linkaxes(ax,'x')
xlim([0,4.5]); 
xlabel('Event rate (Hz)')
    
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig_4a_event_rate'],'-dpdf'); 

%% Silent cells (using Peron et al 2014/O'Connor et al 2010 0.0083Hz definition) 
% O'Connor 2010 found ~13% silent cells across cortical layers. 
% 44 of cells had FRs < 1 Hz. L2/3 this figure was closer to 60%
% L2/3 - mean 3Hz, median 0.18Hz. ~26% silent cells

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 10 10]);

for j = 1:nmeths
    subplot(3,3,j)
    not_silent = zeros(ncells,1);
    not_silent(find(ER_g(j,:)>0.0083)) = 1;
    
    pie(hist(not_silent,2))
    colormap([0,0,0;.5,.5,.5])
    title(meth_names_paper{meths(j)})
end
axis off
%  suptitle('Silent (black) vs active (gray) cells')
 
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);
print([exportpath,'Fig_4b_silent_cells'],'-dsvg')