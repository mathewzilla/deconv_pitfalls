function myeb(Y,varargin);
%
% myeb(Y,varargin);
%
% This function makes nice coloured, shaded error bars. Exactly what
% it does depends on Y, and on whether you give it one or more inputs.
%
% If you only pass it Y, and no other arguments, it assuemd you're
% giving it raw data.
%
%		myeb(Raw_Data)
%
% 	.) if Y is 2D array, it will then plot mean(Y) with errorbars given
% 	by std(Y). In this case there is only one mean vector with its
% 	errorbars.
%
%	.) if Y is 3D array, it will plot size(Y,3) lines with the
%	associated errorbars. Line k will be mean(Y(:,:,k)) with errorbars
%	given by std(Y(:,:,k))
%
% If you pass it 2 arguments, each has to be at most 2D.
%
%		myeb(mu,std)
%
% 	.) if mu and std are 1D, it just plots one line given by mu with a
% 	shaded region given by std.
%
%	.) if mu and std are 2D, it will plot size(Y,2) lines in the
%	standard sequence of colours; each line mu(:,k) will have a shaded
%	region in the same colour, but less saturated given by std(:,k)
%
%
% Quentin Huys, 2007
% Center for Theoretical Neuroscience, Columbia University
% Email: qhuys [at] n e u r o theory [dot] columbia.edu
% (just get rid of the spaces, replace [at] with @ and [dot] with .)
%
% Edited by Mathew H. Evans, 2015
%
% If you pass it 3 arguments,
%
%       myeb(xtick,mu,std)
%
%   the function works like errorbar, using the
%
% If you pass it 4 arguments,
%
%       myeb(xtick,mu,std,color)
%
%   Specify colour with 3 element rgb code. e.g. [0 0 0] for black. Add a
%   4th element to set the shading as a lighter version of
%   the line color e.g. [0 0 0 0] as gray
%
% If you pass it 5 arguments,
%
%       myeb(xtick,mu,std,line_color,shading_color)
%
%   Specify line and shading colours with 3 element rgb code. e.g. [0 0 0] for black. [.5 .5 .5] for gray.

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;


if length(varargin)==0;
    
    if length(size(Y))==2
        if size(Y,2)>=8;
            col = distinguishable_colors(size(Y,2));
            ccol=col+.8; ccol(ccol>1)=1;
        end
        m=mean(Y);
        s=std(Y);
        ind1=1:length(m);
        ind2=ind1(end:-1:1);
        hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
        set(h,'edgecolor',.6*ones(1,3))
        plot(ind1,m,'linewidth',2)
        hold off
    elseif length(size(Y))>2
        if size(Y,2)>=8;
            col = distinguishable_colors(size(Y,2));
            ccol=col+.8; ccol(ccol>1)=1;
        end
        cla; hold on;
        ind1=1:size(Y,2);
        ind2=ind1(end:-1:1);
        if size(Y,3)>8; col=jet(size(Y,3));ccol=col+.8; ccol(ccol>1)=1;end
        for k=1:size(Y,3)
            m=mean(Y(:,:,k));
            s=std(Y(:,:,k));
            h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],ccol(k,:));
            set(h,'edgecolor',ccol(k,:))
        end
        for k=1:size(Y,3)
            m=mean(Y(:,:,k));
            s=std(Y(:,:,k));
            plot(ind1,m,'linewidth',1,'color',col(k,:))
        end
        hold off
    end
    
elseif length(varargin)==1;
    if size(Y,2)>=8;
        col = distinguishable_colors(size(Y,2));
        ccol=col+.8; ccol(ccol>1)=1;
    end
    m=Y;
    s=varargin{1};
    if length(size(Y))>2;
    elseif min(size(Y))==1;
        if size(m,1)>1; m=m';s=s';end
        ind1=1:length(m);
        ind2=ind1(end:-1:1);
        hold on; h=fill([ind1 ind2],[m-s m(ind2)+s(ind2)],.6*ones(1,3));
        set(h,'edgecolor',.6*ones(1,3))
        plot(ind1,m,'linewidth',2)
        hold off
    else
        ind1=(1:size(Y,1));
        ind2=ind1(end:-1:1);
        cla; hold on;
        if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
        for k=1:size(Y,2)
            mm=m(:,k)';
            ss=s(:,k)';
            h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
            set(h,'edgecolor',ccol(k,:))
        end
        for k=1:size(Y,2);
            mm=m(:,k)';
            ss=s(:,k)';
            plot(ind1,mm,'linewidth',1,'color',col(k,:))
        end
        hold off
    end
    
    
elseif length(varargin)==2;
    
    xt = Y;
    m=varargin{1};
    s=varargin{2};
    if length(size(m))>2;
    elseif min(size(m))==1;
        if size(m,1)>1; m=m';s=s';end
        ind1=xt;
        ind2=fliplr(xt);
        ind1b = 1:length(xt);
        ind2b=ind1b(end:-1:1);
        hold on; h=fill([ind1 ind2],[m-s m(ind2b)+s(ind2b)],.6*ones(1,3));
        set(h,'edgecolor',.6*ones(1,3))
        plot(ind1,m,'k','linewidth',2)
        hold off
    else
        % This part hasn't been changed from above, so won't work
        ind1=xt
        ind2=ind1(end:-1:1);
        cla; hold on;
        if size(m,2)>8; col=jet(size(m,2));ccol=col+.8; ccol(ccol>1)=1;end
        for k=1:size(Y,2)
            mm=m(:,k)';
            ss=s(:,k)';
            h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
            set(h,'edgecolor',ccol(k,:))
        end
        for k=1:size(m,2);
            mm=m(:,k)';
            ss=s(:,k)';
            plot(ind1,mm,'linewidth',2,'color',col(k,:))
        end
        hold off
    end
    
elseif length(varargin)==3;
    xt = Y;
    m=varargin{1};
    s=varargin{2};
    color = varargin{3};
    colour = color(1:3);
    if numel(color) == 4;
        barcolour = colour + 0.8;
        barcolour(barcolour>1) = 1;
    else
        barcolour = colour;
    end
    if length(size(m))>2;
    elseif min(size(m))==1;
        if size(m,1)>1; m=m';s=s';end
        ind1=xt;
        ind2=fliplr(xt);
        ind1b = 1:length(xt);
        ind2b=ind1b(end:-1:1);
        hold on; h=fill([ind1 ind2],[m-s m(ind2b)+s(ind2b)],barcolour);
        set(h,'edgecolor',barcolour)
        plot(ind1,m,'color',colour,'linewidth',2)
        hold off
    else
        % This part hasn't been changed from above, so won't work
        ind1=(1:size(Y,1));
        ind2=ind1(end:-1:1);
        cla; hold on;
        if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
        for k=1:size(Y,2)
            mm=m(:,k)';
            ss=s(:,k)';
            h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
            set(h,'edgecolor',ccol(k,:))
        end
        for k=1:size(Y,2);
            mm=m(:,k)';
            ss=s(:,k)';
            plot(ind1,mm,'linewidth',2,'color',col(k,:))
        end
        hold off
    end
    
elseif length(varargin)==4;
    
    xt = Y;
    m=varargin{1};
    s=varargin{2};
    color = varargin{3};
    colour = color(1:3);
    barcolour = varargin{4};
    
    if length(size(m))>2;
    elseif min(size(m))==1;
        if size(m,1)>1; m=m';s=s';end
        ind1=xt;
        ind2=fliplr(xt);
        ind1b = 1:length(xt);
        ind2b=ind1b(end:-1:1);
        hold on; h=fill([ind1 ind2],[m-s m(ind2b)+s(ind2b)],barcolour);
        set(h,'edgecolor',barcolour)
        plot(ind1,m,'color',colour,'linewidth',2)
        hold off
    else
        % This part hasn't been changed from above, so won't work
        ind1=(1:size(Y,1));
        ind2=ind1(end:-1:1);
        cla; hold on;
        if size(Y,2)>8; col=jet(size(Y,2));ccol=col+.8; ccol(ccol>1)=1;end
        for k=1:size(Y,2)
            mm=m(:,k)';
            ss=s(:,k)';
            h=fill([ind1 ind2],[mm-ss mm(ind2)+ss(ind2)],ccol(k,:));
            set(h,'edgecolor',ccol(k,:))
        end
        for k=1:size(Y,2);
            mm=m(:,k)';
            ss=s(:,k)';
            plot(ind1,mm,'linewidth',2,'color',col(k,:))
        end
        hold off
    end
end