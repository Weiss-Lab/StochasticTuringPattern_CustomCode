close all
clear
% maxk=5; %!!!!


load('processed_entropy__1e-2_0');
scalefactor=0.625;
allk=out.allk;

data.colors={...
    'r';...
    'r--';...
    'y';...
    'y--';...
    'g';...
    'g--';...
    'b';...
    'b--';...
    'c';...
    'c--';...
    'k';
    'k--';
    'ko'
    }

data.allfiles={...
    '1e-7'; ...
    '1e-6'; ...
    '1e-5.5'; ...
    '1e-5'; ...
    '1e-4.5'; ...
    '1e-4'; ...
    '1e-3.5'; ...
    '1e-3'; ...
    '1e-2.5'; ...
    '1e-2'; ...
    'Mix-red';...
    'pINV5';...
    'pASK201';...
    };

numfiles=size(data.allfiles,1);

collect_allS=zeros(size(data.allfiles,1),size(out.allS,2));

for whichfile=1:size(data.allfiles,1)
    load(char(input.savestring(whichfile)));
    collect_allS(whichfile,:)=out.allS;
end

xscale=scalefactor.*allk;


% 
% for filtind=1:size(collect_allS,1)
%     filt_collect_allS(filtind,:)=medfilt1(collect_allS(filtind,:));
% end

% %find k that gives max entropic measure
% [maxS,kmax]=max(collect_allS');
% figure;semilogx(IPTGconcs,kmax,'o')
% xlabel('IPTG concentration (M)');
% ylabel('scale of highest entropic measure (pixels)');
% figure;semilogx(IPTGconcs,kmax.*scalefactor,'o')
% xlabel('IPTG concentration (M)');
% ylabel('scale of highest entropic measure (\mum)');
% 
% %find k that gives max entropic measure
% [maxS,kmax]=max(filt_collect_allS');
% figure;semilogx(IPTGconcs,kmax,'o')
% xlabel('IPTG concentration (M)');
% ylabel('scale of highest entropic measure (pixels)');
% figure;semilogx(IPTGconcs,kmax.*scalefactor,'o')
% xlabel('IPTG concentration (M)');
% ylabel('scale of highest entropic measure (\mum)');

figure;plot(allk,collect_allS)
xlabel('scale (pixels)')
ylabel('entropic measure')
legend(data.allfiles);
% 
% figure;plot(1:maxk,filt_collect_allS)
% xlabel('scale (pixels)')
% ylabel('entropic measure')
% legend(data.allfiles);
% title('median filtered')
% 

% figure;plot(xscale,collect_allS)
% xlabel('scale (\mum)')
% ylabel('entropic measure')
% legend(data.allfiles);


    

    %calculate smooth fits
    xx=log(IPTGconcs);
    sxx=min(xx):((max(xx)-min(xx))/40):max(xx);
    
%     allsyy=zeros(numfiles,length(xx));
    for ifit=1:length(allk)
        yy=smooth(xx,collect_allS(1:10,ifit),'lowess');
        
        %     syy=interp1(xx,yy,sxx,'cubic');
        allsyy(ifit,:)=spline(xx,yy,sxx);
        %     figure;semilogx(IPTGconcs,collect_allS(1:10,1),'o',exp(sxx),syy)
    end
    
    
    
    
    
    
    
figure;
% hold on;
% for thisk=[]
    subplot(2,3,1);semilogx(IPTGconcs,collect_allS(1:10,1),'o',exp(sxx),allsyy(1,:),'k')%,...
%         'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],...
%         'XScale','log');title('scale=5.6 \mum')
    xlim([3.162e-008 0.0316]);
    subplot(2,3,2);semilogx(IPTGconcs,collect_allS(1:10,2),'o',exp(sxx),allsyy(2,:),'k');title('scale=10 \mum')
    xlim([3.162e-008 0.0316]);
    subplot(2,3,3);semilogx(IPTGconcs,collect_allS(1:10,3),'o',exp(sxx),allsyy(3,:),'k');title('scale=22 \mum')
    xlim([3.162e-008 0.0316]);
    subplot(2,3,4);semilogx(IPTGconcs,collect_allS(1:10,4),'o',exp(sxx),allsyy(4,:),'k');title('scale=54 \mum')
    xlim([3.162e-008 0.0316]);
    ylim([50 210]);
    subplot(2,3,5);semilogx(IPTGconcs,collect_allS(1:10,5),'o',exp(sxx),allsyy(5,:),'k');title('scale=72 \mum')
    xlim([3.162e-008 0.0316]);
    
    
    
    
    figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1,...
    'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],...
    'XScale','log',...
    'XMinorTick','off');
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim([5 15]);
box('on');
hold('all');

% Create semilogx
semilogx(IPTGconcs,collect_allS(1:10,1),'Parent',subplot1,'Marker','o','LineStyle','none');

% Create semilogx
semilogx(exp(sxx),allsyy(1,:),'Parent',subplot1,'Color',[0 0 0]);

% Create title
title('scale=5.6 \mum');
ylabel('entropy measure')

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1,'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],'XScale','log','XMinorTick','off');
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
box('on');
hold('all');

% Create semilogx
semilogx(IPTGconcs,collect_allS(1:10,2),'Parent',subplot2,'Marker','o','LineStyle','none');

% Create semilogx
semilogx(exp(sxx),allsyy(2,:),'Parent',subplot2,'Color',[0 0 0]);


% Create title
title('scale=10 \mum');

% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1,'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],'XScale','log','XMinorTick','off');
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
box('on');
hold('all');

% Create semilogx
semilogx(IPTGconcs,collect_allS(1:10,3),'Parent',subplot3,'Marker','o','LineStyle','none');

% Create semilogx
semilogx(exp(sxx),allsyy(3,:),'Parent',subplot3,'Color',[0 0 0]);

% Create title
title('scale=22 \mum');

xlabel('IPTG (M)')

% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1,'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],'XScale','log','XMinorTick','off');
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim([50 210]);
box('on');
hold('all');

% Create semilogx
semilogx(IPTGconcs,collect_allS(1:10,4),'Parent',subplot4,'Marker','o','LineStyle','none');

% Create semilogx
semilogx(exp(sxx),allsyy(4,:),'Parent',subplot4,'Color',[0 0 0]);

% Create title
title('scale=54 \mum');
xlabel('IPTG (M)')
ylabel('entropy measure')

% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1,'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],'XScale','log','XMinorTick','off');
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
box('on');
hold('all');

% Create semilogx
semilogx(IPTGconcs,collect_allS(1:10,5),'Parent',subplot5,'Marker','o','LineStyle','none');

% Create semilogx
semilogx(exp(sxx),allsyy(5,:),'Parent',subplot5,'Color',[0 0 0]);

% Create title
title('scale=72 \mum');
xlabel('IPTG (M)')
% 
% figure;
% % hold on;
% % for thisk=[]
%     subplot(2,3,1);plot(1:10,collect_allS(1:10,5),'o');title('scale=5')
%     subplot(2,3,2);plot(1:10,collect_allS(1:10,10),'o');title('scale=10')
%     subplot(2,3,3);plot(1:10,collect_allS(1:10,15),'o');title('scale=15')
%     subplot(2,3,4);plot(1:10,collect_allS(1:10,20),'o');title('scale=20')
%     subplot(2,3,5);plot(1:10,collect_allS(1:10,25),'o');title('scale=25')
%     subplot(2,3,6);plot(1:10,collect_allS(1:10,30),'o');title('scale=30')
% %     
% %     
% figure;
% % hold on;
% % for thisk=[]
%     subplot(2,3,1);plot(1:10,collect_allS(1:10,35),'o');title('scale=35')
%     subplot(2,3,2);plot(1:10,collect_allS(1:10,40),'o');title('scale=40')
%     subplot(2,3,3);plot(1:10,collect_allS(1:10,45),'o');title('scale=45')
%     subplot(2,3,4);plot(1:10,collect_allS(1:10,50),'o');title('scale=50')
%     subplot(2,3,5);plot(1:10,collect_allS(1:10,55),'o');title('scale=55')
%     subplot(2,3,6);plot(1:10,collect_allS(1:10,60),'o');title('scale=60')