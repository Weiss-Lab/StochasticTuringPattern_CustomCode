close all
clear

%options
zscores=0; %set to 1 to calculate E(moran), std(moran) and Z-scores
runs=1000;

input.allims={...
    '1e-7.0_c3.tif'; ...
    '1e-6.0_c3.tif'; ...
    '1e-5.5_c3.tif'; ...
    '1e-5.0_c3.tif'; ...
    '1e-4.5_c3.tif'; ...
    '1e-4.0_c3.tif'; ...
    '1e-3.5_c3.tif'; ...
    '1e-3.0_c3.tif'; ...
    '1e-2.5_c3.tif'; ...
    '1e-2.0_c3.tif'; ...
%     'Mix-red.tif';...
%     'pINV5.tif';...
%     'pASK201.tif';...
    };

numims=size(input.allims,1)

IPTGconcs=[10^-7 10^-6 10^-5.5 10^-5 10^-4.5 10^-4 10^-3.5 10^-3 10^-2.5 10^-2];
allmoran=zeros(1,numims)
allrandmoran=zeros(numims,runs);
meansofar=zeros(numims,runs);
stdsofar=zeros(numims,runs);
z=zeros(1,numims);
meansuppmoran=zeros(1,numims);

for whichim=1:numims
    fname=char(input.allims(whichim));
    I=imread(fname);
    %     It=I;%(50:120,50:120);
%    imbw=im2bw((2^16/2^12).*I,1*graythresh((2^16/2^12).*I));

    imbw=im2bw(I);
    
%     figure,
%     imshow(imbw);
%     
    Ifinal=imbw;
    
    momo=moranI(Ifinal);
    allmoran(whichim)=momo;
    
    %     figure;imshow(Ifinal)%!!!!!!!!!
    % end %!!!!!!!!!!!!
    
    
    %calculate E(moran) and std(moran) for each image hist
    if zscores==1
        for thisrun=1:runs
            randim=impermute(Ifinal);
            allrandmoran(whichim,thisrun)=moranI(randim);
            lastmean=mean(allrandmoran(whichim,1:thisrun));
            laststd=std(allrandmoran(whichim,1:thisrun));
            meansofar(whichim,thisrun)=lastmean;
            stdsofar(whichim,thisrun)=laststd;
        end
        z(whichim)=(momo-lastmean)/laststd;
        meansuppmoran(whichim)=momo-lastmean;
    end
end
if zscores==1
    figure;plot(meansofar');title('mean vs. number of runs')
    figure;plot(stdsofar');title('std vs. number of runs')
    figure;bar(meansuppmoran);title('moran I - E(moran I)')
    figure;bar(z);title('Z score')
end
%
% %calculate E(moran) and std(moran) for each image hist
% runs=100;
% allrandmoran=zeros(1,runs);
% meansofar=zeros(1,runs);
% stdsofar=zeros(1,runs);
% for thisrun=1:runs
%     randim=impermute(I);
%     allrandmoran(thisrun)=moranI(randim);
%     meansofar(thisrun)=mean(allrandmoran(1:thisrun));
%     stdsofar(thisrun)=std(allrandmoran(1:thisrun));
%     waitbar(thisrun/runs)
% end
% figure;plot(meansofar)
% figure;plot(stdsofar)

figure;semilogx(IPTGconcs,allmoran(1:10))
figure;bar(allmoran);title('moran I')
% legend(input.allims)

%fit Moran vs. IPTG curve
xx=log(IPTGconcs);
sxx=min(xx):((max(xx)-min(xx))/40):max(xx);
yy=smooth(xx,allmoran(1:10),'lowess');

yy2=smooth(xx,yy,'lowess');
yy3=smooth(xx,yy2,'lowess');

%syy=interp1(xx,yy,sxx,'cubic');
syy=interp1(xx,yy3,sxx,'cubic');
% figure;semilogx(IPTGconcs,allmoran(1:10),'o',exp(sxx),syy)

    
    
% Create figure
figure1 = figure;
% Create axes
axes('Parent',figure1,'XScale','log','XMinorTick','on');
% Uncomment the following line to preserve the X-limits of the axes
% xlim([3.162e-008 0.0316]);
box('on');
hold('all');
% Create semilogx

semilogx(IPTGconcs,allmoran(1:10),'o',exp(sxx),syy,'k');


% semilogx(IPTGconcs,allmoran(1:10),'Marker','o','LineStyle','none');
% semilogx(exp(sxx),syy,'Color',[0 0 0],'LineStyle','none');
% Create xlabel
xlabel({'IPTG (M)'},'FontSize',16);
% Create ylabel
ylabel({'Moran I'},'FontSize',16);










% Create figure
figure1 = figure;
% Create axes
axes('Parent',figure1,'XTick',[1e-007 1e-006 1e-005 0.0001 0.001 0.01],...
    'XScale','log',...
    'XMinorTick','on',...
    'Position',[0.1845 0.1938 0.7205 0.7313],...
    'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
xlim([3.162e-008 0.0316]);
box('on');
hold('all');
% Create semilogx
semilogx(IPTGconcs,allmoran(1:10),'MarkerSize',10,'Marker','o','LineWidth',2,...
    'LineStyle','none');
% Create semilogx
semilogx(exp(sxx),syy,'LineWidth',2,'Color',[0 0 0]);
% Create xlabel
xlabel({'IPTG (M)'},'FontSize',16);
% Create ylabel
ylabel({'Moran I'},'FontSize',16);





