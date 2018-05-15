% a code to characterize pattern domains
tic
clear all;
close all;

disp([' ']);
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
disp(['new processing ...']);
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);

x1 = 1 ;
chan = 3 ;
mst = 1; %1--mean, 2--std. 3--total. of red spots.
IPTG=0;
%IPTG = [ -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0]' ;

%IPTG = [ -10 -9 -8 -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0]' ;

num_images = size(IPTG,1) ;

%num_images = 10 ;


%dirname1 = ['C:\Users\Ting\Pictures\Exported Images\IPTG_pub\'];
dirname1 = ['.\FIG2v3.bmp'];


Mat_Dom_mean = [] ;
Mat_Dom_std  = [] ;
Mat_Dom_sum  = [] ;

Mat_WDom_mean = [] ;
Mat_WDom_std  = [] ;
Mat_WDom_sum  = [] ;

Moran_all = [] ;
Area_all = [] ;
wxx_sum_all = [];
wij_sum_all = [];

fid = fopen('rfp_area.tif', 'wt');
fclose(fid); 

for img_dex = 1: num_images
   % close all;
    %load the image
    x2 = IPTG(img_dex);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    x2 = -105 ; % high red
%    x2 = -2 ; % low red 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %disp(['img_dex=' num2str(img_dex)]);
    
   
    
    %PicX = imread([dirname1, dirname2]);
    PicX = imread([dirname1]);
     
%    figure, imshow(PicX,'InitialMagnification', 'fit');
    %convert to a binary one.
    Pic_bw = im2bw(PicX, 0.1) ; % attention: choosing the right threshold is critical!

    %count and label the domains
    [labeled,numObjects] = bwlabel(Pic_bw,4);
    numObjects; % Count all distinct objects in the image.
    %disp(['original numObjects= ' num2str(numObjects)]);

    %show the domains using pseudocolor image
    RGB_label = label2rgb(labeled, @spring, 'c', 'shuffle');
    figure, imshow(RGB_label,'InitialMagnification', 'fit');

    %analyze domain characters
    domaindata = regionprops(labeled,'basic'); %which returns: area, centroid, and 'bounding box'

    alldomains = [domaindata.Area];
    max_area = max(alldomains); % Find the maximum area of all the grains.
    disp(['max_area= ' num2str(max_area)]);

    % %show domains with centroids
    % figure, imshow(RGB_label,'InitialMagnification', 'fit');
     centroids = cat(1, domaindata.Centroid);
%     hold on;
%     plot(centroids(:,1), centroids(:,2), 'b.', 'linewidth', 2.0);

    %filter small domains
    thres_area = 50 ;
    idx = find([domaindata.Area] > thres_area );
    %idx = find([domaindata.Area] < thres_area );
    Judge_domain = ismember(labeled, idx);

    %show the image after filtering
    re_labeled = labeled .* Judge_domain ;
    figure, imshow(re_labeled,'InitialMagnification', 'fit');

    Pic_bwFTD = Pic_bw .* Judge_domain ;
    [labeled3,numObjects3] = bwlabel(Pic_bwFTD,8);

    disp(['new numObjects= ' num2str(numObjects3)]);

    
    %show the domains in the filtered image using pseudocolor image
    RGB_label3 = label2rgb(labeled3, @spring, 'c', 'shuffle');
%    figure, imshow(RGB_label3,'InitialMagnification', 'fit');
%    nametitle = [num2str(x1), 'e', num2str(x2)];
%    title(nametitle);
    
    %analyze domain characters for filtered image
    domaindata3 = regionprops(labeled3,'basic'); %which returns: area, centroid, and 'bounding box'

    %show the centroids in the filtered image
  %figure, imshow(RGB_label3,'InitialMagnification', 'fit');

    centroids = cat(1, domaindata3.Centroid);
  %hold on;
  %plot(centroids(:,1), centroids(:,2), 'k.', 'linewidth', 1.5);

  
  
  
%[domaindata3.Area]  
RSpot = [domaindata3.Area];
RSpot1 = sum(RSpot);
RSpot2 = sum(RSpot.*RSpot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X_bar
MI_m = mean2(Pic_bwFTD) ;
%N
MI_n  = size(Pic_bwFTD,2).*size(Pic_bwFTD,1) ;
%reshape the matrix to vector
MI_FP1 = reshape(Pic_bwFTD, [],1);
MI_LB3 = reshape(labeled3, [],1);
%denominator
MI_I_dom = sum(sum( (MI_FP1-MI_m).*(MI_FP1-MI_m) ));

wei1=1.0;
wei2=1.0;

Wij_sum = wei1*RSpot2 + wei2*((MI_n-RSpot1).^2) ;
Wxx_sum = wei1*RSpot2.*((1-MI_m).^2) + wei2*((MI_n-RSpot1).^2) .* (MI_m.^2) ;

Moran_I = MI_n.*Wxx_sum./Wij_sum./MI_I_dom;

Moran_I ;

Moran_all = [ Moran_all, Moran_I ];
Area_all = [ Area_all, RSpot1 ];
wij_sum_all = [wij_sum_all, Wij_sum ];
wxx_sum_all = [wxx_sum_all, Wxx_sum ];

    % %power spectrum for the filtered binary image
    % tic 
    % F1=fft2(labeled3,2048,2048); 
    % F2=fftshift(F1); 
    % PS=abs(F2).*abs(F2);
    % figure, imagesc(log(PS));
    % colorbar;
    % set(gca,'clim',[5 30]);
    % toc

    %%%%%%%%%%%%%%%%%%%%%%%to plot domain distribution
  figure(99), hist( [domaindata3.Area], 20) ;
  title('Area distribution');
    
    %figure(200),
    %[num1,xout1] = hist([domaindata3.Area],20);
    %ntot1=sum(num1);
    %barhand1 = bar(xout1, num1./ntot1);
    
    
%     figure, hist(log(sqrt([domaindata3.Area])),[1.5:0.1:3.8]); 
%     title('Log(sqrt(area)) distribution');
%     
%     %area =r^2. ==> log(area) = 2*log(r)
%     figure, hist(log([domaindata3.Area]),20); %similar to log(effective radius)
%     title('Log(area) distribution');
%     figure, hist(sqrt([domaindata3.Area]),20);
%     title('sqrt(area) dist');
%     figure, hist( ([domaindata3.Area].*[domaindata3.Area]),20);
%     title('(area.*area) distribution');
%     figure, hist(log([domaindata3.Area].*[domaindata3.Area]),20);
%     title('log(area.*area) distribution');
    %%%%%%%%%%%%%%%%%%%%%%%to plot weighted domain distribution
    size_lab3 = size(labeled3,1)*size(labeled3,2);
    vec_lab3 = reshape(labeled3,[],size_lab3);
 
 %!!!   
 %  RSpot = [domaindata3.Area];
 
    vec_area3 = [domaindata3.Area];
    index_3 = find(vec_lab3>0);
    vec_lab3_nonzero = vec_lab3(index_3) ;

    
    weighted_area = vec_area3(vec_lab3_nonzero) ;
    %area histogram is weighted by their area. 
    
    %figure(99), hist(weighted_area, 75)

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % added_area3 = [0,vec_area3];
    % %scaning over the vector and count the corresponding area.
    % weighted_area = added_area3(1+vec_lab3) ;
    % %remove all 'zero' area.
    % index_3 = find(weighted_area);
    % weighted_area_corr = weighted_area(index_3) ;
    % 
    % %weighted_area(weighted_area == 0) = [] ;

%     %logarithmic weighted area
%     figure,
%     hist(log(weighted_area+1e-5),3:0.2:7.5);
%     [num2,xout2] = hist(log(weighted_area+1e-5),3:0.1:7.5);
%     ntot2=sum(num2);
%     barhand2 = bar(xout2, num2./ntot2);
%     axis([2.9 7.5 0 0.15]);
%     
%     
%     
%     %%%%%%weighted area
%     %figure(110),
%     figure,
%     histep = 60;
%     nstep = 25 ;
%     xhis=0.5*histep:histep:nstep*histep;
%     
%     [num7,xout7] = hist(weighted_area,xhis);
%     ntot=sum(num7);   
%     if (img_dex~=1)
%         hold on;
%     end
%         
% %   barhand = bar(xout, num./ntot);
% %   axis([0 1500 0 0.22])
%    barhand = bar(xout7, num7);
%    axis([0 1500 0 10000])
%    set(barhand,'EdgeColor','w','FaceColor', [.02 .839 .98], 'BarWidth',1);
%    title('Area-weighted Domain Size histogram','fontsize',16);
% 
     mean_weighted_area = mean(weighted_area) ;
     std_weighted_area = std(weighted_area) ;
%     disp(['weighed_area: mean, std']);
%     [mean_weighted_area, std_weighted_area]
%     disp([' ']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %domain size statistics
    Dom_mean = mean([domaindata3.Area]);
    Dom_std = std([domaindata3.Area]);
    Dom_sum = sum([domaindata3.Area]);
    
    Mat_Dom_mean = [Mat_Dom_mean, Dom_mean];
    Mat_Dom_std  = [Mat_Dom_std,  Dom_std] ;
    Mat_Dom_sum  = [Mat_Dom_sum,  Dom_sum] ;

    %!!!   
 %  RSpot2 = RSpot1.^2;
 
    SW_area = sum([domaindata3.Area].*[domaindata3.Area]);
    
    %equivalent to: SW_area = sum(weighted_area);
       
    Mat_WDom_mean = [Mat_WDom_mean, mean_weighted_area] ;
    Mat_WDom_std  = [Mat_WDom_std,  mean_weighted_area] ;
    Mat_WDom_sum  = [Mat_WDom_sum,  SW_area];
    
    disp(['[Dom_mean, std, sum]= ' num2str(Dom_mean), ', ' num2str(Dom_std), ', ' num2str(Dom_sum)]);
    
    disp(['lined:']);
    disp([num2str(x1), 'e', num2str(x2), '  ' num2str(Dom_mean), ' ' num2str(Dom_std), ' ' num2str(Dom_sum)]);

    fid = fopen('rfp_area.txt', 'a+');
    fprintf(fid, '%.1fe%d %12.5f %12.5f %12.5f\n', x1, x2, Dom_mean, Dom_std, Dom_sum);
    fclose(fid);

%     %domain size histogram
%     %figure,
%     %hist([domaindata3.Area], 25);
%     %NDh = findobj(gca,'Type', 'patch');
%     %set(NDh,'EdgeColor','w','FaceColor', [.02 .839 .98]);
% 
%     histep = 120;
%     nstep = 30 ;
%     xhis=0.5*histep:histep:nstep*histep;
%     %figure,hist(log(allgrains+1e-5),xhis);
%     [num,xout] = hist([domaindata3.Area],xhis);
%     ntot=sum(num);
%     figure, barhand = bar(xout, num./ntot);
% 
%     set(barhand,'EdgeColor','w','FaceColor', [.02 .839 .98], 'BarWidth',1);
%     title('Domain size statistics','fontsize',16);
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %nearest distance statistics
%     CoMaxis = cat(1, domaindata3.Centroid);
% 
%     CMLen = size(CoMaxis,1) ;
% 
%     minDD = zeros(CMLen,1);
%     for si = 1: CMLen
%         [minDD(si),Ds] = MinDomainDist(CoMaxis(si,:),CoMaxis);
%     end
% 
%     %nearest distance histogram
%     figure
%     hist(minDD, 25, 'BarWidth',1);
% 
%     minDD_mean = mean(minDD);
%     minDD_std  = std(minDD);
%     disp(['[minDD_mean, std]= ', num2str(minDD_mean), ', ' num2str(minDD_std)]);
% 
%     NDh = findobj(gca,'Type', 'patch');
%     set(NDh,'EdgeColor','w','FaceColor', [.02 .839 .98]);
%     title('Nearest distance','fontsize',16);

end

%plot figures
%domain area dist.
Mat_Dom_mss=[Mat_Dom_mean', Mat_Dom_std', Mat_Dom_sum'];

%weighed domain area dist.
Mat_WDom_mss=[Mat_WDom_mean', Mat_WDom_std', Mat_WDom_sum'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to plot weighted area statistics.
%figure(55),
figure,
%subplot(2,1,2)
if chan==2,
    lincha='g';
else
    lincha='r';
end
    
plotobj = plot(IPTG, Mat_WDom_mss(:,3),'o');
 title('Total size of red domains vs IPTG','fontname','arial','fontsize',16);
% xlabel('IPTG','fontname','arial', 'fontsize',16);
% ylabel('Total Red spot size', 'fontname','arial','fontsize',16);
set(gca, 'fontsize',14);
set(plotobj,...
                'LineWidth',2,...
                'MarkerEdgeColor',lincha,...
                'MarkerFaceColor',lincha,...
                'MarkerSize',8)

%%%%%%%%%%%%fitting Fig(55)
a = 4.0e7 ;
b = 1 ;
c = 1.0*10^(-4.0) ;
d = -1.2;
f = 0.7e7 ;

x0=-7.2:0.1:-1.8;
x = 10.^(x0);

y = a - a./(b+(x./c).^d)+ f;
    

hold on ;
plot(x0,y,'g-','linewidth',1.0);

% hold on ;
% errorbar(IPTG, XFP2(:,mors),XFP2(:,2),'k','linestyle','none','linewidth',1.0);
% 

axis([-7.1 -1.9 0.7*min(y) 1.5*max(y)]);
            
            
            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot area dist.
% %figure(3),
% figure,
% 
% %subplot(2,1,2)
% if chan==2,
%     lincha='g';
% else
%     lincha='r';
% end
%     
% plotobj = plot(IPTG, Mat_Dom_mss(:,mst),'o');
%  title('mean red domain size vs IPTG','fontname','arial','fontsize',16);
% % xlabel('IPTG','fontname','arial', 'fontsize',16);
% % ylabel('Total Red spot size', 'fontname','arial','fontsize',16);
% set(gca, 'fontsize',14);
% set(plotobj,...
%                 'LineWidth',2,...
%                 'MarkerEdgeColor',lincha,...
%                 'MarkerFaceColor',lincha,...
%                 'MarkerSize',8)
%             
%             
%             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %curve fitting
% if mst==3
% a = 6.0e4 ;
% b = 1 ;
% c = 1.0*10^(-4.4) ;
% d = -1.5;
% f = 3.7e4 ;
% elseif mst==1
% a = 2e2 ;
% b = 1 ;
% c = 1.0*10^(-4.4) ;
% d = -1.5;
% f = 1e2 ;
% else
%     ;
% end    
%     
% a = 1.5e2 ;
% b = 1 ;
% c = 1.0*10^(-3.7) ;
% d = -1.5;
% f = 1e2 ;
% 
% x0=-7.2:0.1:-1.8;
% x = 10.^(x0);
% 
% y = a - a./(b+(x./c).^d)+ f;
%     
% 
% hold on ;
% plot(x0,y,'k-','linewidth',1.0);
% 
% % hold on ;
% % errorbar(IPTG, XFP2(:,mors),XFP2(:,2),'k','linestyle','none','linewidth',1.0);
% % 
% 
% axis([-7.1 -1.9 0.7*min(y) 1.5*max(y)]);
% 
% %axis([-7.1 -1.9 0.0 110000]);

figure(101);
subplot(2,2,1),
plot(IPTG, Moran_all,'g-o');
axis([-7.1 -1.9 0.0 1.2*max(Moran_all)]);
title('Moran I vs. IPTG');
subplot(2,2,2),
plot(IPTG, Area_all,'b-o');
axis([-7.1 -1.9 0.0 1.2*max(Area_all)]);
title('Total Red Area vs. IPTG');
subplot(2,2,3),
plot(IPTG, wij_sum_all,'k-o');
axis([-7.1 -1.9 0.0 1.2*max(wij_sum_all)]);
title('wij-sum(collectivity) vs. IPTG');
subplot(2,2,4),
plot(IPTG, wxx_sum_all,'r-o');
axis([-7.1 -1.9 0.0 1.2*max(wxx_sum_all)]);
title('wxx-sum vs. IPTG');

