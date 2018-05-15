
% To get spatially averaged mean and variance for GFP and RFP
%happy new year!

clear all;
close all ;

disp(['new run starts:']);
%sample size
Isize = 1 ;

%channel number
x1 = 1 ;
chan = 2 ;
IPTG = [ -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0]' ;
%IPTG = [-10.0 -9.0 -8.0 -7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0]' ;


num_images = size(IPTG,1) ;


%dirname1 = ['C:\Users\Ting\Pictures\Exported Images\IPTG_pub\Absolute\'];

%2016 new directory:
dirname1 = ['D:\Turing_project\Figure4\IPTG_pub\Absolute\'];


%dirname2 = [num2str(x1), 'e', num2str(x2), '_c', num2str(x3), '.TIF' ];
%filename = [dirname1, dirname2];


XFP = [];
std_XFP=[];

for img_dex = 1: num_images
    
    m_FP=[];
    v_FP=[];

    x2 = IPTG(img_dex);
    
    disp(['img_dex=' num2str(img_dex)]);
 
    %read in the file
    %Imgx = imread(fname);
    
    if ( mod(x2,1) == 0 )
        dirname2 = [num2str(x1), 'e', num2str(x2), '.0', '_c', num2str(chan), '.TIF' ];
    else
        dirname2 = [num2str(x1), 'e', num2str(x2), '_c', num2str(chan), '.TIF' ];
    end
    fname = [dirname1, dirname2];
    
    Imgx = imread(fname);
    %convert the matrix to a vector
    Im_vec = reshape(Imgx,size(Imgx,1)*size(Imgx,2),[]);
            
    %convert to double type
    Im_vec = double(Im_vec);

    %background substraction
%    if chan==2
%        background = 180.9862;
%    elseif chan==3
%        background = 184.6638;
%    end
%    
%   Im_vec = Im_vec - background ;
    
    %calculate mean and stdandard deviation
    m_FP    = mean(Im_vec);
    std_FP  = std(Im_vec);
    
    disp(fname);
    [m_FP std_FP]
  
    XFP = [XFP, m_FP];
    std_XFP = [std_XFP, std_FP];
    
end

%subtract background noise
background = 40 ; 
XFP = XFP - background ;

XFP_ms=[XFP',std_XFP'];


figure(3), 
%subplot(2,1,1)
if chan==2,
    lincha='g';
else
    lincha='r';
end
    
plotobj = plot(IPTG, XFP_ms(:,1),'o');
%title('Averaged GFP vs IPTG','fontname','arial','fontsize',16);
% xlabel('IPTG','fontname','arial', 'fontsize',16);
% ylabel('Spatially averaged GFP', 'fontname','arial','fontsize',16);
set(gca, 'fontsize',14);
set(plotobj,...
                'LineWidth',2,...
                'MarkerEdgeColor',lincha,...
                'MarkerFaceColor',lincha,...
                'MarkerSize',8);
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5            

%curve fitting
%curve fitting
a = 55 ;
b = 1 ;
c = 1.0*10^(-3.8) ;
d = 1.85;
f = 58 - background ;

%x0=-7.2:0.1:-1.8;
x0=-10.2:0.1:-1.8;
x = 10.^(x0);

y = a - a./(b+(x./c).^d)+ f;


hold on ;
plot(x0,y,'k-','linewidth',1.0);
axis([-7.1 -1.9 0.9*min(y) 1.1*max(y)])
axis([-7.1 -1.9 0.0 1.1*max(y)])
%axis([-7.1 -1.9 0.0 82.5])

% 
% hold on ;
% errorbar(IPTG, XFP2(:,mors),XFP2(:,2),'k','linestyle','none','linewidth',1.0);
% 
% axis([-7.1 -1.9 0.9*min(Y) 1.1*max(Y)])
