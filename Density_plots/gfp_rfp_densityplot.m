% a code to show gfp-rfp density plot (log scale)
clear all;
close all;

tic ;
disp([' ']);
disp(['new processing ...']);
disp(['']);

itype = 3 ; %1-tif, 2-jpg %speical tif

fpick = 5 ;

togchan = 3 ;

%file folder
if itype == 1
    %dirname=['D:\Turings data\EXPO\Control_exp\']; %my laptop
    dirname=['C:\Research\Turings data\EXPO\Control_exp\']; % littleorange
elseif itype == 2
    %dirname=['D:\Turings data\EXPO\Data4pub\Data4pub2\']; %my laptop
    dirname=['C:\Research\Turings data\EXPO\Data4pub\Data4pub2\']; %littleorange
else
    %dirname=['C:\Research\Turings data\EXPO\Exported Images\10X_cell_layer\np_col-4\']; %littleorange
    dirname=['C:\Users\Ting\Pictures\Exported Images\10X_cell_layer3\np_col-', num2str(togchan), '\linear_gray16\'];
end

%files
switch fpick
    case 1
    filename = ['1OD--Mix-pASK201_and_pINV5-3e-3 M IPTG-25hrs-10X-MG1655-expo150-200_v1_'];
    case 2
    filename = ['1OD--pASK201-3e-3 M IPTG-25hrs-10X-MG1655-expo150-200_v1_'];
    case 3
    filename = ['1OD--pINV5-3e-3 M IPTG-25hrs-10X-MG1655-expo150-200_v5_'];
    case 4
    filename = ['2OD-1e-2 M IPTG-25hrs-10X-pFNK805.10-pFNK512-MG1655-expo150-2k_v4_'];
    
    case 5
    %filename = ['np_col-4_c45_gray_'];
    filename = ['np_col-', num2str(togchan),'_'];
    
    
    otherwise
        disp('wrong file-choosing\n');
end

file = [dirname, filename];

if (itype ==1)
    img = imread([file 'c2.TIF']) ;
    imr = imread([file 'c3.TIF']) ;
else if itype==3
        if togchan == 4
            img = imread([file 'c5.TIF']) ;
            imr = imread([file 'c6.TIF']) ;
        end
        if togchan == 3
            img = imread([file 'c8.TIF']) ;
            imr = imread([file 'c9.TIF']) ;
        end
    end
    
    gfp = double(img) ;
    rfp = double(imr) ;
end


if itype == 2
    if fpick ==4
        img = imread([file 'c2-2.JPG']) ;
        imr = imread([file 'c3-3.JPG']) ;
    else
        img = imread([file 'c2.JPG']) ;
        imr = imread([file 'c3.JPG']) ;
    end
    
    %for jpg
    imgray_g = rgb2gray(img);
    imgray_r = rgb2gray(imr);
    gfp = double(imgray_g) ;
    rfp = double(imgray_r) ;
end

gfp1 = reshape(gfp,size(gfp,1)*size(gfp,2),[]);
rfp1 = reshape(rfp,size(rfp,1)*size(rfp,2),[]);

m_gfp = mean(gfp1);
m_rfp = mean(rfp1);

gfp1n = gfp1./m_gfp ;
rfp1n = rfp1./m_rfp ;

if fpick==4
    if itype == 2
        rfp1n = rfp1n./6 ;
    else
        rfp1n = rfp1n./6 ;
    end
end


FP = [gfp1n,rfp1n];

%FP =[gfp1, rfp1] ;


%generate linearly spaced edges:
%Xn = 50 ;
%Yn = 50 ;

%Xn = 32 ;
%Yn = 32 ;

Xn = 20 ;
Yn = 20 ;
%xedges = linspace(0,1.6,Xn)'; 
%yedges = linspace(0,0.3,Yn)';

xedges = linspace(0, 1.6, Xn)'; 
yedges = linspace(0, 1.6, Yn)';

%logspace
n = length(FP) ;
%n = 5e4;
H = zeros(Yn,Xn) ;
for i = 1:n
    x = dsearchn(xedges,FP(i,1)) ;
    y = dsearchn(yedges,FP(i,2)) ;
    H(y,x) = H(y,x) + 1 ;
end ;

%set '0' as NaN
% for i=1:size(H,1)
%     for j=1:size(H,2)
%         if H(i,j)==0
%             H(i,j)=NaN;
%         end
%     end
% end

figure , h=pcolor(xedges,yedges,H) ;
set(h,'Linestyle','none');
% 
% tic ;
% hist2d(FP);
% %hist2d(D,Xn,Yn,[Xlo Xhi],[Ylo Yhi])
% hist2d(FP,32,32,[0 2.0],[0 2.0])
% toc ;

figure, h2=pcolor(xedges,yedges,log10(H)) ;
set(h2,'Linestyle','none');
colorbar;
set(gca,'Clim', [0 4.7]);

toc ;