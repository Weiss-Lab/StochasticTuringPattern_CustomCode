close all
clc
clear

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
allavgareas=zeros(1,numims);
allstdareas=zeros(1,numims);


for whichim=1:numims
    fname=char(input.allims(whichim));
    %     I=imresize(imread(fname),0.25);%shrink image 4x
    I=imread(fname);
    
%     Iseg=I(197:699,402:904);  %503 by 503 image
    Iseg=I(98:798,302:1002);  %701 by 701 image
    
%     Ifinal=im2bw((2^16/2^12).*Iseg,1*graythresh((2^16/2^12).*Iseg));
    meanIseg=mean(reshape(Iseg,1,numel(Iseg)));
    stdIseg=std(double(reshape(Iseg,1,numel(Iseg))));
%    Ifinal=im2bw(Iseg.*16,(meanIseg+2*stdIseg)*16/(2^16));

Ifinal=im2bw(Iseg);
    
    
%     figure;
%     subplot(1,2,1);imshow(imadjust(Iseg))
%     subplot(1,2,2);imshow(Ifinal)
    se = strel('disk',2);
    Ier = imerode(Ifinal,se);
    Ierdi=imdilate(Ier,se);
    figure;
    subplot(1,3,1);imshow(Ifinal)
    subplot(1,3,2);imshow(Ier)
    subplot(1,3,3);imshow(Ierdi)   

    cc=bwlabel(Ierdi,8);
    rr1=regionprops(cc,'Area');
    rr2=regionprops(cc,'centroid');
    centroids=cat(1,rr2.Centroid);
    areas=cat(1,rr1.Area);
    x=centroids(:,1);
    y=centroids(:,2);
    allavgareas(whichim)=mean(areas);
    allstdareas(whichim)=std(areas);
    
    for k=1:length(x)
        for j=1:length(x)
            if (k==j)
                dist(k,j)=1000000;
            else
                
            dist(k,j)=sqrt(((x(k)-x(j))^2+(y(k)-y(j))^2));
            end
        end
    end
    for k=1:length(x)
        minlength(k)=min(dist(k,:));
    end
    meanlength(whichim)=mean(minlength);
    stdlength(whichim)=std(minlength);
    maxminlength(whichim)=max(minlength);
end

figure;errorbar(1:numims,allavgareas,allstdareas./2,'o')
figure;errorbar(1:numims,meanlength,stdlength./2,'o')
figure;plot(1:numims,maxminlength,'o')

disp(['min average sqrt(area)= ' num2str(sqrt(min(allavgareas(1:10))))])
disp(['max average sqrt(area)= ' num2str(sqrt(max(allavgareas(1:10))))])
disp(['average min distance between centroids= ' num2str(max(meanlength(1:10)))])
disp(['maximum min distance between centroids= ' num2str(max(maxminlength(1:10)))])
disp(['mean min distance between centroids= ' num2str(mean(maxminlength(1:10)))])
