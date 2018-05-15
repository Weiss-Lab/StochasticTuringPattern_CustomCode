% This code is used to analyze the statistics of the spots from Turing
% pattern images. 

% Input: images from experiments
% Output: Statistics of the spots, including the max, min, mean, and the total
% number of the spots 
% Filter is applied to discard tiny spots below a threshold 


clear all; 
close all;

%Key parameters for the code:
para_bw_threshold = 0.13; % threshold for imaging conversion (im2bw) to black-and-white. The threshold above is very critical.
para_mini_ditected_area = 100; %minimal spot sizes to be detected (areas smaller than that will be filtered out)



%read the image
    FigX2 = imread('FIG2v3','bmp');

%show the image and convert it to black-and-white with right threshold
    %figure(1),
    %subplot(1,2,1),
    figure(1),
    imshow(FigX2);  

    Fig_bw = im2bw(FigX2,para_bw_threshold);

    %figure(1),
    %subplot(1,2,2);
    %figure(2),
    imshow(Fig_bw);

%show the domains using pseudocolor image
    PseudoFig = label2rgb(Fig_bw, @spring, 'c', 'shuffle'); 
    %figure(2),
    %subplot(1,2,1),
    figure(3)
    imshow(PseudoFig,'InitialMagnification', 'fit');  

%Analyze domain characters
    %Count and label the domains
    [labeledImage,numObjects] = bwlabel(Fig_bw,8);
    %numObjects; % Count all distinct objects in the image.
    %disp(['original numObjects= ' num2str(numObjects)]);

    domaindata = regionprops(labeledImage,'basic'); %which returns: area, centroid, and 'bounding box'

    alldomains = [domaindata.Area];
    disp(['Original']);
    num_area=length(alldomains);
    disp(['num_area(ori)= ' num2str(num_area)]);
    mean_area=mean(alldomains);  % Find the mean area size for the filtered image.
    disp(['mean_area(ori)= ' num2str(mean_area)]);
    max_area = max(alldomains); % Find the maximum area of all the grains.
    disp(['max_area(ori)= ' num2str(max_area)]);
    min_area = min(alldomains); % Find the minimum area of all the grains.
    disp(['max_area(ori)= ' num2str(min_area)]);
   
    %Collect all of the actual (non-zero)areas:
    alldomainsR=alldomains(alldomains~=0); %get all the non zero values
    disp(['Original Real']);
    num_area=length(alldomains);
    disp(['num_area(oriR)= ' num2str(num_area)]);
    mean_area=mean(alldomains);  % Find the mean area size for the filtered image.
    disp(['mean_area(oriR)= ' num2str(mean_area)]);
    max_area = max(alldomains); % Find the maximum area of all the grains.
    disp(['max_area(oriR)= ' num2str(max_area)]);
    min_area = min(alldomains); % Find the minimum area of all the grains.
    disp(['max_area(oriR)= ' num2str(min_area)]);
    
    %filter small domains
    thres_area = para_mini_ditected_area;
    idx = find([domaindata.Area] > thres_area );
    Judge_domain = ismember(labeledImage, idx);

    %show the filtered image
    Re_labeledImg = labeledImage .* Judge_domain ;
    
    %figure(2), subplot(1,2,2), 
    figure(4),
    imshow(Re_labeledImg,'InitialMagnification', 'fit');

    %Analyze the filtered image
    domaindataF = regionprops(Re_labeledImg,'basic'); %which returns: area, centroid, and 'bounding box'

    alldomainsF = [domaindataF.Area];
    
    disp(['Filtered']);
    num_areaF=length(alldomainsF);  % Find area number for the filtered image.
    disp(['num_area(fil)= ' num2str(num_areaF)]);
    mean_areaF=mean(alldomainsF);  % Find the mean area size for the filtered image.
    disp(['mean_area(fil)= ' num2str(mean_areaF)]);
    max_areaF = max(alldomainsF);  % Find the maximal area in the filtered image.
    disp(['max_area(fil)= ' num2str(max_areaF)]);
    min_areaF = min(alldomainsF);  % Find the minimal area in the filtered image.
    disp(['max_area(fil)= ' num2str(min_areaF)]);
    
    %Collect all of the actual (non-zero)areas:
    alldomainsFR=alldomainsF(alldomainsF~=0); %get all the non zero values

    disp(['Filtered Real']);
    num_areaFR=length(alldomainsFR);  % Find area number for the filtered image.
    disp(['num_area(filR)= ' num2str(num_areaFR)]);
    mean_areaFR=mean(alldomainsFR);  % Find the mean area size for the filtered image.
    disp(['mean_area(filR)= ' num2str(mean_areaFR)]);
    max_areaFR = max(alldomainsFR);  % Find the maximal area in the filtered image.
    disp(['max_area(filR)= ' num2str(max_areaFR)]);
    min_areaFR = min(alldomainsFR);  % Find the minimal area in the filtered image.
    disp(['min_area(filR)= ' num2str(min_areaFR)]);

    %Nomalize the image according to the scale bar
    %image size: 1206*896 pixels; 7.46*Scale_bar = 1206 pixels; 
    %Scale_bar = 1206/7.46=161.6622 pixels
    % Scale_bar=100 um
    um_per_pix = 100*7.46/1206; 
    
    disp(['Filtered Real Scaled']);
    alldomainsFRS = alldomainsFR * um_per_pix * um_per_pix;
    num_areaFRS=length(alldomainsFRS);  % Find area number for the filtered image.
    disp(['num_area(filRS)= ' num2str(num_areaFRS)]);
    mean_areaFRS=mean(alldomainsFRS);  % Find the mean area size for the filtered image.
    disp(['mean_area(filRS)= ' num2str(mean_areaFRS)]);
    max_areaFRS = max(alldomainsFRS);  % Find the maximal area in the filtered image.
    disp(['max_area(filRS)= ' num2str(max_areaFRS)]);
    min_areaFRS = min(alldomainsFRS);  % Find the minimal area in the filtered image.
    disp(['min_area(filRS)= ' num2str(min_areaFRS)]);

    
