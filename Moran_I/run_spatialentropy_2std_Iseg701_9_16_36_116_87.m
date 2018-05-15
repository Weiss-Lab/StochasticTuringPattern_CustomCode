function run_spatialentropy
%#function spatialentropy5
%#function logfactdiv

% close all
% clc
% clear



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

input.savestring={
    'processed_entropy__1e-7_0'; ...
    'processed_entropy__1e-6_0'; ...
    'processed_entropy__1e-5_5'; ...
    'processed_entropy__1e-5_0'; ...
    'processed_entropy__1e-4_5'; ...
    'processed_entropy__1e-4_0'; ...
    'processed_entropy__1e-3_5'; ...
    'processed_entropy__1e-3_0'; ...
    'processed_entropy__1e-2_5'; ...
    'processed_entropy__1e-2_0'; ...
%     'processed_entropy__Mix-red';...
%     'processed_entropy__pINV5';...
%     'processed_entropy__pASK201';...
    };

numims=size(input.allims,1)

IPTGconcs=[10^-7 10^-6 10^-5.5 10^-5 10^-4.5 10^-4 10^-3.5 10^-3 10^-2.5 10^-2];

for whichim=1:numims
    fname=char(input.allims(whichim));
    %     I=imresize(imread(fname),0.25);%shrink image 4x
    I=imread(fname);
    %     Iseg=I(197:699,402:904);
    Iseg=I(98:798,302:1002);  %701 by 701 image
    %     Ifinal=im2bw((2^16/2^12).*Iseg,1*graythresh((2^16/2^12).*Iseg));
    meanIseg=mean(reshape(Iseg,1,numel(Iseg)));
    stdIseg=std(double(reshape(Iseg,1,numel(Iseg))));
%    Ifinal=im2bw(Iseg.*16,(meanIseg+2*stdIseg)*16/(2^16)); %!
    
Ifinal=im2bw(Iseg); %!    
    
    %     figure;imshow(Ifinal)
    %     options.maxk=70;
        options.setk=[9 16 36 87 116];  %!
    out=spatialentropy5(Ifinal,options);
    save(char(input.savestring(whichim)));
end
