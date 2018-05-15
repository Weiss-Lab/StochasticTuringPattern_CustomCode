% A program for loading and analyzing fluorescence image from Diffusion experiments.

% Challenges: (1) images may have different sizes; (2) the fluo intensity
% is not centered. 

%NOTE: 
%timeid=11-14: 12 hrs diffusion set
%timeid=21-24: 20 hrs diffusion set
%timeid=31-34: 28 hrs diffusion set
% subfold=503, plasmid responding to 3OC12HSL
% subfold=202, plasmid responding to C4HSL 


for timeid=12:1:12
    subfold = 202 ;
    
    %row=510; %col=700;
    if mod(timeid,10)==1
        close all;
    end

    %time_series =[0, 2.9, 5.5, 9.25]; % experiment from my previous studies
    time_series = [0, 150, 300, 480]; % experiment from Nicholas studies

    % Frame shift
    if (subfold==503)
        shift(11:14)=[0, 0 , 0, 30];
        shift(21:24)=[0, 0 , 0, 0];
        shift(31:34)=[0, 0 , 0, -20];    
    end

    if (subfold==202)
        shift(11:14)=[0, 3, 20, 10];
        shift(21:24)=[0, 3, 4, 3];
        shift(31:34)=[0, 3 , 3, 6];
    end

    %Record exposure time
    if subfold ==202
        %row=510;
        expo = 0.1 ; 
    end
    if subfold ==503 
        %row=480;
        expo = 0.2 ;
    end
    
    % Specify left, right and middle lines to the intensity figure.
    if subfold==202
        left1=120; % left line
        right1=1600-150; % right line
    elseif subfold==503 
        left1=120; % left line
        right1=1600-230; % right line
    else
        ;
    end
    mid1=(left1+right1)./2; % middle line
    X_length = (right1-left1) ; % length of the space 

    %%%%%%%%%%%%%%%%% Open the image files %%%%%%%%%%%%%%%%%%%%%%%%
    
    %Pathx=['.\Old_images\', num2str(subfold) '\'];
    Pathx=['.\ActualRun20151101\'];

    %pFNK202
    if (subfold == 202)
        switch timeid
            case 11
                name=['pFNK202_12hours_0min_100ms_analyzed2.tif'];
            case 12
                name=['pFNK202_12hours_150min_100ms_analyzed2.tif'];
            case 13
                name=['pFNK202_12hours_300min_100ms_analyzed2.tif'];
            case 14
                 name=['pFNK202_12hours_480min_100ms_analyzed2.tif'];
            case 21
                name=['pFNK202_20hours_480min_100ms_analyzed2.tif'];
            case 22
                name=['pFNK202_20hours_630min_100ms_analyzed2.tif'];
            case 23
                name=['pFNK202_20hours_780min_100ms_analyzed2.tif'];
            case 24
                name=['pFNK202_20hours_960min_100ms_analyzed2.tif'];
            case 31
                name=['pFNK202_28hours_960min_100ms_analyzed2.tif'];
            case 32
                name=['pFNK202_28hours_1110min_100ms_analyzed2.tif'];
            case 33
                name=['pFNK202_28hours_1260min_100ms_analyzed2.tif'];
            case 34
                name=['pFNK202_28hours_1440min_100ms_analyzed2.tif'];
            otherwise
                ;
        end
    end

    %pFNK503
    if (subfold == 503)
        switch timeid
            case 11
                name=['pFNK503_12hours_0min_200ms_analyzed2.tif'];
            case 12
                name=['pFNK503_12hours_150min_200ms_analyzed2.tif'];
            case 13
                name=['pFNK503_12hours_300min_200ms_analyzed2.tif'];
            case 14
                name=['pFNK503_12hours_480min_200ms_analyzed2.tif'];
            case 21
                name=['pFNK503_20hours_480min_200ms_analyzed2.tif'];
            case 22
                name=['pFNK503_20hours_630min_200ms_analyzed2.tif'];
            case 23
                name=['pFNK503_20hours_780min_200ms_analyzed2.tif'];
            case 24
                %name=['pFNK503_20hours_960min_200ms_analyzed2.tif'];
                %name=['pFNK503_20hours_960min_200ms_analyzed2-1.tif'];
                %name=['pFNK503_20hours_960min_200ms_analyzed2-2.tif'];
                name=['pFNK503_20hours_960min_200ms_analyzed2-3.tif'];
            case 31
                name=['pFNK503_28hours_960min_200ms_analyzed2.tif'];
            case 32
                name=['pFNK503_28hours_1110min_200ms_analyzed2.tif'];
            case 33
                name=['pFNK503_28hours_1260min_200ms_analyzed2.tif'];
            case 34
                name=['pFNK503_28hours_1440min_200ms_analyzed2.tif'];            
            otherwise
                ;
        end
    end


    % Read the image file
    filedex = [Pathx, name];
    Pic = imread(filedex);

    Pd = double(Pic) ;
    
    size(Pd)
    
    % Convert the image to gray scale
    PicGray=rgb2gray(Pic);

    figure(102);
    imshow(PicGray);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure,
    % imshow(Pic,'InitialMagnification', 'fit');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % row=480;
    % col=700;
    Nxy=size(Pd);

    color=['r', 'm', 'g', 'b', 'k', 'c', 'y'];
    %       1    2    3    4    5    6    7 

    %space= [1:size(Pd,2)]-shift(timeid);
    %space= [1:size(PicGray,2)]-shift(timeid);

    ShiA = 3*2;
    Shi_size = ShiA/3;
    col_ind = 0;
    tot_color = [];
    tot_PicGray = zeros(1, length(PicGray(1,:)));
    
    for shifter = -ShiA:Shi_size:ShiA
        col_ind = col_ind + 1;  
        space= [1:size(PicGray,2)]; % space as the column length of the image 
        % If need a normalization in space, "space" shall be scaled. 
        
        row0 = floor(size(PicGray,1)/2); % half of the row length of the image (i.e., middle of the horizontal direction of the imge)
        row_shifted = row0 + shifter; % scan around the middle horizontal line
        
        figure(10); hold on; % plot the scanned lines around the horizontal middle
        plot(space,PicGray(row_shifted,:),[color(col_ind)]);
        
%         PicGray_corrected0 = PicGray(row_shifted,:);
%         figure(33)
%         hold on;
%         PicGray_corrected = PicGray(row_shifted,:);
%         plot(space,PicGray_corrected,[color(col_ind)]);
%         %space_correct = space(25:length(space)-25);
        

        tot_color = [tot_color, sum(PicGray(row_shifted,:))./length(PicGray(row_shifted,:))];
        tot_PicGray = tot_PicGray + double(PicGray(row_shifted,:));

        figure(20);
        subplot(3,3,col_ind);
        plot(space,PicGray(row_shifted,:),[color(col_ind)]);
        Ylimit = max(PicGray(row_shifted,:));
        %axis([0 max(space) 0 200]);
    
        figure(25);
        hold on;
        plot(space,PicGray(row_shifted,:),[color(col_ind)]);
        Ylimit = max(PicGray(row_shifted,:));
        
        grid on;
    end
    %hold on; 
    box on;

    %Plot the averaged fluo intensity
    figure(30)
    hold on;
    space2= [1:size(PicGray,2)]+shift(timeid);
    meanPicGray = tot_PicGray./(2*ShiA/Shi_size +1);
    mean_color = ['r', 'm', 'g', 'b'];
    color_ind = mean_color(mod(timeid,10));
    plot(space2,meanPicGray,color_ind);
    grid on; 
    box on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Jleft = find(space2==left1);
    Jright = find(space2==right1);
    Xnorm1 = space2(Jleft:Jright);
    Ynorm1 = meanPicGray(Jleft:Jright); 
    figure (40)
    hold on; plot(Xnorm1, Ynorm1,color_ind);
    grid on; box on;
    
    figure (41)
    Xnorm2 = Xnorm1-mid1;
    hold on; plot(Xnorm2,Ynorm1,color_ind);
    grid on; box on;
    axis([(left1-mid1) (right1-mid1) 0 300]);
    
    figure (42)
    Xnorm3 = Xnorm2/(X_length/2);
    hold on; plot(Xnorm3,Ynorm1,color_ind);
    grid on; box on;
    axis([(left1-mid1)/(X_length/2) (right1-mid1)/(X_length/2) 0 300]);
    
    figure (43)
    % Find the scaling factor from pixals to dimension (cm)
    if (subfold == 202)
        plate_pixal = 1525;
    elseif (subfold == 503)
        plate_pixal = 1465;
    end
    plate_diameter = 9 ; % 9 cm plate
    Xnorm4 = Xnorm2*plate_diameter/plate_pixal;
    hold on; plot(Xnorm4,Ynorm1,color_ind);
    grid on; box on;
    %axis([(left1-mid1)/(X_length/2) (right1-mid1)/(X_length/2) 0 300]);
    axis([-4 4 0 260])

  
%     figure (44)
%     plate_pixal = 1465;
%     plate_diameter = 9 ; 
%     Xnorm4 = Xnorm2*plate_diameter/plate_pixal;
%     hold on; plot(Xnorm4,Ynorm1,color_ind);
%     grid on; box on;
%     %axis([(left1-mid1)/(X_length/2) (right1-mid1)/(X_length/2) 0 300]);
%     axis([-4 4 0 260])   %Record exposure time
%     if subfold ==202
%         expo = 0.1 ; 
%     end
%     if subfold ==503 
%         expo = 0.2 ;
%     end

    
    %size(PicGray)
    [M,I] = max(tot_color);

    title(['AHL diffusion. Exp.:  pFNK=', ...
        num2str(subfold)  ', expo=' num2str(expo) ] );

    %Save data into a file. 
    fp=fopen(['.\ActualRun20151101\' num2str(subfold), '-', num2str(timeid), '.txt'],'wt');
    for j=1:size(Ynorm1,2)
        %fprintf(fp, '%d\t%d\n',space(j),Pd(row,j));
        %fprintf(fp, '%d\t%d\n',space(j),Pd(row,j));
        fprintf(fp, '%d\t%d\n',Xnorm4(j),Ynorm1(j));
    end
    fclose(fp);
    
    size(PicGray);

end



figure(30);
x1=[left1, left1];  
y1=[0, 300];
hold on;    plot(x1, y1, 'k','linewidth', 2.0);
x1=[right1, right1];  
y1=[0, 300];
hold on;    plot(x1, y1, 'k','linewidth', 2.0);
x1=[mid1, mid1];    
y1=[0, 300];
hold on;    plot(x1, y1, 'y','linewidth', 2.0);






