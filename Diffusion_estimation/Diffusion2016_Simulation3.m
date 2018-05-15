%This code is used to simulate fluo evolution over time for a given initial AHL amount.
%Mathematical model: 
% Variables: S: AHL; F: Fluorescence; N: cell density  
% dS/dt = Diff_s * d^2/dx^2 (S) - gamma_s * S ;
% dF/dt = afa_f * N * S/(K_s + S) - gmma_f * F ;
% dN/dt = afa_n * N * (1 - N/N_l);

clear all;
close all;

% There are two ways to normalize:
% Divide everything by the mean fluo at T=0;
% Minus everything by the mean fluo at T=0;

subfold = 503 ; %202 or 503;

timeidi = 11;
timeidf = timeidi + 3;


%Load pamameters;
%===============================================
%========== Original parameters ====================
% gamma_f = 0.04;
% gamma_n = 0.00;
% beta_f = 0.0;
% K_s = 3.0e0; %uM
% N_l = 5.0;

%========== new parameters for 503 ====================
if (subfold == 503)
    switch timeidi
        case 11 % (subfold == 503)
            afa_n = 0.15; %Ori: 0.15
            N_l = 1.3; %Ori: 5.0
            gamma_f = 0.04; %Ori: 0.04
            gamma_n = 0.00; %Ori: 0.0
            beta_f = 0.0; %Ori: 0.0
            K_s = 3.0e0; %uM %Ori: 3.0
            afa_f=0.18; %0.2; %0.15; %Ori: 
            Diff_s = 0.01.*1/3; %Ori: 0.003
            gamma_s = 0.02; %Ori: 
        case 21 % (subfold == 503)
            afa_n = 0.15; %Ori: 0.15
            N_l = 1.3; %Ori: 5.0
            gamma_f = 0.04; %Ori: 0.04
            gamma_n = 0.00; %Ori: 0.0
            beta_f = 0.0; %Ori: 0.0
            K_s = 3.0e0; %uM %Ori: 3.0
            afa_f=0.12; %0.15; %Ori: 
            Diff_s = 0.008.*1/3; %Ori: 0.003
            gamma_s = 0.02; %Ori: 
  
        case 31 % (subfold == 503)
            afa_n = 0.05; %Ori: 0.15
            N_l = 1.2; %Ori: 5.0
            gamma_f = 0.04; %Ori: 0.04
            gamma_n = 0.00; %Ori: 0.0
            beta_f = 0.0; %Ori: 0.0
            K_s = 3.0e0; %uM %Ori: 3.0
            afa_f=0.07; %0.15; %Ori: 
            Diff_s = 0.005.*1/3; %Ori: 0.003
            gamma_s = 0.02; %Ori: 
 
        otherwise
            ;
    end

%     afa_n = 0.15; %Ori: 0.15
%     N_l = 1.3; %Ori: 5.0
%     
%     gamma_f = 0.04; %Ori: 0.04
%     gamma_n = 0.00; %Ori: 0.0
%     beta_f = 0.0; %Ori: 0.0
%     K_s = 3.0e0; %uM %Ori: 3.0
%     afa_f=0.2; %0.15; %Ori: 
%     Diff_s = 0.01; %Ori: 0.003
%     gamma_s = 0.02; %Ori: 
elseif (subfold == 202)
%========== new parameters for 202 ====================
    switch timeidi
        case 11  % (subfold == 202)
            afa_n = 0.15; %Ori: 
            N_l = 1.3; %Ori: 
            gamma_f = 0.04; %Ori: 
            gamma_n = 0.00; %Ori: 
            beta_f = 0.0; %Ori: 
            K_s = 3.0e0; %Ori: 
            afa_f=1.15; %Ori: 
            Diff_s = 0.25.*1/3; %0.065; %Ori: 0.065
            gamma_s = 0.2; %Ori: 
        case 21 % (subfold == 202)
            afa_n = 0.15; %Ori: 
            N_l = 1.3; %Ori: 
            gamma_f = 0.04; %Ori: 
            gamma_n = 0.00; %Ori: 
            beta_f = 0.0; %Ori: 
            K_s = 3.0e0; %Ori: 
            afa_f=1.2; %Ori: 
            Diff_s = 0.25.*1/3; %0.065; %Ori: 0.065
            gamma_s = 0.2; %Ori: 
        case 31 % (subfold == 202)
            afa_n = 0.05; %Ori: 
            N_l = 1.2; %Ori: 
            gamma_f = 0.04; %Ori: 
            gamma_n = 0.00; %Ori: 
            beta_f = 0.0; %Ori: 
            K_s = 3.0e0; %Ori: 
            afa_f=0.40; %Ori: 
            Diff_s = 0.16.*1/3; %0.065; %Ori: 0.065
            gamma_s = 0.2; %Ori: 
        otherwise
            ;
    end
end

% 60 mm plates instead of 90 mm plates
%Diff_s = Diff_s .*1/3;

if subfold ==503
%     %========== Original parameters ====================
%     afa_f=0.12;
%     Diff_s = 0.003; %0.4;
%     gamma_s = 0.02; %0.2;
%     afa_n = 0.075;
    %========== New parameters ===========================
end
if subfold ==202
%     %========== Original parameters ====================
%     afa_f=0.9;
%     %Diff_s = 0.065; %0.4;
%     Diff_s = 0.065; %0.4;
%     gamma_s = 0.2; %0.2;
%     afa_n = 0.11;
    %========== New parameters ====================
end

if subfold ==503
    switch timeidi
        case 11
            mean_fluo0 = 95; 
        case 21 
            mean_fluo0 = 117;
            %mean_fluo0 = 95; 
        case 31
            mean_fluo0 = 145;
            %mean_fluo0 = 95; 
        otherwise
            ;
    end
elseif subfold ==202
    switch timeidi
        case 11
            mean_fluo0 = 31.94;
        case 21 
            mean_fluo0 = 30;
            %mean_fluo0 = 31.94;
        case 31
            %mean_fluo0 = 35;
            %mean_fluo0 = 31.94;
            mean_fluo0 = 70;
        otherwise
            ;
    end
end

%===============================================

%Load initial conditions
%===============================================
U_0 = 10.0e0; %mM

if subfold ==503
%     %========== Original parameters ====================
%    x0=0.3; %x0>=0
     %========== New parameters ====================
    %x0=1.2; 
    switch timeidi
        case 11
            x0=1.0; %1.2;  %11
        case 21
            x0=0.9;  %21
        case 31
            x0=0.83;
        otherwise
            ;
    end
end
if subfold ==202
%     %========== Original parameters ====================
%    x0=0.25;
     %========== New parameters ====================
     %x0=0.3;
    switch timeidi
        case 11
            x0=0.3;  %11
        case 21
            x0=1.0; %1.2;  %21
        case 31
            x0=0.7;
        otherwise
            ;
    end
end

% 60mm plate instead of 90 mm plates
x0 = x0 .*2/3;

U_ini = U_0./(x0.*2) ;
%===============================================

%Boundary conditions
%===============================================
%===============================================

%Save parameter sets to a file 
Para = [gamma_f, gamma_n, afa_f, afa_n, Diff_s, gamma_s, K_s, N_l, ...
    U_ini, x0, beta_f]; 
%Para = [(1)gamma_f, (2)gamma_n, (3)afa_f, (4)afa_n, (5)Diff_s, (6)gamma_s, 
%(7)K_s, (8)N_l, (9)U_ini, (10)x0, (11)beta_f]; 

fid = fopen('Para.txt', 'wt');
fprintf(fid, '%f\n', Para);
fclose(fid);

xl = -4.0 ;
xr = -xl ;
t0 = 0 ;
tf = 10 ;
%Nt_step = 100 ;
Nt_step = 1000 ;

%x = linspace(xl,xr,61);
x = linspace(xl,xr,61);
t = linspace(t0,tf,Nt_step+1);

m = 0;
%sol = pdepe(m,@pdex1pde3,@pdex1ic3,@pdex1bc3,x,t);
sol = pdepe(m,@pdex1pde3,@pdex1ic3,@pdex1bc3,x,t);

    S = sol(:,:,1);
    F = sol(:,:,2);
    N = sol(:,:,3);

%Total Fluo = Real Fluo + AutoFluorescence    
    %F = F + N.*1.0;
    %sol(:,:,2) = sol(:,:,2) + sol(:,:,3).*1.0 ;
    F_tot = F + sol(:,:,3).*1.0 ;

    
%Plot the simulated results
    figure(100)
    subplot(2,2,1)
    h1=surf(x,t,S);
    title(['Numerical solution computed: S (AHL)']);
    xlabel('Distance x');
    ylabel('Time t');
    set(h1,'linestyle', 'none');
    %colorbar
    view(2) % view in 2D
    
    subplot(2,2,2)
    h1=surf(x,t,F_tot);
    title(['Numerical solution computed: F (Fluores)']);
    xlabel('Distance x');
    ylabel('Time t');
    set(h1,'linestyle', 'none');
    view(2)
    
    subplot(2,2,3)
    h1=surf(x,t,N);
    title(['Numerical solution computed: N (Cell Den)']);
    xlabel('Distance x');
    ylabel('Time t');
    set(h1,'linestyle', 'none');
    view(2)
    
    %Plot the simulated time evolution of fluorescence intensity
    figure(200)
    %figure,
    %subplot(2,2,4)
    %figure,
    Nk = 5 ;
    Nsize = floor(Nt_step/Nk) ;
    %F.P.

    color=['r', 'm', 'g', 'b', 'k', 'ko', 'kx'];

    %Time points where images were taken      
    %X_val=[0; 29; 55; 92]+1; % Time points for the 2009 experiments
    
 %! it seems that for strain 202, X_val=[0; 15; 50; 80]+1 is better than X_val=[0; 25; 50; 80]+1;  
 % while for Strain 503, X_val=[0; 25; 50; 80]+1 is good. 
 X_val=[0; 25; 50; 80]+1; 
    
 if subfold == 503   
%     X_val=[0; 25; 50; 80]+1; 
     X_val=[0; 250; 500; 800]+1; 
 elseif subfold == 202
%    X_val=[0; 15; 50; 80]+1; 
     X_val=[0; 150; 500; 800]+1; 
end;
   
%X_val=[0; 15; 50; 80]+1; 

    %X_val=[0; 25; 50; 80]+1; % Time points for the 2015 experiments
    % one time step = 0.1 hr.
    
    for k=1:4
        plot(x,F_tot(X_val(k),:),[color(k) '-']);
        hold on;
    end
    title(['Fluoresence in different time: ' ...
        'D_s=' num2str(Diff_s) ', \gamma_s=' num2str(gamma_s)  ',K_s=' num2str(K_s) ] );
    xlabel('Distance x')
    ylabel('Fluoscence (a.u.)')
    axis([-4 4 -0.2 max(max(F_tot(:,:))) ]);
    box on;
    
    %Plot the comparison between experiments and simulations
    figure('Position', [200, 300, 1500, 600]);

        %experimental data
        subplot(1,2,1)
        %    timeidi = 11;
        %    timeidf = 14;
            for timeid = timeidi:1:timeidf
                %Read & plot the experimental data. 
                formatSpec = '%f %f\n';
                sizeA = [2 Inf];
                
                fp=fopen(['.\ActualRun20151101\' num2str(subfold), '-', num2str(timeid), '.txt'],'r');
                    A = fscanf(fp,formatSpec,sizeA);
                fclose(fp);
               
                A = A';
               
                mean_color = ['r', 'm', 'g', 'b'];
                color_ind = mean_color(mod(timeid,10));
                title(['Original Experimental Fluorescence over Time: ', '[Strain-', num2str(subfold), ']']);                
                xlabel('Distance X (cm)');
                ylabel('Fluoscence (a.u.)');
                plot(A(:,1),A(:,2),color_ind); 
                A=[];
                grid on; hold on;
                
  %               if (timeid == timeidi)
  %                   mean_fluo0 = mean(A(:,2));
  %               end
                
            end
            axis([-4 4 0 300 ]);
            box on;  

        %simulated results
        subplot(1,2,2)
            for k=1:4
                plot(x,F_tot(X_val(k),:),[color(k) '-']);
                hold on;
            end
            title(['Fluoresence in different time: ' ...
                'D_s=' num2str(Diff_s) ', \gamma_s=' num2str(gamma_s)  ',K_s=' num2str(K_s) ] );
            xlabel('Distance X (cm)')
            ylabel('Fluoscence (a.u.)')
            axis([-4 4 -0.2 max(max(F_tot(:,:))) ]);
            box on; grid on;   
        
%        fig = gcf;
        %FigHandle = figure('Position', [200, 300, 1500, 600]);

 
   %Plot the comparison between experiments and simulations
    %figure('Position', [200, 300, 1500, 600]);
    figure('Position', [200, 300, 750, 300]);

        %experimental data
        subplot(1,2,1)
            for timeid = timeidi:1:timeidf
                %Read & plot the experimental data. 
                formatSpec = '%f %f\n';
                sizeA = [2 Inf];
                
                fp=fopen(['.\ActualRun20151101\' num2str(subfold), '-', num2str(timeid), '.txt'],'r');
                    A = fscanf(fp,formatSpec,sizeA);
                fclose(fp);
               
                A = A';
               
                mean_color = ['r', 'm', 'g', 'b'];
                color_ind = mean_color(mod(timeid,10));
                
                plot(A(:,1).*2/3,A(:,2)./mean_fluo0,color_ind, 'linewidth', 1.0);
                
                if subfold == 503
                    maxY = 2.6;
                elseif subfold == 202
                    maxY = 8.0;
                end
                %axis([-4 4 0 maxY]);
                axis([-3 3 0 maxY]);
                %axis([-4 4 0 max(max(F_tot(:,:)))]);
                %title(['Experimental Fluorescence over Time: ', '[Strain-', num2str(subfold), ']']);                
                %title(['Experimental Fluorescence over Time: ', '[Strain-', ... 
                %    num2str(subfold), '; Time=', num2str(timeidi),'-', num2str(timeidf) ']']);                
                title(['Experiment']);                

                xlabel('Distance X (cm)')
                ylabel('Normalized Fluoscence (a.u.)')               
                %grid on; 
                hold on;
                
            end
      



        %simulated results
        subplot(1,2,2)
            for k=1:4
                plot(x,F_tot(X_val(k),:),[color(k) '-'],'linewidth', 1.5);
                hold on;
            end
            %title(['Simulated Fluoresence over Time ', '[Strain-', num2str(subfold), '], ' ...
            %    'D_s=' num2str(Diff_s) ', \gamma_s=' num2str(gamma_s)  ',K_s=' num2str(K_s) ] );
            title(['Simulation']);                
            xlabel('Distance x')
            ylabel('Fluoscence (a.u.)')
            axis([-3 3 0 maxY]);
            %axis([-4 4 0 maxY]);
            %axis([-4 4 0 max(max(F_tot(:,:))) ]);
            box on; 
            %grid on;   
        
   %Plot the comparison between experiments and simulations
    figure('Position', [600, 300, 700, 600]);
           for timeid = timeidi:1:timeidf
                %Read & plot the experimental data. 
                formatSpec = '%f %f\n';
                sizeA = [2 Inf];
                
                fp=fopen(['.\ActualRun20151101\' num2str(subfold), '-', num2str(timeid), '.txt'],'r');
                    A = fscanf(fp,formatSpec,sizeA);
                fclose(fp);
               
                A = A';
               
                mean_color = ['r', 'm', 'g', 'b'];
                color_ind = mean_color(mod(timeid,10));
                
                %experimental data
 %               plot(A(:,1),A(:,2)./mean_fluo0,color_ind);
 % 60 mm plates instead of 90 mm.
                plot(A(:,1).*2/3,A(:,2)./mean_fluo0,color_ind);

                if subfold == 503
                    maxY = 2.6;
                elseif subfold == 202
                    maxY = 8.0;
                end
                axis([-4 4 0 maxY]);
                %axis([-4 4 0 max(max(F_tot(:,:)))]);
                title(['Experimental Fluorescence over Time: ', '[Strain-', ... 
                    num2str(subfold), '; Time=', num2str(timeidi),'-', num2str(timeidf) ']']);                
                xlabel('Distance X (cm)')
                ylabel('Normalized Fluoscence (a.u.)')               
                grid on; hold on;
     
           end

           %simulated data
           for k=1:4
                hold on; plot(x,F_tot(X_val(k),:),[color(k) '-']);
           end
