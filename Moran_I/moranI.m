function moran=moranI(im) 
%
% Returns MoranI statistic for an image im
%
% David Karig, September 10, 2009
%

% %%%%for debugging%%%%%%%%
% clc
% im=rand(20,20);
% %%%%%%%%%%%%%%%%%%%%%%%%%

lr=size(im,1);
lc=size(im,2);

%reshape image array into col vector
y=reshape(im',numel(im),1);
y=double(y-mean(y)); %mean suppress

N=numel(im); %sample size

% if N<250 %matlab sucks at big things 
%     % weight matrix for distance=1
%     w=(makedistmap(lr,lc)==1);
%     S0=sum(sum(w));%sum of weights
%     moran=N/S0*(y'*w*y)/(y'*y);
% else %need to use the stupid way 
    
    imms=double(im-sum(sum(im))/N); %mean suppress
    sumit=0;
    S0=0;
    for ii=1:lr
        for jj=1:lc
            if (ii+1)<=lr
                sumit=sumit+imms(ii,jj)*imms(ii+1,jj);
                S0=S0+1;
            end
            if (ii-1)>0
                sumit=sumit+imms(ii,jj)*imms(ii-1,jj);
                S0=S0+1;
            end
            if (jj+1)<=lc
                sumit=sumit+imms(ii,jj)*imms(ii,jj+1);
                S0=S0+1;
            end
            if (jj-1)>0
                sumit=sumit+imms(ii,jj)*imms(ii,jj-1);
                S0=S0+1;
            end
        end
%         if N>10000
%             waitbar(ii/lr);
%         end 
    end
    %divide by variance and normalize by sample size and sum of weights 
    moran=sumit*N/(S0*y'*y);
    
    
    
%     %really slow implementation
%     r=repmat((1:lr)',1,lc);
%     r=reshape(r',numel(r),1);
%     c=repmat((1:lc)',lr,1);
%     c=reshape(c,numel(c),1);
%     sumit=0;
%     S0=0;
%     for ii=1:N
%         for jj=1:N
%             d=sqrt((r(ii)-r(jj))^2+(c(ii)-c(jj))^2);%distance
%             if d==1
%                 w=1;
%             else
%                 w=0;
%             end
%             sumit=sumit+w*y(ii)*y(jj);
%             S0=S0+w;
%         end
%     end
%     %divide by variance and normalize by sample size and sum of weights
%     moran=sumit*N/(S0*y'*y)

end 
