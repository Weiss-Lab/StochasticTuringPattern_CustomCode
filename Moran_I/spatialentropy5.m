function out=spatialentropy5(I,options)
%#function logfactdiv
% Implementation of "Entropic measure of spatial disorder for systems of
% finite-sized objects" by R. Piasecki.  Physica A 277 (2000) 157-173
%
% inputs:  I = image matrix
%          options.maxk = ignore scales larger than this value
%          options.setk = use only k values in this vector
%
% outputs:  
%           out.Isq = binary square image processed
%           out.n = number of 'on' pixels
%           out.L = length of Isq
%           out.allS = entropic measures at each k
%           out.allSnorm = an attempt at normalizing allS to remove 
%                   dependence on percentage of "on" pixels
%                          
%
% David Karig
% August 30, 2009
% September 13, 2009 :  Support max scale option & specified kvalues


%% handle options

if exist('options')
    if isfield(options,'maxk')
        maxk=options.maxk;
        nomaxk=0;
    else
        nomaxk=1;
    end
    if isfield(options,'setk')
        setk=options.setk;
        nosetk=0;
    else
        nosetk=1;
    end
else
    nomaxk=1;
    nosetk=1;
end

%% import the image and convert to binary, square image
Ibw=im2bw(I); %make sure image is binary
lenx=size(I,2);
leny=size(I,1);
%chose largest possible square region of image
L=min(lenx,leny);
Isq=Ibw(1:L,1:L);
% figure;imshow(Isq);

%% calculate entropy for different spatial scales
if and(nomaxk,nosetk)
    allk=1:L;
elseif and(nomaxk==0,nosetk==1)
    allk=1:maxk;
elseif nosetk==0
    allk=setk;
end


allS=zeros(1,length(allk)); %store entropy for each scale
h=waitbar(0,'Microsoft style...extra slow at the end');
% figure;
for  scaleind=1:length(allk) %calculate entropy at each scale
    
    k=allk(scaleind);
    %find m = number of periodic repititions necessary (point (6) page 6)
    m=1;
    while mod(L*m,k)~=0
        m=m+1;
    end
    
    %calculate parameters
    n=sum(sum(Isq))*m^2; %number of 'on' pixels in image
    chi=(L*m/k)^2; %total number of cells
    r0=mod(n,chi);
    n0=(n-r0)/chi;
    fixedpart=(r0/chi)*log((k^2-n0)/(n0+1));
    
    %divide image into cells and calculate entropy for each
    sumsofar=0;
    for xx=1:(L*m/k)
        for yy=1:(L*m/k)
            ll=1+mod((xx-1)*k,L);
            rr=1+mod(xx*k-1,L);
            uu=1+mod((yy-1)*k,L);
            dd=1+mod(yy*k-1,L);
            
            %grab correct horizontal part
            if rr>=ll
                cellxx=Isq(:,ll:rr);
            else
                cellxx=[Isq(:,ll:L) Isq(:,1:rr)];
            end
            
            %grab correct vertical part
            if dd>=uu
                thiscell=cellxx(uu:dd,:);
            else
                thiscell=[cellxx(uu:L,:) ; cellxx(1:dd,:)];
            end
            
            ni=sum(sum(thiscell));
            %             if and(n0<100,(k^2-n0)<100)
            %                 sumsofar=sumsofar+log( (factorial(ni)*factorial(k^2-ni)) / (factorial(n0)*factorial(k^2-n0) ));
            %             else %need compact factorial division
            
            num1=ni;
            num2=(k^2-ni);
            den1=n0;
            den2=(k^2-n0);
            
            if and(num1<=num2,den1<=den2)
                logsum=logfactdiv(num1,den1)+logfactdiv(num2,den2);
            elseif and(num1<=num2,den1>den2)
                logsum=logfactdiv(num1,den2)+logfactdiv(num2,den1);
            elseif and(num1>num2,den1<=den2)
                logsum=logfactdiv(num2,den1)+logfactdiv(num1,den2);
            elseif and(num1>num2,den1>den2)
                logsum=logfactdiv(num2,den2)+logfactdiv(num1,den1);
            else
                error('oh no.')
            end
            
            sumsofar=sumsofar+logsum;
            %             end
        end
    end
    waitbar(scaleind/length(allk));
    allS(scaleind)=fixedpart+(1/chi)*sumsofar;
%     allSnorm(scaleind)=allS(scaleind)/log( nchoosek(L^2,sum(sum(Isq))) );
%     plot(1:k,allS(1:k)))
end
out.Isq=Isq;
out.n=n;
out.allk=allk;
out.L=L;
out.allS=allS;
% out.allSnorm=allSnorm;
close(h);
