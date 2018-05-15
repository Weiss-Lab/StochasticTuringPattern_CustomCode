function ratio=factdiv(num,den)
%  Finds log of ratio of two factorials
%  output ratio is calculated to give ln(num!/den!)
%  
%  David Karig, August 29, 2009
partb=0;
if den>num
    for prods=1:(den-num)
        partb=partb-log(den-prods+1);
    end
elseif num>den
    for prods=1:(num-den)
        partb=partb+log(num-prods+1);
    end
elseif den==num
    partb=0;
else
    error('Whatchoo talkin bout Willis?!');
end
ratio=partb;