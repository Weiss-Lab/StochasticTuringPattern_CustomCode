function F = modelFuncT(g1,g2)

the1=1;
the2=2;
the3=1;
the4=4;
the5=2;
the6=1;

D=g2(1);
D2=g2(2);

gu =g2(3);
gv =g2(4);
gc =g2(5);
giu=g2(6);
giv=g2(7);

au =g2(8);
av =g2(9);

aiu=g2(10);
aiv=g2(11);
ac =g2(12);


lam_u=g2(13);
lam_v=g2(14);

fod_1 =g2(15);
fod_2 =g2(16);
fod_3 =g2(17);
fod_4 =g2(18);
fod_5 =g2(19);
fod_6 =g2(20);

Kd_1 = g2(21);
Kd_2 = g2(22);
Kd_3 = g2(23);

Kd_5 = g2(24);

Kc_3 = g2(25);



Kd_4 = g2(26);
Kd_6 = g2(27);


Lac  = g2(28);
IPTG = g2(29);



U  = g1(1);
V  = g1(2);
Iu = g1(3);
Iv = g1(4);
C  = g1(5);


%prefactors
Ru = lam_u * Iu ;
Rv = lam_v * ( 1+1/fod_5*((C/Kd_5)^the5) )/( 1+((C/Kd_5)^the5) ) ; %Model 804

Leff= Lac*( 1+1/fod_6*((IPTG/Kd_6)^the6) )/( 1+((IPTG/Kd_6)^the6) ); %effective Lac

X1 = Ru*U ;			
X2 = Rv*V/( 1+U/Kc_3 ) ;

HFun_uv = ( 1+fod_1*((X1/Kd_1)^the1) )/( 1+((X1/Kd_1)^the1) )*( 1+1.0/fod_2*((C/Kd_2)^the2) )/( 1+((C/Kd_2)^the2) );

HFun_c  = ( 1+fod_3*((X2/Kd_3)^the3) )/( 1+((X2/Kd_3)^the3) )*( 1+1.0/fod_4*((Leff/Kd_4)^the4) )/( 1+((Leff/Kd_4)^the4) );

F=[(au*Iu-gu*U);
(av*Iv-gv*V);
(aiu*HFun_uv-giu*Iu);
(aiv*HFun_uv-giv*Iv);
(ac*HFun_c-gc*C)];
end