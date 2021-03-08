function [E] = Eval_E_GaAs_111A(p_total,t)
%evaluation of surface energy of each surface
A=54.4974;
kb=1.38*10^(-23);
E_As2=-8.09781;
ZPE_As2=0.0277;
E_As4=-19.5184;
ZPE_As4=0.0985;
stand_chemp_As2=-0.007572-6.409*10^(-5)*t-0.0003641*t.*log(t);
stand_chemp_As4=-0.04046+0.001852*t-0.0007858*t.*log(t);
a=exp((2*(E_As2+ZPE_As2+stand_chemp_As2)-(E_As4+ZPE_As4+stand_chemp_As4))./(kb*t/1.602*10^(19)));
p_As2=(-1+sqrt(1+4*a.*p_total))./(2*a);

%evlauation of As chemp as a function of PAs and T
u_As_gas=(E_As2+ZPE_As2+stand_chemp_As2+kb/1.6*10^(19).*t.*log(p_As2))/2+5.3911;
u_As_solid=0.0382*t+0.002359*t-0.0004568*t.*log(t);
u_As=min(u_As_gas,u_As_solid);

%evaluation of (111)A surface energy(meV/A^2) as a function of As chemp(eV)
As_trimer_100=51.31+(51.31-74.95)/(0+0.6441).*u_As;
As_trimer_101=46.50+(46.50-93.77)/(0+0.6441).*u_As;
Ga_vacancy=52.76+(52.76-52.76)/(0+0.6441).*u_As;

E=min(As_trimer_100,As_trimer_101);
E=min(E,Ga_vacancy);

end