function [ideal Ga_vacancy Ga_vacancy_sf Ga_vacancy_sf_TW Ga_vacancy_sf_WZ] = Eval_Ci_GaAs_111A(p_total,t,Fvib_111A_surface,Fvib_111A_surface_sf,Fvib_111A_surface_sf_TW,Fvib_111A_surface_sf_WZ,Fvib_GaAs,Fvib_As)
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
u_As_gas=(E_As2+ZPE_As2+stand_chemp_As2+kb/1.6*10^(19).*t.*log(p_As2))/2+5.3915;
u_As_solid=0.0382*t+0.002359*t-0.0004568*t.*log(t);
u_As=min(u_As_gas,u_As_solid);

%evaluation of (111)A surface energy(meV) as a function of As chemp(eV)
ideal=103.10+(103.10-91.28)/(0+0.6444).*u_As;
%As_adatom=72.17+(72.17-72.17)/(0+0.6444).*u_As;
%As_trimer_010=53.63+(53.63-77.28)/(0+0.6444).*u_As;
%As_trimer_011=48.09+(48.09-95.38)/(0+0.6444).*u_As;
%As_trimer_100=51.42+(51.42-75.07)/(0+0.6444).*u_As;
%As_trimer_100=51.42+(Fvib_111A_surface(:,1)-Fvib_GaAs*32-Fvib_As*(-4))*1000/(54.4974*4)+(51.42-75.07)/(0+0.6444).*u_As;
%As_trimer_101=46.61+(46.61-93.91)/(0+0.6444).*u_As;
%As_trimer_101=46.61+(Fvib_111A_surface(:,2)-Fvib_GaAs*28-Fvib_As*(+4))*1000/(54.4974*4)+(46.61-93.91)/(0+0.6444).*u_As;
%Ga_vacancy=52.83+(52.83-52.83)/(0+0.6444).*u_As;
Ga_vacancy=52.83+(Fvib_111A_surface-Fvib_GaAs*28-Fvib_As*(-12))*1000/(54.4974*4)+(52.83-52.83)/(0+0.6444).*u_As;
%Ga_vacancy_sf=55.98+(55.98-55.98)/(0+0.6444).*u_As;
Ga_vacancy_sf=55.98+(Fvib_111A_surface_sf-Fvib_GaAs*28-Fvib_As*(-12))*1000/(54.4974*4)+(55.98-55.98)/(0+0.6444).*u_As;
Ga_vacancy_sf_TW=54.44+(Fvib_111A_surface_sf_TW-Fvib_GaAs*28-Fvib_As*(-12))*1000/(54.4974*4)+(54.44-54.44)/(0+0.6444).*u_As;
Ga_vacancy_sf_WZ=56.86+(Fvib_111A_surface_sf_WZ-Fvib_GaAs*28-Fvib_As*(-12))*1000/(54.4974*4)+(56.86-56.86)/(0+0.6444).*u_As;

%g=[1 4 4 4 4 4 4];

%evaluation partition function Z
%Z=0;
%Z=Z+g(1)*exp(-(ideal*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(2)*exp(-(As_adatom*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(3)*exp(-(As_trimer_010*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(4)*exp(-(As_trimer_011*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(5)*exp(-(As_trimer_100*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(6)*exp(-(As_trimer_101*A)./(kb/1.6*10^(22).*t));
%Z=Z+g(7)*exp(-(Ga_vacancy*A)./(kb/1.6*10^(22).*t));

%evaluation concentraion
%C_ideal=g(1)*exp(-(ideal*A)./(kb/1.6*10^(22).*t))./Z;
%C_As_adatom=g(2)*exp(-(As_adatom*A)./(kb/1.6*10^(22).*t))./Z;
%C_As_trimer_010=g(3)*exp(-(As_trimer_010*A)./(kb/1.6*10^(22).*t))./Z;
%C_As_trimer_011=g(4)*exp(-(As_trimer_011*A)./(kb/1.6*10^(22).*t))./Z;
%C_As_trimer_100=g(5)*exp(-(As_trimer_100*A)./(kb/1.6*10^(22).*t))./Z;
%C_As_trimer_101=g(6)*exp(-(As_trimer_101*A)./(kb/1.6*10^(22).*t))./Z;
%C_Ga_vacancy=g(7)*exp(-(Ga_vacancy*A)./(kb/1.6*10^(22).*t))./Z;

%E=As_trimer_100.*C_As_trimer_100+As_trimer_101.*C_As_trimer_101+Ga_vacancy.*C_Ga_vacancy;

%E=ideal.*C_ideal+As_adatom.*C_As_adatom+As_trimer_010.*C_As_trimer_010+As_trimer_011.*C_As_trimer_011;
%E=E+As_trimer_100.*C_As_trimer_100+As_trimer_101.*C_As_trimer_101+Ga_vacancy.*C_Ga_vacancy;

%E=min(As_trimer_100,As_trimer_101);
%E=min(E,Ga_vacancy);

%E=min(ideal,As_adatom);
%E=min(E,As_trimer_010);
%E=min(E,As_trimer_011);
%E=min(E,As_trimer_100);
%E=min(As_trimer_100,As_trimer_101);
%E=min(E,Ga_vacancy);
%E=min(E,Ga_vacancy_sf);

end