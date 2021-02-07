function [trans, vib, zpe, rot] = Eval_chemp(mass,coordination,freq,sym,temp)
%mass(kg), coordination(m), freq(Hz), temp(K)
h=6.626*10^(-34);
haba=h/2/pi;
kb=1.38*10^(-23);
%evaluate total mass of ideal gas
m_total=0;
a=length(mass);
for i=1:a
    m_total = m_total+mass(i);
end
%evaluate principal moments of inertia
x_cm=0;
y_cm=0;
z_cm=0;
for i=1:a
    x_cm = x_cm + mass(i)*coordination(i,1)/m_total;
    y_cm = y_cm + mass(i)*coordination(i,2)/m_total;
    z_cm = z_cm + mass(i)*coordination(i,3)/m_total;
end
Inertia=zeros(3,3);
for i=1:a
    Inertia(1,1) = Inertia(1,1) + mass(i)*((coordination(i,2)-y_cm).^2+(coordination(i,3)-z_cm).^2);
    Inertia(2,2) = Inertia(2,2) + mass(i)*((coordination(i,1)-x_cm).^2+(coordination(i,3)-z_cm).^2);
    Inertia(3,3) = Inertia(3,3) + mass(i)*((coordination(i,1)-x_cm).^2+(coordination(i,2)-y_cm).^2);
    Inertia(1,2) = Inertia(1,2) + mass(i)*(coordination(i,1)-x_cm).*(coordination(i,2)-y_cm);
    Inertia(1,3) = Inertia(1,3) + mass(i)*(coordination(i,1)-x_cm).*(coordination(i,3)-z_cm);
    Inertia(2,3) = Inertia(2,3) + mass(i)*(coordination(i,2)-y_cm).*(coordination(i,3)-z_cm);
end
Inertia(2,1)=Inertia(1,2);
Inertia(3,1)=Inertia(1,3);
Inertia(3,2)=Inertia(2,3);
Principal_Inertia=eig(Inertia);

%evaluate chemical potential of translational motion of ideal gas
trans = -kb/1.6*10^19*temp.*log((2*pi*m_total/h^2)^(1.5).*(kb.*temp).^(2.5)/101325);

%evaluate chemical potential of vibrational motion of ideal gas
vib=0;
zpe=0;
a=length(freq);
for i=1:a
    vib = vib+kb/1.6*10^(19)*temp.*log((1-exp(-h*freq(i)/kb./temp)));
    zpe = zpe+1/2*h*freq(i)/1.6*10^(19);
end

%evaluate chemical potential of rotational motion of ideal gas
rot = -kb/1.6*10^(19)*temp.*log(pi^(0.5)/sym.*sqrt(2*Principal_Inertia(1)*kb*temp/haba^2).*sqrt(2*Principal_Inertia(2)*kb*temp/haba^2).*sqrt(2*Principal_Inertia(3)*kb*temp/haba^2));
end