% test derivatives


clear
i = 1;
j = 2;
bus_angle(i) = 0.0;
bus_angle(j) = 45*(pi/180);
bus_angle
v(i) = 1.0;
v(j) = 0.9;
%v(j) = 1.0;
v

ibranch = 1;
R(ibranch) = 0.1;
%R(ibranch) = 0.0;
X(ibranch) = 0.2;
Bcap(ibranch) = 0.04;
%Bcap(ibranch) = 0.0;
Gline =  R(ibranch)/(R(ibranch)^2 + X(ibranch)^2)
Bline = -X(ibranch)/(R(ibranch)^2 + X(ibranch)^2)
Bcapline = Bcap(ibranch)/2

P = Gline*v(i)^2 - Gline*v(i)*v(j)*cos(bus_angle(i)-bus_angle(j))-Bline*v(i)*v(j)*sin(bus_angle(i)-bus_angle(j))

Q = Bline*v(i)*v(j)*cos(bus_angle(i) - bus_angle(j)) - (Bcapline*v(i)^2) - Bline*v(i)^2 - Gline*v(i)*v(j)*sin(bus_angle(i) - bus_angle(j))


yline = 1/( R(ibranch) + sqrt(-1)*X(ibranch) )

ycap = sqrt(-1)*Bcap(ibranch)/2

vfrom = v(i)*( cos(bus_angle(i)) + sqrt(-1)*sin(bus_angle(i)) )
vto   = v(j)*( cos(bus_angle(j)) + sqrt(-1)*sin(bus_angle(j)) )

Iline = (vfrom - vto)*yline + vfrom*ycap

PjQ = vfrom*conj(Iline)

% >> syms Gline Bline v(i) v(j) angle(i) angle(j) Bcapline
% 
% >> Q = -1*v(i)*v(j)*sin(angle(i)-angle(j))*Gline - v(i)^2*Bline + v(i)*v(j)*cos(angle(i)-angle(j))*Bline - (v(i)^2)*Bcapline
%  
%  
% >> diff(Q,angle(i))
%  
% - Gline*v(i)*v(j)*cos(angle(i) - angle(j)) - Bline*v(i)*v(j)*sin(angle(i) - angle(j))
%  
% >> diff(Q,v(i))
%  
% Bline*v(j)*cos(angle(i) - angle(j)) - 2*Bcapline*v(i) - 2*Bline*v(i) - Gline*v(j)*sin(angle(i) - angle(j))
%  
% >> diff(Q,angle(j))
%  
% Gline*v(i)*v(j)*cos(angle(i) - angle(j)) + Bline*v(i)*v(j)*sin(angle(i) - angle(j))
%  
% >> diff(Q,v(j))
%  
% Bline*v(i)*cos(angle(i) - angle(j)) - Gline*v(i)*sin(angle(i) - angle(j))
% 
% 
