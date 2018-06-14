
function [] = test 
clc 
clear all 
close all 
%%%%
%NARRATIVE
input = 2 ;
%%%%
%INITIAL CONSTANTS 
Re = 6378e3;         %Earth Radius 
B0 = 3.07e-5;        % Magnetic equator field strength in T
Md = [0,0,-B0*Re^3]; % Magnetic moment
Rm = [0,0,0]; 
Rtesting = [ 4*Re, 0, 0];     %Original Position
mp = 1.6726231e-27;  % proton mass 
qp = 1.60217733e-19; % proton charge
eV = 1.60217733e-19; % conversion factor eV to J
c  = 299792458.0;    % speed of light
QpoverMp = qp/mp; 
Ep = 10e6*eV;       % 10 MeV proton
vp = c*sqrt(Ep*(Ep+2*mp*c^2))/(Ep+mp*c^2); % relativistic velocity
V2 = vp^2;
gamma = 1/sqrt(1-V2/c^2); % relativistic factor
alpha = 15; %Pitch angle defined as arctan(vpara/vperpen)
%%---%%

%%%%
%DEDUCTIONS FROM THE INITIAL PARAMETERS 
%Calling auxiliary functions 
[B, B2, b, dB] = GiveMagDip(Md,Rm, Rtesting);
[Tg, Tb, Td] = GivePeriods(Rtesting);
lm = GiveMirrorLat(alpha);
t =[0, input*Tb] ;

v = [0 vp*sind(alpha) vp*cosd(alpha)]; %velocity vector assuming x component to be 0 
vpara = v*b'; 
vperpen = sqrt(V2-vpara^2); 
input = [Rtesting, v];
mu = gamma^2*mp*vperpen^2/(2*sqrt(B2)); %First invariant 

rho = abs(gamma/QpoverMp*vperpen/(sqrt(B2))); %Gyro radius 
Rgc = [4*Re- sign(QpoverMp*Md(3))*rho; 0; 0]  ;  %Guiding centre position and velocity
Vgc = vpara; 
inputgc = [Rgc' , Vgc];

L = norm(Rtesting)/Re ; %L shell
%%%%
%Loss Cone alert
alphaloss = asind(sqrt(1/sqrt(4*L^6-3*L^5))); %angle (Latitude) at which the mirror point falls on the Earth. 
fprintf(1,'|alpha| %.2f |pi-alpha| %.2f loss cone %.2f deg\n', ...
        abs(alpha),abs(180-alpha),alphaloss);
if abs(alpha)<alphaloss || abs(180-alpha)<alphaloss,
  fprintf(1,'warning: particle in loss cone!\n');
end



%%%%
%PLOTING 
trace(input,inputgc)

%%%
%MOTION TRACING
options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','off','vectorized','on');


function trace(input,inputgc)
%%Define planet surface 3d sphere 
%%I have had to add the options in a differnet way 

[t,output] = ode15s(@motion, t, input,odeset('AbsTol',1e-12,'RelTol',1e-12,'Stats','off','vectorized','on'));
Z = output(3);
Rg = sqrt(output(1)^2+output(2)^2);
Rtot = sqrt(output(1)^2+output(2)^2+output(3)^2);

[tgc, outputgc] = ode15s(@motiongc, t, inputgc);
Zgc = outputgc(3)
Rggc = sqrt(outputgc(1)^2+outputgc(2)^2)
Rtotgc = sqrt(outputgc(1)^2+outputgc(2)^2+outputgc(3)^2)


end         

%%%%
%AUXILIARY FUNCTIONS 

 function [B, B2, b, dB] = GiveMagDip(m, rm, r )
%Preliminaires 
R = r-rm ;
R2= R*R' ;
R3 = (sqrt(R2))^1.5;
iR3 = 1/R3;
r = R/sqrt(R2);
mr  = m*r' ;

%%%
%B 
B = iR3.*(3*mr.*r - m);
B2 = B*B' ;
b = B./sqrt(B2);

%%%%
%dB

f1 = 3./R3.^2;
u = sqrt(3*mr.^2 + m*m');
f2 = R2.*mr./u ;
f3 = R.*u ;
db1 = f1.*(f2.*m(1).*(1-r(1).^2) - R(1).*f3(1));
db2 = f1.*(f2.*m(2).*(1-r(2).^2) - R(2).*f3(2));
db3 = f1.*(f2.*m(3).*(1-r(3).^2) - R(3).*f3(3));
dB = [db1 db2 db3];

end

function [Tg, Tb, Td] = GivePeriods(r)
%https://farside.ph.utexas.edu/teaching/plasma/lectures1/node22.html
%https://farside.ph.utexas.edu/teaching/plasma/lectures1/node23.html
%Preliminaries 
[B, B2]= GiveMagDip(Md, Rm, r);
[B, B2e]= GiveMagDip(Md, Rm, [Re 0 0] );
B = sqrt(B2);
Be = sqrt(B2e);
        
%%%%
%Gyro calculations 
Wg = QpoverMp/gamma*B ;
Tg = 2*pi./Wg; 

%%%%
%Bounce calculations 
Tb = 0.117*(norm(r)/Re)*(c/vp)*(1-0.4635*sind(alpha)^0.75);

%%%%
%Drift calculations 
Td = 2*pi*QpoverMp*Be*Re^3/vp^2/norm(r)*(1-1/3*(sind(alpha))^0.62);               
end 

function lm = GiveMirrorLat(alpha)
%The latitude at which Bm is reached according to the first invariant mu
eqn = @(x,a) cosd(x).^6-sind(a)^2.*sqrt(1+3*sind(x).^2);
lm = fzero(@(x) eqn(x,alpha),45);    
end 


%%%%
%ODE'S 
function out = motion(t,input)
r = [input(1),input(2),input(3)];
v = [input(4),input(5),input(6)];
B = GiveMagDip(Md, Rm, r);
fac = gamma./ QpoverMp/gamma; 
dv = fac*cross(v, B);
out = [v'; dv'];
end 

function out = motiongc(t,input)
r = [input(1),input(2),input(3)];
v = input(4); 
[B, B2, b, dB] = GiveMagDip(Md, Rm, r);
facdv = mu/gamma^2/mp; %eqn 23) Ozturk
fac = gamma./(2*QpoverMp*B2).*(V2 + v.^2);
dr = fac.*cross(b, dB) + v.*b;
dv = -facdv*(b*dB');  %%Yields 0 might be wrong
out = [dr dv]';


end 
            
end 
