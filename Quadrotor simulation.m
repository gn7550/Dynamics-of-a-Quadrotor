% prob 5.a
%     %(*single rotor will cause only a single direction of torque and 
%   counter-torque by the 3rd property of Newton's law(action and reaction), so the 
% single rotor will rotate itself, which means impossible to control it's orientation. *)

% Prob 5.b
% f = f1 + f2 + f3 + f4; 
% torque = [d*(f2 - f4);
%           d*(f3 - f1);
%           c*(-f1 + f2 - f3 + f4)];
% 		(* if a pair of tow1 and tow3 and a pair of tow2 and tow4 have the 
% same direction of torque,the quadrotor will keep rotating about z axis no matter what 
% thrust force rotors have, for the same reason in prob 5.a. So, each pair has to 
% have opposite direction of torque so that you can control its orientation 
% by adjusting thrust forcec of the pairs. *)

% prob 5.d
clear all
clc
close all
I=[0.082,0,0; 0,0.0845,0 ; 0,0,0.1377];
m=4.34; d=0.3; c=8*10^-4; g=9.81; p=[0;0;0];
x=[0;0;0]; xdot=[0;0;0]; w=[0;0;0]; theta=[0;0;0]; thetadot=[0;0;0];
t=0:0.01:3 ; dt=0.01;
R(:,:,1)=eye(3);
for i=1:length(t)
t=i*0.01-0.01;
f1= abs(sin(2*pi*t))+1/4*m*g; f2= abs(cos(pi*t))+1/4*m*g;
f3= abs(cos(2*pi*t))+1/4*m*g; f4= abs(sin(pi*t))+1/4*m*g;
f=[0;0;f1+f2+f3+f4];  %thrust vector

taw= [d*(f2-f4);d*(f3-f1);c*(-f1+f2-f3+f4)];


  
 wdot= inv(I)*(taw-cross(w(:,i),I*w(:,i)));
 w(:,i+1)= w(:,i) + dt*wdot;
 omega(:,:,i+1)= skew( w(:,i+1) );
 Rdot(:,:,i)= R(:,:,i)*omega(:,:,i+1);
  
 R(:,:,i+1)= R(:,:,i)*expm(dt*omega(:,:,i+1));
 M=[m,0,0;0,m,0;0,0,m];
 a= inv(M)*( f-cross(w(:,i+1),m*xdot) );
 
 xdot = xdot + dt*a;
 x(:,i+1) = x(:,i) + dt * xdot; 
 p(:,i+1)= R(:,:,i)* x(:,i+1); %center of mass position
end
t=0:0.01:3.01;
X=p(1,:); Y=p(2,:); Z=p(3,:);
 
figure
plot3(X,Y,Z,'r')
xlim([0 3])
ylim([-7 0])
zlim([0 50])
grid on
hold on
for i=1:11
X11(i)=X(30*i-29);
Y11(i)=Y(30*i-29);
Z11(i)=Z(30*i-29);
q11(:,i)=[X11(i),Y11(i),Z11(i)]';
R11(:,:,i)=R(:,:,30*i-29);
    xi=q11(:,i)+R11(:,1,i)/2;
    xg=q11(:,i)-R11(:,1,i)/2;
    yi=q11(:,i)+R11(:,2,i)/2
    yg=q11(:,i)-R11(:,2,i)/2
    zi=q11(:,i)+R11(:,3,i)/2;
    zg=q11(:,i)-R11(:,3,i)/2;
plot3([xi(1),xg(1)],[xi(2),xg(2)],[xi(3),xg(3)],[yi(1),yg(1)],[yi(2),yg(2)],[yi(3),yg(3)])
end
t=0:0.1:3.01;
scatter3(X11,Y11,Z11,'r')
xlim([0 3])
ylim([-7 0])
zlim([0 50])
grid on
% its orientaion in R11 matrix


function N=skew(n)
N=[0 -n(3) n(2) ; n(3) 0 -n(1); -n(2) n(1) 0 ] ;
end
