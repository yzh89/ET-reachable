function [sdot] = sysQuad(t,s,u)

sdot = zeros(13,1);
x = s(1);
y = s(2);
z = s(3);
xdot = s(4);
ydot = s(5);
zdot = s(6);
qW = s(7);
qX = s(8);
qY = s(9);
qZ = s(10);
p = s(11);
q = s(12);
r = s(13);

m = 1;
g = 9.81;
L = 0.25;
I = [20*m*L^2,0,0;0,20*m*L^2,0;0,0,20*m*L^2];
params.mass = m;
params.I = I;
params.invI = inv(I);
params.grav = g;

Rot = quattoRot([qW,qX,qY,qZ]');
[phi,theta,yawangle] = RotToRPY_ZXY(Rot);

BRW = [ cos(yawangle)*cos(theta) - sin(phi)*sin(yawangle)*sin(theta), cos(theta)*sin(yawangle) + cos(yawangle)*sin(phi)*sin(theta), -cos(phi)*sin(theta);...
    -cos(phi)*sin(yawangle),  cos(phi)*cos(yawangle),  sin(phi);...
    cos(yawangle)*sin(theta) + cos(theta)*sin(phi)*sin(yawangle), sin(yawangle)*sin(theta) - cos(yawangle)*cos(theta)*sin(phi),  cos(phi)*cos(theta)];

WRB = BRW';
sdot(1) = xdot;
sdot(2) = ydot;
sdot(3) = zdot;

accel = 1/params.mass *(WRB * [0;0;u(1)+params.mass*g] - [0;0;params.mass*params.grav]);

sdot(4) = accel(1);
sdot(5) = accel(2);
sdot(6) = accel(3);

omega = [p;q;r];

%this enforces the magnitude 1 constraint for the quaternion
K_quat = 2;
quaterror = 1 - (qW^2 + qX^2 + qY^2 + qZ^2); 

qdot = -1/2*[0,-p,-q,-r;...
    p,0,-r,q;...
    q,r,0,-p;...
    r,-q,p,0]*[qW,qX,qY,qZ]' + K_quat*quaterror*[qW,qX,qY,qZ]';

sdot(7) = qdot(1);
sdot(8) = qdot(2);
sdot(9) = qdot(3);
sdot(10) = qdot(4);

pqrdot = params.invI * ([u(2);u(3);u(4)]*L - cross(omega,params.I*omega));

sdot(11) = pqrdot(1);
sdot(12) = pqrdot(2);
sdot(13) = pqrdot(3);