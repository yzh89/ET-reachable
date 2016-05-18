
clear

figure (1)

phi0 = 0.0; theta0 = 0.0; psi0 = 0.0;
Rot0 = RPYtoRot_ZXY(phi0,theta0,psi0);
Quat0 = RotToQuat(Rot0);

u = [-0.1818, 0, 0,0;
    0.3874, 0, 0,0;
    0.1078, 0.145, 0, 0;
    0.1078, -0.145, 0, 0;
    0.1078, 0, 0.1448, 0;
    0.1078, 0, -0.1448,0]';
[r,c]=size(u);

for i=1:c
    [tout,yout] = ode45(@(s,t) sysQuad(s,t,u(:,i)),0:0.1:4,[0;0;0 ; 0.2 ; 0; 0; Quat0; zeros(3,1)]);
    plot3(tout,yout(:,1),yout(:,2),'k--');
end

%% aircraft b

phi0 = 0.0; theta0 = 0.0; psi0 = 0;
Rot0 = RPYtoRot_ZXY(phi0,theta0,psi0);
Quat0 = RotToQuat(Rot0);

for i=1:c
    [tout2,yout2] = ode45(@(s,t) sysQuad(s,t,u(:,i)),0:0.1:4,[1.6;0.5;0 ; -0.2 ; 0; 0; Quat0; zeros(3,1)]);
    plot3(tout2,yout2(:,1),yout2(:,2),'r--');
end
