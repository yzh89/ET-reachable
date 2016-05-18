
clear

figure (1)

phi0 = 0.0; theta0 = 0.0; psi0 = 0.0;
Rot0 = RPYtoRot_ZXY(phi0,theta0,psi0);
Quat0 = RotToQuat(Rot0);

u = [0.1898, 0.09979, -1.9e-15,0;
    0.0258, 0.09979, -1.59e-15,0;
    0.1078, 0.054, -1.84e-15, 0;
    0.1137, 0.1443, -0.009618, 0;
    0.1078, 0.09979, 0.0819, 0;
    %0.1078,0.0933907,-1.35284e-15,0;
    0.0986, 0.1034, -0.0812,0]';
[r,c]=size(u);

for i=1:c
    [tout,yout] = ode45(@(s,t) sysQuad(s,t,u(:,i)),0:0.1:4,[0;0;0 ; 0.2 ; 0; 0; Quat0; zeros(3,1)]);
    plot3(tout,yout(:,1),yout(:,2),'k--');
end

%% aircraft b

phi0 = 0.0; theta0 = 0.0; psi0 = 0;
Rot0 = RPYtoRot_ZXY(phi0,theta0,psi0);
Quat0 = RotToQuat(Rot0);

u2 = [0.1078 -0.0626, 0,0;
    0.1078 -0.1448,0, 0;
    0.1857, -0.1037, 0, 0;
    0.105, -0.1037, 0,0;
    0.0299, -0.1037, 0, 0;
    0.1078, -0.1037, 0.078, 0;
    0.1078, -0.1037, -0.0778,0]';
[r,c]=size(u2);

for i=1:c
    [tout2,yout2] = ode45(@(s,t) sysQuad(s,t,u2(:,i)),0:0.1:4,[1.6;0.5;0 ; -0.2 ; 0; 0; Quat0; zeros(3,1)]);
    plot3(tout2,yout2(:,1),yout2(:,2),'r--');
end
