%function quadreach()
clear
elltool.setconf('accurate')

use_old=true;
use_old_inter = false;
plotting = false;
if use_old
    load result_fixedwing_0
else
    
    % g = 9.81;       % gravity (m/s^2)
    % m = 1;          % mass (kg)
    % L = 0.25;       % length of rotor arm (m)
    % J = 20*m*L^2;   % moment of inertia (kg m^2) (est)
    param_chap5 % this script generate a trim linearization of fixed wing dynamics
    
    % Nd = 7;    % discretization for angle;
    T_end = 10;
    
    
    P_lon = zeros(6,12);
    P_lon(1,1)=1; % x
    P_lon(2,4)=1; % u 
    P_lon(3,6)=1; % w
    P_lon(4,11)=1;% q
    P_lon(5,8)=1; % theta
    P_lon(6,3)=1; % z
    
    %x_0_lat = [0 0 16 0 -0.0298 0.1073 0.1957 0]';
    %u_lon_cent = [1 0 0 0; 0 0 0 1]*u_trim;
    %u_lat_cent = [0.0316, -0.0425]';
    Bc_ = B_lon*[1 0 0 0; 0 0 0 1]';
    Bc = [zeros(1, size(Bc_,2)); Bc_];
    
    Ac = [0 1 0 0 0 0; zeros(5,1), A_lon];

    % control bounds
    %Ku = [pi/10+u_lon_cent(1),0.01+u_lon_cent(2);...
    %      pi/10-u_lon_cent(1),0.01-u_lon_cent(2)];
    Ku = [0.01,0.01;...
          0.01,0.01];
    centVec = -diff(Ku)/2;
    shMat = diag((Ku(1,:)-centVec).^2,0);
    uBoundsEllObj = ellipsoid(centVec', shMat);
    
    % initial directions (some random vectors in R^4):
    %dirsMat = [diag([1 1 pi/6]),rand(3,6)];
    dirsMat=[eye(6,6),rand(6,3)];
    %load('result_fixedwing_0.mat','dirsMat')
    
    % uBoundsEllObj=uBoundsEllObj.getShape(1.5);
    % plot(uBoundsEllObj)
    
    % linear system for system A

    lsys_A = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj);
    % time interval
    timeVec = [0 T_end];
    
    % initial conditions:
    %x0EllObj = [x_0_lon] + ellipsoid(diag([0.1,0.01, 0.01, 0.001, 0.001, 0.1]));
    x0EllObj_A =  ellipsoid(diag([0.1,0.01, 0.01, 0.001, 0.001, 0.1]));
        
    % reach set
    rsObj_A = elltool.reach.ReachContinuous(lsys_A, x0EllObj_A, dirsMat, timeVec,...
        'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-5);
    rsObj_B = rsObj_A.getCopy();
    
    % %rsObj_B = elltool.reach.ReachContinuous(lsys, x0EllObj_B, dirsMat, timeVec,...
    % %    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);
    save result_fixedwing_0;
    
end
%%
simOut = sim('mavsim_chap5_no_plot','SaveState','on','StateSaveName','xout');
x_star = get(simOut,'xout');

simOut_B = sim('mavsim_chap5_no_plot', 'InitFcn', 'param_chap5B', 'SaveState','on','StateSaveName','xout_B');
x_B_star = get(simOut_B, 'xout_B');

basisMat = [1 zeros(1,5); zeros(1,5) 1]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
psObj_B = rsObj_B.projection(basisMat);  % reach set projection

psObj_copy = psObj_A.getCopyWithCenterModified(x_star(1:end-1,[1,3]));
psObj_copy_B = psObj_B.getCopyWithCenterModified(x_B_star(1:end-1,[1,3]));

plotA=psObj_copy.cut(T_end);
plotB=psObj_copy_B.cut(T_end);

plObj=plotA.plotByEa('g'); % to have the use of plObj isn't necessary
hold on
plotA.plotByIa('r',plObj);
plotB.plotByEa('y',plObj);
plotB.plotByIa('b',plObj);

% to have the use of plObj isn't necessary
plObj=psObj_copy.plotByEa('g');  % external apprx. of reach set 1 (red)
hold on
psObj_copy.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)
psObj_copy_B.plotByEa('y',plObj);
psObj_copy_B.plotByIa('b',plObj);

plot3(0:0.1:10,x_star(:,1)',x_star(:,3)','k');
plot3(0:0.1:10,x_B_star(:,1)',x_B_star(:,3)','k--');

%%
[x0,x0shMat]=x0EllObj_A.double();
[ctrMat, ttVec] = rsObj_B.get_center();
xB = [0, zeros(1,4), 2]';
[qc,Qc]=findControlSet_nonlinear(x0,x0shMat,xB,Ac,Bc,centVec',shMat,100);

%%
iter=1;
q=qc{iter};
Q=Qc{iter};
uBoundsEllObj_cu = ellipsoid(q, Q);

lsys_cu = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj_cu);
rsObj_A = elltool.reach.ReachContinuous(lsys_cu, x0EllObj_A, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-5);
basisMat = [1 zeros(1,5); zeros(1,5) 1]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection

psObj_copy = psObj_A.getCopyWithCenterModified(x_star(1:end-1,[1,3]));

plotA=psObj_copy.cut(T_end);
plObj=plotA.plotByEa('g'); % to have the use of plObj isn't necessary
hold on
plotA.plotByIa('r',plObj);
plotB.plotByEa('y',plObj);
plotB.plotByIa('b',plObj);

% to have the use of plObj isn't necessary
plObj=psObj_copy.plotByEa('g');  % external apprx. of reach set 1 (red)
hold on
psObj_copy.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)
psObj_copy_B.plotByEa('y',plObj);
psObj_copy_B.plotByIa('b',plObj);
plot3(0:0.1:10,x_star(:,1)',x_star(:,3)');
plot3(0:0.1:10,x_B_star(:,1)',x_B_star(:,3)','k--');

save result_fixwing_bkup


% %%
%
% [x0,x0shMat]=x0EllObj_A.double();
% [ctrMat, ttVec] = rsObj_B.get_center();
% xB=ctrMat(:,end);
% [qc,Qc]=findControlSet_v2(x0,x0shMat,xB,Ac,Bc,centVec',shMat);
%
%
% %%
% %for iter=1:length(qc)
% iter=6;
%     q=qc{iter};
%     Q=Qc{iter};
%  uBoundsEllObj_cu = ellipsoid(q, Q);
%  %uBoundsEllObj=uBoundsEllObj.getShape(1.5);
% % plot(uBoundsEllObj)
% % linear system
% lsys_cu = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj_cu);
% %dirsMat = [[eye(3);zeros(7,3)],rand(10,7)];
% rsObj_A = elltool.reach.ReachContinuous(lsys_cu, x0EllObj_A, dirsMat, timeVec,...
%     'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-8);
%
% basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
% psObj_A = rsObj_A.projection(basisMat);  % reach set projection
% psObj_B = rsObj_B.projection(basisMat);  % reach set projection
% % plot projection of reach set external approximation:
%
% % to have the use of plObj isn't necessary
% plObj=psObj_A.plotByEa('r');  % external apprx. of reach set 1 (red)
% %firstPsObj.plotByEa('r');
% hold on
% psObj_A.plotByIa('g',plObj);  % internal apprx. of reach set 1 (green)
% psObj_B.plotByEa('y',plObj);
% psObj_B.plotByIa('b',plObj);
% %end
%
% save result_bkup