%function quadreach()
clear
elltool.setconf('accurate')

use_old=true;
if use_old
    load result
else

% g = 9.81;       % gravity (m/s^2)
% m = 1;          % mass (kg)
% L = 0.25;       % length of rotor arm (m)
% J = 20*m*L^2;   % moment of inertia (kg m^2) (est)

v0 = 1;     % speed m/s
Nd = 7;    % discretization for angle;
T_end = 4;

th0_A=linspace(-pi/2,pi/2,Nd+1)+pi/2/Nd;
thg_A=linspace(-pi/2,pi/2,Nd+1);

th0_B=linspace(pi/2,3*pi/2,Nd+1)+pi/2/Nd;
thg_B=linspace(pi/2,3*pi/2,Nd+1);

Ac_A = cell(Nd,1);
for i=1:Nd
    Ac_A{i} = [0, 0, -v0*sin(th0_A(i));
        0, 0,  v0*cos(th0_A(i));
        0, 0, 0];
end
Ac_B = cell(Nd,1);
for i=1:Nd
    Ac_B{i} = [0, 0, -v0*sin(th0_B(i));
        0, 0,  v0*cos(th0_B(i));
        0, 0, 0];
end
Bc = [0;0;1];

% control bounds
Ku = 0.1*[1; 1];
centVec = -diff(Ku)/2;
shMat = diag((Ku(1,:)-centVec).^2,0);
uBoundsEllObj = ellipsoid(centVec', shMat);

% initial directions (some random vectors in R^4):
%dirsMat = [diag([1 1 pi/6]),rand(3,6)];
%dirsMat=[[eye(3);zeros(7,3)],rand(10,7)];
load('result.mat','dirsMat')

% uBoundsEllObj=uBoundsEllObj.getShape(1.5);
% plot(uBoundsEllObj)

% linear system for system A
for i=1:Nd
    lsys_A{i} = elltool.linsys.LinSysContinuous(Ac_A{i}, Bc, uBoundsEllObj);  
end
% time interval
timeVec = [0 T_end];

% initial conditions:
th0_A=0;
th0_B=-pi;
x0EllObj_A = [zeros(2,1); th0_A] + ellipsoid(diag([0.01,0.01,0.0001]));
x0EllObj_B= [20*ones(2,1); th0_B] + ellipsoid(diag([0.01,0.01,0.0001]));

% reach set
startI_=find(thg_A<=th0_A);
startI=startI_(end);
rsObj_A = elltool.reach.ReachContinuous(lsys_A{startI}, x0EllObj_A, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);  
grdHypObj_1 = hyperplane([0; 0; 1], thg_A(startI+1));
grdHypObj_2 = hyperplane([0; 0; -1], -thg_A(startI));
% %rsObj_B = elltool.reach.ReachContinuous(lsys, x0EllObj_B, dirsMat, timeVec,...
% %    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7); 

plotA=rsObj_A.cut(T_end);
plObj=plotA.plotByEa('g'); % to have the use of plObj isn't necessary 
hold on
plotA.plotByIa('r',plObj);

%project to see the whoe reachable set for x y dimension
basisMat = [1 zeros(1,2); 0 1 zeros(1,1)]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
%psObj_B = rsObj_B.projection(basisMat);  % reach set projection

% to have the use of plObj isn't necessary 
plObj=psObj_A.plotByEa('g');  % external apprx. of reach set 1 (red)
hold on
psObj_A.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)
%psObj_B.plotByEa('y',plObj);
%psObj_B.plotByIa('b',plObj);
save result;
end
%%
[exEllMat_A, t_axis]= rsObj_A.get_ea();

exEllMat_A_tmp = exEllMat_A(1,:);
for i=2:size(exEllMat_A,1);
    exEllMat_A_tmp = exEllMat_A_tmp.intersection_ea(exEllMat_A(i,:));
end

intersectEllVec = exEllMat_A_tmp.hpintersection(grdHypObj_1);
indNonEmptyVec = ~isEmpty(intersectEllVec);

indNonEmptyVec = find(indNonEmptyVec); 
min(indNonEmptyVec)
max(indNonEmptyVec)
%%
crsObjVec = [];
for iInd = 1:size(indNonEmptyVec, 2)
    curTimeLimVec=[t_axis(indNonEmptyVec(iInd)-1) T_end];
     rsObj = elltool.reach.ReachContinuous(lsys_A{startI+1},...
         intersectEllVec(indNonEmptyVec(iInd)), ...
             dirsMat, curTimeLimVec,'isRegEnabled',true, 'isJustCheck', false, 'regTol', 1e-7);
     crsObjVec = [crsObjVec rsObj];
end

%%
% save result

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