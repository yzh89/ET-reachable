%function quadreach()
clear
elltool.setconf('accurate')

use_old=false;
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
    
    v0 = P.Va0;     % speed m/s
    Nd = 7;    % discretization for angle;
    T_end = 4;
    
    P_lat = zeros(8,12);
    P_lat(1,1)=1; % x
    P_lat(2,2)=1; % y
    P_lat(3,4)=1; % u 
    P_lat(4,5)=1; % v
    P_lat(5,10)=1; %p
    P_lat(6,12)=1; %r
    P_lat(7,7)=1; % phi
    P_lat(8,9)=1; % psi
    
    
    x_0_lat = P_lat*x_trim;
    %x_0_lat = [0 0 16 0 -0.0298 0.1073 0.1957 0]';
    u_lat_cent = [0 1 0 0; 0 0 1 0]*u_trim;
    %u_lat_cent = [0.0316, -0.0425]';
    Bc_ = B_lat*[0 1 0 0; 0 0 1 0]';
    Bc = [zeros(3, size(Bc_,2)); Bc_];
    
    
    phi0=x_0_lat(7);
    psi0_=x_0_lat(8);
    psi0=linspace(-pi/2+psi0_,pi/2+psi0_,Nd+1)+pi/2/Nd;
    psig=linspace(-pi/2+psi0_,pi/2+psi0_,Nd+1);

    Ac = cell(Nd,1);
    for i=1:Nd
        Ac_ = [ 0 0 cos(psi0(i)) -cos(phi0)*sin(psi0(i)); ...
                0 0 sin(psi0(i)) cos(phi0)*cos(psi0(i)); ...
                0 0 0 0];
        Ac{i} = [Ac_ zeros(3, size(A_lat,1)-1);zeros(size(A_lat,1),3),A_lat];
    end
    
    % control bounds
    Ku = [pi/4+u_lat_cent(1),pi/4+u_lat_cent(2);...
          pi/4-u_lat_cent(1), pi/4-u_lat_cent(2)];
    centVec = -diff(Ku)/2;
    shMat = diag((Ku(1,:)-centVec).^2,0);
    uBoundsEllObj = ellipsoid(centVec', shMat);
    
    % initial directions (some random vectors in R^4):
    %dirsMat = [diag([1 1 pi/6]),rand(3,6)];
    dirsMat=[eye(8,8),rand(8,8)];
    %load('result_fixedwing_0.mat','dirsMat')
    
    % uBoundsEllObj=uBoundsEllObj.getShape(1.5);
    % plot(uBoundsEllObj)
    
    % linear system for system A
    Amat = [0 0 1 0 0 0 0 0 ;
0 0 0 1 0 0 0 0;
0 0 0 0 0 0 0 0;
0 0 0 -0.43 4.55 -16.37 9.27 0;
0 0 0 -2.1 -7.8 3.51 0 0;
0 0 0 2.29 -0.23 -4.7 0 0;
0 0 0 0 1 0.26 0 0;
0 0 0 0 0 1.01 0 0;
];

Bmat = [0 0;0 0; 0 0; 0 -1.26; 30 36.7; 12 -2.8; 0 0; 0 0];
    for i=1:Nd
        lsys_A{i} = elltool.linsys.LinSysContinuous(Ac{i}, Bc, uBoundsEllObj);
    end
    % time interval
    timeVec = [0 T_end];
    
    lsys = elltool.linsys.LinSysContinuous(Amat, Bmat, uBoundsEllObj);
    
    % initial conditions:
    x0EllObj = [x_0_lat] + ellipsoid(diag([0.01,0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]));
    
    % reach set
    startI_=find(psig<=psi0_);
    startI=startI_(end);
    rsObj_A = elltool.reach.ReachContinuous(lsys, x0EllObj, dirsMat, timeVec,...
        'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);
    grdHypObj_1 = hyperplane([0; 0; 1], psig(startI+1));
    grdHypObj_2 = hyperplane([0; 0; -1], -psig(startI));
    % %rsObj_B = elltool.reach.ReachContinuous(lsys, x0EllObj_B, dirsMat, timeVec,...
    % %    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);
    save result_fixedwing_0;
    
end

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

basisMat = [1 zeros(1,2); 0 0 1]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
%psObj_B = rsObj_B.projection(basisMat);  % reach set projection

% to have the use of plObj isn't necessary
plObj=psObj_A.plotByEa('g');  % external apprx. of reach set 1 (red)
hold on
psObj_A.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)

%%
if use_old_inter
    load result_fixedwing
else
    [exEllMat_A, t_axis]= rsObj_A.get_ea();
    
    intersectEllVec = exEllMat_A.hpintersection(grdHypObj_1);
    indNonEmptyVec = all(~isEmpty(intersectEllVec));
    
    indNonEmptyVec = find(indNonEmptyVec);
    min(indNonEmptyVec)
    max(indNonEmptyVec)
    
%     crsObjVec = [];
%     for iInd = 1:size(indNonEmptyVec, 2)
%         intersectEllVec_tmp = intersectEllVec(1,indNonEmptyVec(iInd));
%         for i=2:size(intersectEllVec,1);
%             intersectEllVec_tmp = intersectEllVec_tmp.intersection_ea(intersectEllVec(i,indNonEmptyVec(iInd)));
%         end
%         curTimeLimVec=[t_axis(indNonEmptyVec(iInd)-1) T_end];
%         rsObj = elltool.reach.ReachContinuous(lsys_A{startI+1},...
%             intersectEllVec_tmp, ...
%             dirsMat, curTimeLimVec,'isRegEnabled',true, 'isJustCheck', false, 'regTol', 1e-7);
%         crsObjVec = [crsObjVec rsObj];
%     end
%     save('result_fixedwing','crsObjVec');
end
%%
if plotting
basisMat = [1 0 0; 0 1 0]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
%psObj_B = rsObj_B.projection(basisMat);  % reach set projection

% to have the use of plObj isn't necessary
plObj=psObj_A.plotByEa('g'); 
hold on
% psObj_A.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)

for i=1:5:size(crsObjVec,2)
    psObj_A = crsObjVec(i).projection(basisMat);  % reach set projection
    %psObj_B = rsObj_B.projection(basisMat);  % reach set projection
    
    % to have the use of plObj isn't necessary
    psObj_A.plotByEa('y',plObj);  % external apprx. of reach set 1 (red)
%     psObj_A.plotByIa('b',plObj);  % internal apprx. of reach set 1 (green)
end

end

%%
[x0,x0shMat]=x0EllObj.double();
% [ctrMat, ttVec] = rsObj_B.get_center();
xB=[4,0.5,-pi]';
[qc,Qc]=findControlSet_nonlinear(x0,x0shMat,xB,Ac_A{startI},Bc,centVec',shMat);

plot(ellipsoid(qc{4},Qc{4}'*Qc{4}))

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