clear
elltool.setconf('accurate')
use_old=true;
if use_old
    load time_varying_part_0
else
    g = 9.81;       % gravity (m/s^2)
    M = 1;          % mass (kg)
    L = 0.25;       % length of rotor arm (m)
    J = 20*m*L^2;   % moment of inertia (kg m^2) (est)
    
    T = 6;
    
    load('init_param_time_varying.mat'); % this gives initialization parameter and dirc vector
    
    % Define System Matrix
    tmp = [0 g; -g 0; 0 0];
    Ac = [zeros(3),   eye(3),     zeros(3,2), zeros(3,2);
        zeros(3),   zeros(3),   tmp,        zeros(3,2);
        zeros(2,3), zeros(2,3), zeros(2),   eye(2);
        zeros(2,3), zeros(2,3), zeros(2),   zeros(2)];
    
    tmp = [0; 0; 1/M];
    Bc = [zeros(3,1) zeros(3,2);
        tmp        zeros(3,2);
        zeros(2,1) zeros(2,2);
        zeros(2,1) L/J*eye(2)];
  
    % Define Control Set for ego vehicle
    Ku = 0.04*[9.935 3.62 3.62;   %+ boundary
        4.545 3.62 3.62];         %- boundary

    c_U = -diff(Ku)/2;
    M_U = diag((Ku(1,:)-c_U).^2,0);
    U_bar = ellipsoid(c_U', M_U);

    % plot(U_bar)
    
    % Define System
    U_i = U_i.getShape(0.6);
    [c_U_i, M_U_i] = U_i.double();
    U_i = U_i-c_U_i+c_U'; % remove offset in center from example U
    lsys_e_bar = elltool.linsys.LinSysContinuous(Ac, Bc, U_bar);
    lsys_i = elltool.linsys.LinSysContinuous(Ac, Bc, U_i);
    
    timeVec = [0 T];  % time interval% initial conditions:
    c_x0 = [-0.8; 1; 0; 0.16; -0.12; zeros(5,1)];
    x0_e= c_x0 + ellipsoid(diag([0.005,0.005,0.005,0.0001,0.0001,0.0001, zeros(1,4)]));
    c_x0_i = [1; 0; 0; -0.2; zeros(6,1)];
    x0_i= c_x0_i + ellipsoid(diag([0.01,0.01,0.01,0.001,0.001,0.001, zeros(1,4)]));
    
    % initial directions (some random vectors in R^4):
    % for now use simple thing to test.
    dirsMat=[[eye(3);zeros(7,3)],rand(10,7)];
    
    % reach set
    init_tube_e = elltool.reach.ReachContinuous(lsys_e_bar, x0_e, dirsMat, timeVec,...
        'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);
    tube_i = elltool.reach.ReachContinuous(lsys_i, x0_i, dirsMat, timeVec,...
        'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);
    
    %project to see reachable set at end time
    basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
    plotObj_e = init_tube_e.projection(basisMat);  % reach set projection
    plotObj_i = tube_i.projection(basisMat);
    
    for i = 1:T
        plotObj_e_{i} = plotObj_e.cut(i);
        plotObj_i_{i} = plotObj_i.cut(i);
        
        plObj=plotObj_e_{i}.plotByEa('g'); % to have the use of plObj isn't necessary
        hold on
        plotObj_e_{i}.plotByIa('r',plObj);
        plotObj_i_{i}.plotByEa('y',plObj);
        plotObj_i_{i}.plotByIa('b',plObj);
    end
    
    % plot projection of reach set external approximation:
    
    %project to see the whoe reachable set for x y dimension
    basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
    plotObj_e = init_tube_e.projection(basisMat);  % reach set projection
    plotObj_i = tube_i.projection(basisMat);  % reach set projection
    % plot projection of reach set external approximation:
    
    % to have the use of plObj isn't necessary
    plObj=plotObj_e.plotByEa('g');  % external apprx. of reach set 1 (red)
    hold on
    plotObj_e.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)
    plotObj_i.plotByEa('y',plObj);
    plotObj_i.plotByIa('b',plObj);
    
    save time_varying_part_0
end
%%
% dimension of the problem
m = 3; %input dimension
n = 10; %state dimension
% assign sampling and l(k) vector
N = 20;
dt = T/N;
l = zeros(n,N);

%initialize a sepBall
d=1;
sepBall=ellipsoid(zeros(10,1),diag([d^2,d^2,d^2 zeros(1,7)]));
dirsMatSimple=[[eye(3);zeros(7,3)],rand(10,7)];

for k = 1:N
    f = @(s) expm(Ac*(k*dt-s))*Bc*c_U_i;
    r_i = expm(Ac*k*dt)*c_x0_i + integral(f,0,k*dt,'ArrayValued',1);
    r_e = expm(Ac*k*dt)*c_x0 + integral(f,0,k*dt,'ArrayValued',1);
    l_tmp=[eye(3),zeros(3,n-3);zeros(n-3,n)]*(r_i-r_e);
    l(:,k)=l_tmp/norm(l_tmp);
    
    tube_i_ea_k=tube_i.cut(k*dt).get_ea;
    

    for iter=1:length(tube_i_ea_k)
        temp=[tube_i_ea_k(iter),sepBall];
        tube_i_ea_k_new(iter,:)=temp.minksum_ea(dirsMatSimple);
    end
    [a_,b_] = size(tube_i_ea_k_new);

    tube_i_ea_k_new=reshape(tube_i_ea_k_new,a_*b_,1);
    temp=tube_i_ea_k_new(1);
    for iter=1:length(tube_i_ea_k_new)-1
        temp=temp.intersection_ea(tube_i_ea_k_new(iter+1));
    end
    clear extEllipANew;
    tube_i_ea_k_new=temp;
    
    M_i{k} = tube_i_ea_k_new.getShapeMat;
    c_i{k} = tube_i_ea_k_new.getCenterVec;
    clear tube_i_ea_k_new;
    
end
M_x0 = x0_i.double();
save('time_varying_function_input','N','T','Ac','Bc','l','M_i','c_i','c_x0','M_x0','c_U','M_U');
