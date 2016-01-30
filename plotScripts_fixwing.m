grid off
title('Reachable Tube in x z space for time interval [0,10]s','FontSize',16)
hXLabel = xlabel('time');
hYLabel = ylabel('x (m)');
hZLabel = zlabel('z (m)');
set([hXLabel, hYLabel, hZLabel]  , ...
    'FontSize'   , 16          );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 12);

grid on

%%
% reachable set
grid off
title('Reachable Set at t=4s','FontSize',16)
hXLabel = xlabel('x (m)');
hYLabel = ylabel('y (m)');
hZLabel = zlabel('z (m)');
set([hXLabel, hYLabel, hZLabel]  , ...
    'FontSize'   , 16          );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 16);

grid on
%%
%control set
grid off
title('Constraint Control Sets for aircraft A and B','FontSize',16)
hXLabel = xlabel('\delta_e (N)');
hYLabel = ylabel('\delta_f (rad)');
set([hXLabel, hYLabel, hZLabel]  , ...
    'FontSize'   , 16          );
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'LineWidth'   , 1         , ...
  'FontSize'    , 16);

grid on

