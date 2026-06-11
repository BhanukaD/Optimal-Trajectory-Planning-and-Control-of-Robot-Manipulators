clc; clear; close all;

%% UR5 CB3 Robot Dimensions

d1 = 0.089159;
d4 = 0;
d5 = 0;
d6 = 0.10915;
d7 = 0.09465;
d8 = 0.08230;

a2 = 0.42500;
a3 = 0.39225;

%% Modified Denavit-Hartenberg Parameterization

L(1) = RevoluteMDH('alpha', 0,       'a', 0,  'd', d1);
L(2) = RevoluteMDH('alpha', pi/2,    'a', 0,  'd', d4);
L(3) = RevoluteMDH('alpha', 0,       'a', a2, 'd', d5);
L(4) = RevoluteMDH('alpha', 0,       'a', a3, 'd', d6);
L(5) = RevoluteMDH('alpha', -pi/2,   'a', 0,  'd', d7);
L(6) = RevoluteMDH('alpha', pi/2,    'a', 0,  'd', d8);

%% Serial-Link Robot Construction

robot = SerialLink(L, 'name', 'UR5 Modified DH Robot');

%% Joint Limit Assignment

robot.qlim = deg2rad([
   -180  180
   -120  120
   -150  150
   -180  180
   -120  120
   -360  360
]);

q_home = deg2rad([0 0 0 0 0 0]);

disp(robot);
T_home = robot.fkine(q_home);
disp('End-effector pose at home:');
disp(T_home);

%% Workspace Sampling

N = 60000;
P = zeros(N,3);

qlim = robot.qlim;

for i = 1:N
    q = qlim(:,1)' + rand(1,6).*(qlim(:,2)' - qlim(:,1)');
    T = robot.fkine(q);

    if isa(T, 'SE3')
        P(i,:) = T.t';
    else
        P(i,:) = transl(T);
    end
end

%% Workspace Visualization

figure('Color','w','Position',[100 100 1000 850]);
hold on; grid on; axis equal;

scatter3(P(:,1), P(:,2), P(:,3), ...
    4, P(:,3), 'filled', ...
    'MarkerFaceAlpha', 0.18, ...
    'MarkerEdgeAlpha', 0.18);

colormap turbo;
cb = colorbar;
cb.Label.String = 'Z height / m';

robot.plot(q_home, ...
    'workspace', [-1 1 -1 1 -0.3 1.2], ...
    'scale', 0.55, ...
    'jointdiam', 1.2, ...
    'notiles', ...
    'noname');

xlabel('X / m','FontSize',13);
ylabel('Y / m','FontSize',13);
zlabel('Z / m','FontSize',13);
title('UR5 Reachable Workspace Cloud Using Craig Modified DH','FontSize',15);

view(135,25);
camlight headlight;
lighting gouraud;

set(gca,'FontSize',12,'LineWidth',1.1);
box on;

%% Figure Export

exportgraphics(gcf, 'UR5_Workspace_Cloud.png', 'Resolution', 600);

%% Interactive Joint-Space Visualization

robot.teach(q_home);
