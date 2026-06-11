clc; clear; close all;

%% UR5 Robot Model Definition

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

robot = SerialLink(L, 'name', 'UR5');

%% Joint Limit Assignment

robot.qlim = deg2rad([
   -180  180
   -120  120
   -150  150
   -180  180
   -120  120
   -360  360
]);

%% Home Configuration

q_home = deg2rad([0 0 0 0 0 0]);

%% Forward Kinematics Evaluation

disp(robot);
T_home = robot.fkine(q_home);
disp('End-effector pose at home:');
disp(T_home);

%% Robot Visualization

figure;
robot.plot(q_home, ...
    'workspace', [-1 1 -1 1 -0.2 1.2], ...
    'scale', 0.6, ...
    'jointdiam', 1.2);

title('UR5 Robot Visualization');

%% Interactive Joint-Space Visualization

robot.teach(q_home);
