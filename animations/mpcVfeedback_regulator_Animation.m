%% UR5 Nonlinear Tracking Post-Processing
% This script visualizes nonlinear tracking results for selected planned
% trajectories and tracking controllers using the reduced UR5 model.

clear; clc; close all;

%% Simulation and Plotting Settings
trajectoryMode = 3;

trackerMode = 2;

snapshotFrames = [10 50 90];
saveSelectedAnimationAsGIF = true;
gifFrameSkip = 1;
playbackSpeed = 1.0;

showReferencePath = true;
showTrail = true;
showStartGoal = true;
showCurrentReferenceMarker = true;
showCurrentEndpointMarker = true;

% Camera and axis configuration
cameraAz = 45;
cameraEl = 30;
axisPadding = 0.20;

% UR5 visual appearance
ur5LinkColor      = [0.25 0.25 0.25];
ur5JointColor     = [0.35 0.65 0.90];
robotLineWidth    = 4.0;
robotJointSize    = 10;

% Trajectory visualization styling
referencePathColor = [0 0 0];
trailColor         = [0 0.4470 0.7410];
currentRefColor    = [0 0 0];
currentEEColor     = [0 0.4470 0.7410];

% Animation output directory
gifOutputFolder = '.';

%% Load Simulation Data
assert(exist('loadrobot','file') == 2 || exist('loadrobot','builtin') == 5, ...
    'Robotics System Toolbox is required: loadrobot not found.');

assert(exist('planner_data.mat','file') == 2, ...
    'planner_data.mat not found. Run mpc1_pid_lqr_full.txt first.');

assert(exist('exec_results.mat','file') == 2, ...
    ['exec_results.mat not found. Run trackingInNonLinearManipulator_MPC_LQR_PID_Anim.txt ', ...
     'or the modified tracking script first.']);

Splan = load('planner_data.mat');
assert(isfield(Splan,'planner_data'), 'planner_data.mat does not contain planner_data.');
planner_data = Splan.planner_data;

Sexec = load('exec_results.mat');
assert(isfield(Sexec,'exec_results'), 'exec_results.mat does not contain exec_results.');
exec_results = Sexec.exec_results;

assert(isfield(exec_results,'trajectories'), ...
    ['exec_results.trajectories not found. The tracking script must save ', ...
     'exec_results.trajectories.<trajectory>.trackers.<tracker>.']);

%% Load Common Model and Time Data
robot = loadrobot("universalUR5", ...
    "DataFormat","column", ...
    "Gravity",[0 0 -9.81]);

if isfield(exec_results,'common')
    t         = exec_results.common.t(:).';
    dt        = exec_results.common.dt;
    q0        = exec_results.common.q0(:);
    qg        = exec_results.common.qg(:);
    idx       = exec_results.common.idx(:).';
    q_fix_456 = exec_results.common.q_fix_456(:);
    reducedBody = exec_results.common.reducedBody;
else
    t         = planner_data.common.t(:).';
    dt        = planner_data.common.dt;
    q0        = planner_data.common.q0(:);
    qg        = planner_data.common.qg(:);
    idx       = planner_data.common.idx(:).';
    q_fix_456 = planner_data.common.q_fix_456(:);
    reducedBody = planner_data.common.reducedBody;
end

methodNames  = {'opt','pid','lqr'};
methodLabels = {'LTV-MPC Plan','PID Regulator Plan','LQR Regulator Plan'};

trackerNames  = {'ltv_mpc','lqr','pid_ff','pid'};
trackerLabels = {'MPC Tracker','LQR Tracker','PID + Feedforward Tracker','PID Tracker'};

assert(trajectoryMode >= 1 && trajectoryMode <= numel(methodNames), ...
    'Invalid trajectoryMode. Use trajectoryMode = 1, 2, or 3.');
assert(trackerMode >= 1 && trackerMode <= numel(trackerNames), ...
    'Invalid trackerMode. Use trackerMode = 1, 2, 3, or 4.');

selectedMethod = methodNames{trajectoryMode};
selectedTrajectoryLabel = methodLabels{trajectoryMode};
selectedTracker = trackerNames{trackerMode};
selectedTrackerLabel = trackerLabels{trackerMode};

%% Assemble Trajectory Data
% The reference path is the corresponding planned trajectory from planner_data

traj = struct();
plotBodies = pickReducedUR5PlotBodies(robot, reducedBody);

for i = 1:numel(methodNames)
    md = methodNames{i};
    baseTrajectoryLabel = methodLabels{i};

    assert(isfield(exec_results.trajectories, md), ...
        'exec_results.trajectories.%s is missing. Re-run the modified tracking script.', md);

    trInfo = exec_results.trajectories.(md);
    assert(isfield(trInfo,'trackers'), ...
        'exec_results.trajectories.%s.trackers is missing.', md);
    assert(isfield(trInfo.trackers, selectedTracker), ...
        'Tracker "%s" is missing for trajectory "%s". Re-run tracking with all trajectory-tracker combinations.', ...
        selectedTracker, md);

    simOut = trInfo.trackers.(selectedTracker);

    if isfield(trInfo,'q_ref')
        q_ref = trInfo.q_ref;
    elseif isfield(simOut,'q_ref')
        q_ref = simOut.q_ref;
    elseif isfield(planner_data,md) && isfield(planner_data.(md),'q_ref')
        q_ref = planner_data.(md).q_ref;
    else
        error('Could not find reference q_ref for trajectory %s.', md);
    end

    if isfield(simOut,'ref_endpoint')
        ee_ref = simOut.ref_endpoint;
    elseif isfield(trInfo,'ref_endpoint')
        ee_ref = trInfo.ref_endpoint;
    elseif isfield(trInfo,'planner_endpoint_path')
        ee_ref = trInfo.planner_endpoint_path;
    elseif isfield(planner_data,md) && isfield(planner_data.(md),'endpoint_path')
        ee_ref = planner_data.(md).endpoint_path;
    else
        ee_ref = reducedEndpointSeries(robot, reducedBody, q_ref, q_fix_456);
    end

    assert(isfield(simOut,'q'), ...
        'simOut.q is missing for trajectory "%s" and tracker "%s".', md, selectedTracker);
    q_anim = simOut.q;

    if isfield(simOut,'endpoint')
        ee_anim = simOut.endpoint;
    else
        ee_anim = reducedEndpointSeries(robot, reducedBody, q_anim, q_fix_456);
    end

    q_anim  = ensure3xN(q_anim,  ['q_' md '_' selectedTracker]);
    q_ref   = ensure3xN(q_ref,   ['q_ref_' md]);
    ee_anim = ensure3xN(ee_anim, ['ee_' md '_' selectedTracker]);
    ee_ref  = ensure3xN(ee_ref,  ['ee_ref_' md]);

    Ncommon = min([size(q_anim,2), size(q_ref,2), size(ee_anim,2), size(ee_ref,2), numel(t)]);
    q_anim  = q_anim(:,1:Ncommon);
    q_ref   = q_ref(:,1:Ncommon);
    ee_anim = ee_anim(:,1:Ncommon);
    ee_ref  = ee_ref(:,1:Ncommon);
    t_now   = t(1:Ncommon);

    chain = reducedChainSeries(robot, plotBodies, q_anim, q_fix_456);

    traj.(md).baseLabel    = baseTrajectoryLabel;
    traj.(md).trackerLabel = selectedTrackerLabel;
    traj.(md).label        = sprintf('%s tracked by %s', baseTrajectoryLabel, selectedTrackerLabel);
    traj.(md).sourceLabel  = 'nonlinear executed';
    traj.(md).q            = q_anim;
    traj.(md).q_ref        = q_ref;
    traj.(md).ee           = ee_anim;
    traj.(md).ee_ref       = ee_ref;
    traj.(md).chain        = chain;
    traj.(md).t            = t_now;
    traj.(md).N            = Ncommon;
    traj.(md).plotBodies   = plotBodies;
    traj.(md).methodName   = md;
    traj.(md).trackerName  = selectedTracker;

    if isfield(simOut,'metrics')
        traj.(md).metrics = simOut.metrics;
    end
    if isfield(simOut,'err_norm')
        traj.(md).err_norm = simOut.err_norm(:,1:Ncommon);
    end
end

%% Validate Snapshot Frames
Nmin = inf;
for i = 1:numel(methodNames)
    Nmin = min(Nmin, traj.(methodNames{i}).N);
end

assert(all(snapshotFrames >= 1) && all(snapshotFrames <= Nmin), ...
    'snapshotFrames must be within 1:%d.', Nmin);

%% Compute Common Plot Limits
allPts = [];
for i = 1:numel(methodNames)
    md = methodNames{i};
    allPts = [allPts, traj.(md).ee, traj.(md).ee_ref, reshape(traj.(md).chain, 3, [])];
end

xyzMin = min(allPts, [], 2);
xyzMax = max(allPts, [], 2);

xLim = [xyzMin(1)-axisPadding, xyzMax(1)+axisPadding];
yLim = [xyzMin(2)-axisPadding, xyzMax(2)+axisPadding];
zLim = [min(0, xyzMin(3)-axisPadding), xyzMax(3)+axisPadding];

if range(xLim) < 0.5, xLim = mean(xLim) + [-0.5 0.5]; end
if range(yLim) < 0.5, yLim = mean(yLim) + [-0.5 0.5]; end
if range(zLim) < 0.5, zLim = mean(zLim) + [-0.5 0.5]; end

%% Generate Selected Animation
selectedData = traj.(selectedMethod);

fig1 = figure('Name',sprintf('Figure 1 - %s / %s Animation', ...
    selectedTrajectoryLabel, selectedTrackerLabel), 'Color','w');
ax1 = axes('Parent', fig1);

if saveSelectedAnimationAsGIF
    safeTraj = sanitizeFileTag(selectedMethod);
    safeTracker = sanitizeFileTag(selectedTracker);
    gifFileName = sprintf('UR5_%s_plan_tracked_by_%s.gif', safeTraj, safeTracker);
    gifFilePath = fullfile(gifOutputFolder, gifFileName);
    if exist(gifFilePath,'file') == 2
        delete(gifFilePath);
    end
else
    gifFilePath = '';
end

animateSingleMethod( ...
    ax1, selectedData, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor, ...
    gifFilePath, gifFrameSkip, playbackSpeed);

%% Generate Selected Snapshot Panel
fig2 = createSnapshotFigure1x3( ...
    2, selectedData, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor);

%% Figure 3: All Trajectories, Selected Tracker
fig3 = figure('Name',sprintf('Figure 3 - All Plans tracked by %s', selectedTrackerLabel), ...
    'Color','w', 'Units','pixels', 'Position',[100 80 1500 900]);

for r = 1:numel(methodNames)
    md = methodNames{r};
    data = traj.(md);

    for c = 1:numel(snapshotFrames)
        k = snapshotFrames(c);
        ax = subplot(numel(methodNames), numel(snapshotFrames), ...
            (r-1)*numel(snapshotFrames) + c, 'Parent', fig3);

        drawSnapshotOnAxes( ...
            ax, data, k, ...
            xLim, yLim, zLim, ...
            cameraAz, cameraEl, ...
            showReferencePath, showTrail, showStartGoal, ...
            showCurrentReferenceMarker, showCurrentEndpointMarker, ...
            ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
            referencePathColor, trailColor, currentRefColor, currentEEColor, ...
            true);

        title(ax, sprintf('%s | %s | Frame %d/%d | t = %.3f s', ...
            data.baseLabel, selectedTrackerLabel, k, data.N, data.t(k)), ...
            'FontWeight','bold');
        set(ax,'FontSize',9);
    end
end

sgtitle(fig3, sprintf('Reduced UR5 nonlinear tracking snapshots - %s', selectedTrackerLabel), ...
    'FontWeight','bold');

%% Display Execution Summary
fprintf('\n============================================================\n');
fprintf('UR5 post-process tracking animation completed.\n');
fprintf('Selected trajectoryMode    : %d -> %s (%s)\n', trajectoryMode, selectedTrajectoryLabel, selectedMethod);
fprintf('Selected trackerMode       : %d -> %s (%s)\n', trackerMode, selectedTrackerLabel, selectedTracker);
fprintf('Reduced endpoint body      : %s\n', reducedBody);
if isfield(selectedData,'metrics')
    M = selectedData.metrics;
    if isfield(M,'rms_position_error')
        fprintf('RMS position tracking err  : %.6e rad\n', M.rms_position_error);
    end
    if isfield(M,'final_position_error')
        fprintf('Final position tracking err: %.6e rad\n', M.final_position_error);
    end
end
if saveSelectedAnimationAsGIF
    fprintf('GIF saved to               : %s\n', gifFilePath);
else
    fprintf('GIF saving                 : disabled\n');
end
fprintf('Figures created            : 1 animation, 2 selected snapshots, 3 all-plan panel\n');
fprintf('============================================================\n\n');


%% 
%% LOCAL FUNCTIONS
%% 

% Animate the selected trajectory and optionally export the animation as a GIF.
function animateSingleMethod( ...
    ax, data, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor, ...
    gifFilePath, gifFrameSkip, playbackSpeed)

N = data.N;
t = data.t;
ee = data.ee;
ee_ref = data.ee_ref;
chain = data.chain;

setupAxes(ax, xLim, yLim, zLim, cameraAz, cameraEl);

startEE = ee_ref(:,1);
goalEE  = ee_ref(:,end);

for k = 1:N
    cla(ax);
    hold(ax,'on');

    if showReferencePath
        plot3(ax, ee_ref(1,:), ee_ref(2,:), ee_ref(3,:), '--', ...
            'Color', referencePathColor, 'LineWidth', 1.2);
    end

    if showTrail
        plot3(ax, ee(1,1:k), ee(2,1:k), ee(3,1:k), '-', ...
            'Color', trailColor, 'LineWidth', 1.8);
    end

    if showStartGoal
        plot3(ax, startEE(1), startEE(2), startEE(3), 'rs', ...
            'MarkerSize', 8, 'LineWidth', 1.4);
        plot3(ax, goalEE(1), goalEE(2), goalEE(3), 'gp', ...
            'MarkerSize', 9, 'LineWidth', 1.4);
    end

    P = chain(:,:,k);

    plot3(ax, P(1,:), P(2,:), P(3,:), '-', ...
        'Color', ur5LinkColor, 'LineWidth', robotLineWidth);

    plot3(ax, P(1,:), P(2,:), P(3,:), 'o', ...
        'LineStyle', 'none', ...
        'MarkerSize', robotJointSize, ...
        'MarkerEdgeColor', ur5JointColor, ...
        'MarkerFaceColor', ur5JointColor);

    if showCurrentEndpointMarker
        plot3(ax, ee(1,k), ee(2,k), ee(3,k), 'o', ...
            'MarkerSize', 7, 'LineWidth', 1.3, ...
            'Color', currentEEColor, ...
            'MarkerFaceColor', currentEEColor);
    end

    if showCurrentReferenceMarker
        plot3(ax, ee_ref(1,k), ee_ref(2,k), ee_ref(3,k), 'o', ...
            'MarkerSize', 6, 'LineWidth', 1.2, ...
            'Color', currentRefColor);
    end

    title(ax, sprintf('%s (%s) | Frame %d/%d | t = %.3f s', ...
        data.label, data.sourceLabel, k, N, t(k)), 'FontWeight','bold');

    drawnow;

    if ~isempty(gifFilePath) && mod(k-1, gifFrameSkip) == 0
        frame = getframe(ax.Parent);
        [imind, cm] = rgb2ind(frame2im(frame), 256);
        if k == 1
            imwrite(imind, cm, gifFilePath, 'gif', 'Loopcount', inf, ...
                'DelayTime', computeDelay(t, playbackSpeed, k));
        else
            imwrite(imind, cm, gifFilePath, 'gif', 'WriteMode', 'append', ...
                'DelayTime', computeDelay(t, playbackSpeed, k));
        end
    end

    pause(computeDelay(t, playbackSpeed, k));
end
end

% Create a three-frame snapshot panel for the selected trajectory.
function fig = createSnapshotFigure1x3( ...
    figNumber, data, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor)

fig = figure('Name',sprintf('Figure %d - %s 1x3 Snapshots', figNumber, data.label), ...
    'Color','w', 'Units','pixels', 'Position',[120 140 1450 420]);

for c = 1:numel(snapshotFrames)
    k = snapshotFrames(c);
    ax = subplot(1,3,c,'Parent',fig);
    drawSnapshotOnAxes( ...
        ax, data, k, ...
        xLim, yLim, zLim, ...
        cameraAz, cameraEl, ...
        showReferencePath, showTrail, showStartGoal, ...
        showCurrentReferenceMarker, showCurrentEndpointMarker, ...
        ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
        referencePathColor, trailColor, currentRefColor, currentEEColor, ...
        true);

    title(ax, sprintf('%s | Frame %d/%d | t = %.3f s', ...
        data.label, k, data.N, data.t(k)), 'FontWeight','bold');
end

sgtitle(fig, sprintf('%s (%s) - 1x3 Snapshot Panel', data.label, data.sourceLabel), ...
    'FontWeight','bold');
end

% Draw the UR5 configuration, reference path, and executed trajectory at one frame.
function drawSnapshotOnAxes( ...
    ax, data, k, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor, ...
    showLabels)

setupAxes(ax, xLim, yLim, zLim, cameraAz, cameraEl);
hold(ax,'on');

ee = data.ee;
ee_ref = data.ee_ref;
chain = data.chain;

startEE = ee_ref(:,1);
goalEE  = ee_ref(:,end);

if showReferencePath
    plot3(ax, ee_ref(1,:), ee_ref(2,:), ee_ref(3,:), '--', ...
        'Color', referencePathColor, 'LineWidth', 1.2);
end

if showTrail
    plot3(ax, ee(1,1:k), ee(2,1:k), ee(3,1:k), '-', ...
        'Color', trailColor, 'LineWidth', 1.8);
end

if showStartGoal
    plot3(ax, startEE(1), startEE(2), startEE(3), 'rs', ...
        'MarkerSize', 8, 'LineWidth', 1.4);
    plot3(ax, goalEE(1), goalEE(2), goalEE(3), 'gp', ...
        'MarkerSize', 9, 'LineWidth', 1.4);
end

P = chain(:,:,k);

plot3(ax, P(1,:), P(2,:), P(3,:), '-', ...
    'Color', ur5LinkColor, 'LineWidth', robotLineWidth);

plot3(ax, P(1,:), P(2,:), P(3,:), 'o', ...
    'LineStyle', 'none', ...
    'MarkerSize', robotJointSize, ...
    'MarkerEdgeColor', ur5JointColor, ...
    'MarkerFaceColor', ur5JointColor);

if showCurrentEndpointMarker
    plot3(ax, ee(1,k), ee(2,k), ee(3,k), 'o', ...
        'MarkerSize', 7, 'LineWidth', 1.3, ...
        'Color', currentEEColor, ...
        'MarkerFaceColor', currentEEColor);
end

if showCurrentReferenceMarker
    plot3(ax, ee_ref(1,k), ee_ref(2,k), ee_ref(3,k), 'o', ...
        'MarkerSize', 6, 'LineWidth', 1.2, ...
        'Color', currentRefColor);
end

if showLabels
    xlabel(ax,'X [m]');
    ylabel(ax,'Y [m]');
    zlabel(ax,'Z [m]');
end
end

% Configure three-dimensional plotting axes.
function setupAxes(ax, xLim, yLim, zLim, cameraAz, cameraEl)
cla(ax);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');
axis(ax,'equal');
xlabel(ax,'X [m]');
ylabel(ax,'Y [m]');
zlabel(ax,'Z [m]');
xlim(ax, xLim);
ylim(ax, yLim);
zlim(ax, zLim);
view(ax, cameraAz, cameraEl);
end

% Select the UR5 body frames used for reduced-chain visualization.
function plotBodies = pickReducedUR5PlotBodies(robot, reducedBody)
preferred = {'base_link','shoulder_link','upper_arm_link','forearm_link','wrist_1_link'};
plotBodies = {};
for i = 1:numel(preferred)
    if any(strcmp(robot.BodyNames, preferred{i}))
        plotBodies{end+1} = preferred{i};
    end
    if strcmp(preferred{i}, reducedBody)
        break;
    end
end

if isempty(plotBodies) || ~strcmp(plotBodies{end}, reducedBody)
    plotBodies = {'base_link', reducedBody};
    plotBodies = plotBodies(ismember(plotBodies, robot.BodyNames));
end

assert(~isempty(plotBodies), 'Could not build reduced UR5 plot body list.');
end

% Compute Cartesian coordinates of the reduced UR5 kinematic chain over time.
function chain = reducedChainSeries(robot, plotBodies, q_traj, q_fix_456)
q_traj = ensure3xN(q_traj, 'q_traj');
numBodies = numel(plotBodies);
numSamples = size(q_traj,2);
chain = zeros(3, numBodies, numSamples);

for k = 1:numSamples
    q_full = [q_traj(:,k); q_fix_456(:)];
    for b = 1:numBodies
        T = getTransform(robot, q_full, plotBodies{b});
        chain(:,b,k) = T(1:3,4);
    end
end
end

% Compute the reduced end-effector position trajectory from joint coordinates.
function ee = reducedEndpointSeries(robot, bodyName, q_traj, q_fix_456)
q_traj = ensure3xN(q_traj, 'q_traj');
numSamples = size(q_traj,2);
ee = zeros(3, numSamples);

for k = 1:numSamples
    q_full = [q_traj(:,k); q_fix_456(:)];
    T = getTransform(robot, q_full, bodyName);
    ee(:,k) = T(1:3,4);
end
end

% Validate that trajectory arrays are real, finite 3-by-N matrices.
function X = ensure3xN(X, varName)
validateattributes(X, {'double'}, {'2d','nrows',3,'finite','real'}, mfilename, varName);
end

% Compute playback delay from simulation time and speed scaling.
function d = computeDelay(t, playbackSpeed, k)
if numel(t) >= 2
    if k < numel(t)
        nominalDt = t(k+1) - t(k);
    else
        nominalDt = mean(diff(t));
    end
else
    nominalDt = 0.02;
end
if ~isfinite(nominalDt) || nominalDt <= 0
    nominalDt = 0.02;
end
if ~isfinite(playbackSpeed) || playbackSpeed <= 0
    playbackSpeed = 1.0;
end
d = nominalDt / playbackSpeed;
end

% Return one of two values based on a logical condition.
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

% Convert a string into a safe filename tag.
function tag = sanitizeFileTag(tagIn)
tag = char(string(tagIn));
tag = regexprep(tag, '[^A-Za-z0-9_]+', '_');
tag = regexprep(tag, '_+', '_');
tag = lower(strtrim(tag));
if isempty(tag)
    tag = 'selected';
end
end
