%% UR5 Post-Processing Tracking Animation
% This script generates animations and snapshot figures for reduced UR5 nonlinear tracking simulations.
% It uses previously saved planner and execution data and preserves the selected trajectory and tracker configuration.

clear; clc; close all;

%% User-defined simulation and visualization settings
trajectoryMode = 1;

trackerMode = 3;

snapshotFrames = [10 50 90];
saveSelectedAnimationAsGIF = true;
gifFrameSkip = 1;
playbackSpeed = 1.0;

showReferencePath = true;
showTrail = true;
showStartGoal = true;
showCurrentReferenceMarker = true;
showCurrentEndpointMarker = true;

cameraAz = 45;
cameraEl = 30;
axisPadding = 0.20;

ur5LinkColor      = [0.25 0.25 0.25];
ur5JointColor     = [0.35 0.65 0.90];
robotLineWidth    = 4.0;
robotJointSize    = 10;

referencePathColor = [0 0 0];
trailColor         = [0 0.4470 0.7410];
currentRefColor    = [0 0 0];
currentEEColor     = [0 0.4470 0.7410];

gifOutputFolder = '.';

%% Load required planner and execution data
assert(exist('loadrobot','file') == 2 || exist('loadrobot','builtin') == 5, ...
    'Robotics System Toolbox is required: loadrobot not found.');

assert(exist('planner_data.mat','file') == 2, ...
    'planner_data.mat not found. Run mpc.txt first.');
assert(exist('exec_results.mat','file') == 2, ...
    'exec_results.mat not found. Run trackingInNonLinearManipulator_4trackers__ALL_TRAJ.txt first.');

Splan = load('planner_data.mat');
assert(isfield(Splan,'planner_data'), 'planner_data.mat does not contain planner_data.');
planner_data = Splan.planner_data;

Sexec = load('exec_results.mat');
assert(isfield(Sexec,'exec_results'), 'exec_results.mat does not contain exec_results.');
exec_results = Sexec.exec_results;
assert(isfield(exec_results,'trajectories'), ...
    ['exec_results.trajectories not found. Run the ALL_TRAJ tracking script, not the older ', ...
     'single-reference tracking script.']);

%% Define shared model and timing data
robot = loadrobot("universalUR5", ...
    "DataFormat","column", ...
    "Gravity",[0 0 -9.81]);

if isfield(exec_results,'common')
    t         = exec_results.common.t(:).';
    dt        = exec_results.common.dt;
    idx       = exec_results.common.idx(:).'; %#ok<NASGU>
    q_fix_456 = exec_results.common.q_fix_456(:);
    reducedBody = exec_results.common.reducedBody;
else
    t         = planner_data.common.t(:).';
    dt        = planner_data.common.dt;
    idx       = planner_data.common.idx(:).'; %#ok<NASGU>
    q_fix_456 = planner_data.common.q_fix_456(:);
    reducedBody = planner_data.common.reducedBody;
end

methodNames  = {'opt','cubic','quintic','trap'};
methodLabels = {'MPC Optimal Plan','Cubic','Quintic','Trapezoidal'};

trackerNames  = {'ltv_mpc','lqr','pid_ff','pid'};
trackerLabels = {'LTV MPC Tracker','LQR Tracker','PID + Feedforward Tracker','PID Tracker'};

assert(trajectoryMode >= 1 && trajectoryMode <= numel(methodNames), ...
    'Invalid trajectoryMode. Use 1, 2, 3, or 4.');
assert(trackerMode >= 1 && trackerMode <= numel(trackerNames), ...
    'Invalid trackerMode. Use 1, 2, 3, or 4.');

selectedMethod = methodNames{trajectoryMode};
selectedTrajectoryLabel = methodLabels{trajectoryMode};
selectedTracker = trackerNames{trackerMode};
selectedTrackerLabel = trackerLabels{trackerMode};

%% Assemble trajectory data for post-processing
traj = struct();

for i = 1:numel(methodNames)
    md = methodNames{i};
    lbl = methodLabels{i};

    assert(isfield(exec_results.trajectories, md), ...
        'exec_results.trajectories.%s is missing.', md);
    assert(isfield(exec_results.trajectories.(md),'trackers'), ...
        'exec_results.trajectories.%s.trackers is missing.', md);
    assert(isfield(exec_results.trajectories.(md).trackers, selectedTracker), ...
        'Tracker %s is missing for trajectory %s.', selectedTracker, md);

    trInfo = exec_results.trajectories.(md);
    simOut = trInfo.trackers.(selectedTracker);

    if isfield(trInfo,'q_ref')
        q_ref = trInfo.q_ref;
    else
        q_ref = planner_data.(md).q_ref;
    end

    if isfield(simOut,'ref_endpoint')
        ee_ref = simOut.ref_endpoint;
    elseif isfield(trInfo,'ref_endpoint')
        ee_ref = trInfo.ref_endpoint;
    elseif isfield(trInfo,'planner_endpoint_path')
        ee_ref = trInfo.planner_endpoint_path;
    elseif isfield(planner_data.(md),'endpoint_path')
        ee_ref = planner_data.(md).endpoint_path;
    else
        ee_ref = reducedEndpointSeries(robot, reducedBody, q_ref, q_fix_456);
    end

    assert(isfield(simOut,'q'), 'simOut.q is missing for %s/%s.', md, selectedTracker);
    q_anim = simOut.q;

    if isfield(simOut,'endpoint')
        ee_anim = simOut.endpoint;
    else
        ee_anim = reducedEndpointSeries(robot, reducedBody, q_anim, q_fix_456);
    end

    q_anim = ensure3xN(q_anim, ['q_' md '_' selectedTracker]);
    q_ref  = ensure3xN(q_ref,  ['q_ref_' md]);
    ee_anim = ensure3xN(ee_anim, ['ee_' md '_' selectedTracker]);
    ee_ref  = ensure3xN(ee_ref,  ['ee_ref_' md]);

    Ncommon = min([size(q_anim,2), size(q_ref,2), size(ee_anim,2), size(ee_ref,2), numel(t)]);
    q_anim  = q_anim(:,1:Ncommon);
    q_ref   = q_ref(:,1:Ncommon);
    ee_anim = ee_anim(:,1:Ncommon);
    ee_ref  = ee_ref(:,1:Ncommon);
    t_now   = t(1:Ncommon);

    plotBodies = pickReducedUR5PlotBodies(robot, reducedBody);
    chain = reducedChainSeries(robot, plotBodies, q_anim, q_fix_456);

    traj.(md).label        = lbl;
    traj.(md).trackerLabel = selectedTrackerLabel;
    traj.(md).trackerName  = selectedTracker;
    traj.(md).sourceLabel  = 'executed nonlinear simulation';
    traj.(md).q            = q_anim;
    traj.(md).q_ref        = q_ref;
    traj.(md).ee           = ee_anim;
    traj.(md).ee_ref       = ee_ref;
    traj.(md).chain        = chain;
    traj.(md).t            = t_now;
    traj.(md).N            = Ncommon;
    traj.(md).plotBodies   = plotBodies;

    if isfield(simOut,'metrics')
        traj.(md).metrics = simOut.metrics;
    end
end

%% Validate requested snapshot frames
Nmin = inf;
for i = 1:numel(methodNames)
    Nmin = min(Nmin, traj.(methodNames{i}).N);
end

assert(all(snapshotFrames >= 1) && all(snapshotFrames <= Nmin), ...
    'snapshotFrames must be within 1:%d for the selected tracker results.', Nmin);

%% Compute common plotting limits
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

%% Generate animation for the selected trajectory and tracker
selectedData = traj.(selectedMethod);
fig1 = figure('Name',sprintf('Figure 1 - %s - %s Animation', ...
    selectedTrajectoryLabel, selectedTrackerLabel), 'Color','w');
ax1 = axes('Parent', fig1);

if saveSelectedAnimationAsGIF
    safeTrajName = sanitizeFileTag(selectedMethod);
    safeTrackerName = sanitizeFileTag(selectedTracker);
    gifFileName = sprintf('UR5_tracking_%s_%s.gif', safeTrajName, safeTrackerName);
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

%% Generate snapshot panel for the selected trajectory and tracker
fig2 = createSnapshotFigure1x3( ...
    2, selectedData, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor);

%% Generate comparative snapshot panel for all trajectories
fig3 = figure('Name',sprintf('Figure 3 - All Trajectories - %s Snapshot Panel', selectedTrackerLabel), ...
    'Color','w', 'Units','pixels', 'Position',[100 80 1500 1200]);

for r = 1:numel(methodNames)
    md = methodNames{r};
    data = traj.(md);

    for c = 1:numel(snapshotFrames)
        k = snapshotFrames(c);
        ax = subplot(4,3,(r-1)*3 + c, 'Parent', fig3);
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
            data.label, data.trackerLabel, k, data.N, data.t(k)), 'FontWeight','bold');
        set(ax,'FontSize',9);
    end
end

sgtitle(fig3, sprintf('Reduced UR5 Nonlinear Tracking Snapshots - %s', selectedTrackerLabel), ...
    'FontWeight','bold');

%% Report completion status
fprintf('\n============================================================\n');
fprintf('UR5 post-process tracking animation completed.\n');
fprintf('Selected trajectory mode   : %d -> %s (%s)\n', trajectoryMode, selectedTrajectoryLabel, selectedMethod);
fprintf('Selected tracker mode      : %d -> %s (%s)\n', trackerMode, selectedTrackerLabel, selectedTracker);
fprintf('Reduced endpoint body      : %s\n', reducedBody);
if saveSelectedAnimationAsGIF
    fprintf('GIF saved to               : %s\n', gifFilePath);
else
    fprintf('GIF saving                 : disabled\n');
end
fprintf('Figures created            : 1 animation, 2 selected snapshots, 3 all-trajectory snapshots\n');
fprintf('============================================================\n\n');

%% Local functions

% Animate one trajectory-tracker combination and optionally export it as a GIF.
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

    title(ax, sprintf('%s | %s | Frame %d/%d | t = %.3f s', ...
        data.label, data.trackerLabel, k, N, t(k)), 'FontWeight','bold');

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

% Create a three-frame snapshot figure for one trajectory-tracker combination.
function fig = createSnapshotFigure1x3( ...
    figNumber, data, snapshotFrames, ...
    xLim, yLim, zLim, ...
    cameraAz, cameraEl, ...
    showReferencePath, showTrail, showStartGoal, ...
    showCurrentReferenceMarker, showCurrentEndpointMarker, ...
    ur5LinkColor, ur5JointColor, robotLineWidth, robotJointSize, ...
    referencePathColor, trailColor, currentRefColor, currentEEColor)

fig = figure('Name',sprintf('Figure %d - %s - %s 1x3 Snapshots', figNumber, data.label, data.trackerLabel), ...
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

    title(ax, sprintf('%s | %s | Frame %d/%d | t = %.3f s', ...
        data.label, data.trackerLabel, k, data.N, data.t(k)), 'FontWeight','bold');
end

sgtitle(fig, sprintf('%s - %s - 1x3 Snapshot Panel', data.label, data.trackerLabel), ...
    'FontWeight','bold');
end

% Draw the reduced UR5 configuration, reference path, and executed path at one time step.
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

% Configure axis limits, grid, labels, aspect ratio, and camera view.
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

% Select UR5 body frames used to visualize the reduced kinematic chain.
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

% Compute reduced kinematic chain coordinates over the full trajectory.
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

% Compute reduced end-effector positions over the full trajectory.
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

% Validate that a trajectory array has size 3-by-N.
function X = ensure3xN(X, varName)
validateattributes(X, {'double'}, {'2d','nrows',3,'finite','real'}, mfilename, varName);
end

% Compute frame delay from trajectory timing and playback speed.
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

% Return one of two inputs according to a logical condition.
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

% Convert a string to a file-name-safe tag.
function tag = sanitizeFileTag(str)
tag = lower(char(str));
tag = regexprep(tag, '[^a-z0-9]+', '_');
tag = regexprep(tag, '^_+|_+$', '');
if isempty(tag)
    tag = 'selection';
end
end
