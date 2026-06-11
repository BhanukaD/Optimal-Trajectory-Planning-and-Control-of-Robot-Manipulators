%% Nonlinear Tracking Simulation and Controller Evaluation

clear; clc; close all;

%% User settings
plannerFile = 'planner_data.mat';
outputFile  = 'exec_results.mat';

% Trajectories to execute. These names match mpc1_pid_lqr_full.txt.
trajectoryNames  = {'opt','pid','lqr'};
trajectoryLabels = {'MPC planned optimal trajectory', ...
                    'PID-regulator planned trajectory', ...
                    'LQR-regulator planned trajectory'};

% Controllers to execute on each selected reference trajectory.
controllerNames = {'pid','pid_ff','lqr','ltv_mpc'};
controllerLabels = {'PID only','PID + feedforward','Midpoint LTI-LQR','Receding-horizon LTV MPC'};

% Optional compact plots from this tracking script.
% The animation is generated later by the separate post-process animation code.
enableSummaryPlots = true;
plotTrajectoryName = 'opt';

% PID gains
Kp = diag([120 120 120]);
Ki = diag([ 10  10  10]);
Kd = diag([ 30  30  30]);
ei_limit = [0.50; 0.50; 0.50];   % [rad*s]

% Midpoint LTI-LQR weights for tracking-error regulation
% State: x_tilde = [q-q_ref; qd-qd_ref]
Qlqr = blkdiag(60*eye(3), 10*eye(3));
Rlqr = 1.0*eye(3);
Plqr_terminal = [];          % not used by midpoint LTI-LQR; kept for saved-data compatibility
linearization_eps_x = 1e-6;
linearization_eps_u = 1e-5;
lqr_uses_nominal_feedforward = true;   % standard tracking: delta_tau = u_nom - K*(x-x_ref)

% Receding-horizon online LTV MPC tracking controller settings
QmpcTrack = blkdiag(1200*eye(3), 140*eye(3));
QfmpcTrack = 80*QmpcTrack;
RduMpcTrack = 0.002*eye(3);
mpcTrackPredictionSteps = 25;
mpcTrackControlHorizon = 1;
mpcTrackUseQuadprog = true;
mpcTrackUseInputBounds = true;
mpcTrackUseIncrementBounds = true;
mpcTrackDeltaUMax = [60; 60; 60];
mpcTrackEpsX = 1e-6;
mpcTrackEpsU = 1e-6;

% Gravity compensation convention used by the reduced plant
% tau_abs = G(q) + delta_tau_cmd, therefore the reduced plant is
% D(q)qdd = delta_tau_cmd - C(q,qd).
gravity_compensation_enabled = true;

% Saturation settings
use_delta_tau_saturation = true;
delta_tau_limit = [60; 60; 60];  % [Nm]

use_tau_abs_saturation = false;
tau_limit = [150; 150; 150];     % [Nm]

% Executed tracking cost weights
Qexec = blkdiag(200*eye(3), 20*eye(3));
Rexec = 1e-2*eye(3);

% Settling-time settings
compute_settling_time_flag = true;
settling_threshold = 0.02;       % threshold on ||q_ref-q||_2 [rad]
settling_hold_samples = 5;

% Plot settings
lineWidthMain = 1.8;
fontSizeAxes = 11;

%% Load planner data and validate
assert(exist(plannerFile,'file') == 2, 'Required file "%s" was not found.', plannerFile);
S = load(plannerFile);
assert(isfield(S,'planner_data'), 'The file "%s" does not contain planner_data.', plannerFile);
planner_data = S.planner_data;

requiredCommonFields = {'t','dt','Tf','N','q0','qg','qd0','qdg','idx','q_fix_456','reducedBody'};
for i = 1:numel(requiredCommonFields)
    assert(isfield(planner_data,'common') && isfield(planner_data.common, requiredCommonFields{i}), ...
        'planner_data.common.%s is missing.', requiredCommonFields{i});
end

t            = planner_data.common.t(:);
dt           = planner_data.common.dt;
Tf           = planner_data.common.Tf;
N            = planner_data.common.N;
q0           = planner_data.common.q0(:);
qg           = planner_data.common.qg(:);
qd0          = planner_data.common.qd0(:);
qdg          = planner_data.common.qdg(:);
idx          = planner_data.common.idx(:).';
q_fix_456    = planner_data.common.q_fix_456(:);
reducedBody  = planner_data.common.reducedBody;

assert(numel(t) == N+1, 'planner_data.common.t must contain N+1 samples.');
assert(numel(idx) == 3, 'idx must contain exactly 3 active joints.');

for itraj = 1:numel(trajectoryNames)
    trName = trajectoryNames{itraj};
    assert(isfield(planner_data,trName), ...
        'planner_data.%s is missing. Run mpc1_pid_lqr_full.txt first.', trName);
    assert(isfield(planner_data.(trName),'q_ref'), ...
        'planner_data.%s.q_ref is missing.', trName);
    assert(isfield(planner_data.(trName),'qd_ref'), ...
        'planner_data.%s.qd_ref is missing.', trName);
    assert(isfield(planner_data.(trName),'qdd_ref'), ...
        'planner_data.%s.qdd_ref is missing.', trName);
    checkTrajectoryDimensions(planner_data.(trName).q_ref, ...
        planner_data.(trName).qd_ref, planner_data.(trName).qdd_ref, N+1, trName);
end

%% Load UR5 model
assert(exist('loadrobot','file') == 2 || exist('loadrobot','builtin') == 5, ...
    'Robotics System Toolbox is required: loadrobot not found.');
robot = loadrobot("universalUR5", "DataFormat","column", "Gravity",[0 0 -9.81]);

%% Common result container
exec_results = struct();
exec_results.common = struct();
exec_results.common.source_planner_file = plannerFile;
exec_results.common.description = ['MPC/PID/LQR planner trajectories executed on nonlinear reduced UR5 using PID, ', ...
    'PID+feedforward, midpoint LTI-LQR, and receding-horizon online LTV MPC tracking.'];
exec_results.common.trajectory_names = trajectoryNames;
exec_results.common.trajectory_labels = trajectoryLabels;
exec_results.common.controller_names = controllerNames;
exec_results.common.controller_labels = controllerLabels;
exec_results.common.animation_trajectory_mode_map.mode_1 = 'exec_results.trajectories.opt';
exec_results.common.animation_trajectory_mode_map.mode_2 = 'exec_results.trajectories.pid';
exec_results.common.animation_trajectory_mode_map.mode_3 = 'exec_results.trajectories.lqr';
exec_results.common.animation_tracker_mode_map.mode_1 = 'trackers.pid';
exec_results.common.animation_tracker_mode_map.mode_2 = 'trackers.pid_ff';
exec_results.common.animation_tracker_mode_map.mode_3 = 'trackers.lqr';
exec_results.common.animation_tracker_mode_map.mode_4 = 'trackers.ltv_mpc';
exec_results.common.legacy_animation_field = 'exec_results.opt is copied from exec_results.trajectories.opt.trackers.pid_ff';
exec_results.common.Kp = Kp;
exec_results.common.Ki = Ki;
exec_results.common.Kd = Kd;
exec_results.common.ei_limit = ei_limit;
exec_results.common.Qexec = Qexec;
exec_results.common.Rexec = Rexec;
exec_results.common.Qlqr = Qlqr;
exec_results.common.Rlqr = Rlqr;
exec_results.common.Plqr_terminal = Plqr_terminal;
exec_results.common.lqr_formula = 'Midpoint LTI-LQR: A,B from numerical linearization at trajectory midpoint; solve CARE A''P+P*A-P*B*inv(R)*B''*P+Q=0; K=inv(R)*B''*P; tracking law delta_tau = u_nom - K*(x-x_ref).';
exec_results.common.gravity_compensation_enabled = gravity_compensation_enabled;
exec_results.common.use_delta_tau_saturation = use_delta_tau_saturation;
exec_results.common.delta_tau_limit = delta_tau_limit;
exec_results.common.use_tau_abs_saturation = use_tau_abs_saturation;
exec_results.common.tau_limit = tau_limit;
exec_results.common.dt = dt;
exec_results.common.t = t;
exec_results.common.Tf = Tf;
exec_results.common.N = N;
exec_results.common.q0 = q0;
exec_results.common.qg = qg;
exec_results.common.qd0 = qd0;
exec_results.common.qdg = qdg;
exec_results.common.idx = idx;
exec_results.common.q_fix_456 = q_fix_456;
exec_results.common.reducedBody = reducedBody;
exec_results.common.settling_threshold = settling_threshold;
exec_results.common.settling_hold_samples = settling_hold_samples;
exec_results.common.lqr_uses_nominal_feedforward = lqr_uses_nominal_feedforward;
exec_results.common.QmpcTrack = QmpcTrack;
exec_results.common.QfmpcTrack = QfmpcTrack;
exec_results.common.RduMpcTrack = RduMpcTrack;
exec_results.common.mpcTrackPredictionSteps = mpcTrackPredictionSteps;
exec_results.common.mpcTrackControlHorizon = mpcTrackControlHorizon;
exec_results.common.mpcTrackUseInputBounds = mpcTrackUseInputBounds;
exec_results.common.mpcTrackUseIncrementBounds = mpcTrackUseIncrementBounds;
exec_results.common.mpcTrackDeltaUMax = mpcTrackDeltaUMax;
exec_results.common.mpc_tracking_formula = 'Receding-horizon LTV MPC tracking around each selected nominal plan: e=[q-qref;qd-qdref], c=delta_tau_cmd-u_nom, zeta=[e;c(k-1)], c(k)=c(k-1)+Delta c(k), delta_tau_cmd=u_nom+c(k).';

mpcTracker = designLtvMpcTracker(QmpcTrack, QfmpcTrack, RduMpcTrack, ...
    mpcTrackPredictionSteps, mpcTrackControlHorizon, ...
    mpcTrackUseQuadprog, mpcTrackUseInputBounds, mpcTrackUseIncrementBounds, ...
    mpcTrackDeltaUMax, mpcTrackEpsX, mpcTrackEpsU);
exec_results.common.mpcTracker = mpcTracker;

%% Cost weights for planning-vs-execution comparisons
if isfield(planner_data,'opt') && isfield(planner_data.opt,'Qx') && isfield(planner_data.opt,'Ru')
    Qgoal = planner_data.opt.Qx;
    Rgoal = planner_data.opt.Ru;
else
    Qgoal = blkdiag(250*eye(3), 20*eye(3));
    Rgoal = 0.35*eye(3);
end
exec_results.common.Qgoal_optimality = Qgoal;
exec_results.common.Rgoal_optimality = Rgoal;
exec_results.common.optimality_metric_description = ...
    'For each selected reference: goal-regulation cost using x_goal=[q-qg;qd-qdg] and delta_tau command. Loss = J_exec_goal - J_plan_goal.';

%% Execute all trajectory-controller combinations
fprintf('\n====================================================================\n');
fprintf(' Nonlinear tracking execution for MPC/PID/LQR planner trajectories\n');
fprintf('====================================================================\n');
fprintf('Planner file                  : %s\n', plannerFile);
fprintf('Gravity compensation enabled  : %d\n', gravity_compensation_enabled);
fprintf('Delta-torque saturation       : %d\n', use_delta_tau_saturation);
fprintf('Absolute torque saturation    : %d\n', use_tau_abs_saturation);
fprintf('dt = %.4f s, Tf = %.4f s, N = %d\n', dt, Tf, N);
fprintf('Reduced endpoint body         : %s\n', reducedBody);
fprintf('====================================================================\n\n');

for itraj = 1:numel(trajectoryNames)
    trName = trajectoryNames{itraj};
    trLabel = trajectoryLabels{itraj};

    q_ref   = planner_data.(trName).q_ref;
    qd_ref  = planner_data.(trName).qd_ref;
    qdd_ref = planner_data.(trName).qdd_ref;

    % Nominal/planned delta torque for this selected reference.
    % Prefer saved planner input if available. Otherwise reconstruct from
    % inverse dynamics using the saved q, qd, qdd reference.
    u_nom = zeros(3,N+1);
    if isfield(planner_data.(trName),'u_ref') && isequal(size(planner_data.(trName).u_ref), [3 N])
        u_nom(:,1:N) = planner_data.(trName).u_ref;
        u_nom(:,N+1) = planner_data.(trName).u_ref(:,N);
    else
        for k = 1:N+1
            u_nom(:,k) = reducedInverseDynamicsFeedforward(robot, idx, q_fix_456, ...
                q_ref(:,k), qd_ref(:,k), qdd_ref(:,k));
        end
    end

    % Design midpoint-linearized LTI-LQR tracker about THIS selected trajectory.
    lqrDesign = designMidpointLTILQR( ...
        robot, idx, q_fix_456, q_ref, qd_ref, u_nom, t, ...
        Qlqr, Rlqr, gravity_compensation_enabled, ...
        linearization_eps_x, linearization_eps_u, lqr_uses_nominal_feedforward);

    % Save trajectory-level reference information for animation.
    trStruct = struct();
    trStruct.name = trName;
    trStruct.label = trLabel;
    trStruct.q_ref = q_ref;
    trStruct.qd_ref = qd_ref;
    trStruct.qdd_ref = qdd_ref;
    trStruct.u_nom = u_nom;
    trStruct.lqrDesign = lqrDesign;

    if isfield(planner_data.(trName),'endpoint_path')
        trStruct.planner_endpoint_path = planner_data.(trName).endpoint_path;
        trStruct.ref_endpoint = planner_data.(trName).endpoint_path;
    else
        trStruct.ref_endpoint = reducedEndpointSeries(robot, reducedBody, q_ref, q_fix_456);
        trStruct.planner_endpoint_path = trStruct.ref_endpoint;
    end

    fprintf('Reference trajectory: %s (%s)\n', trName, trLabel);
    fprintf('--------------------------------------------------------------------\n');

    for i = 1:numel(controllerNames)
        cname = controllerNames{i};

        simOut = simulateOptimalReferenceExecution( ...
            cname, robot, idx, q_fix_456, ...
            q0, qd0, q_ref, qd_ref, qdd_ref, u_nom, ...
            t, dt, Kp, Ki, Kd, ei_limit, lqrDesign, ...
            gravity_compensation_enabled, ...
            use_delta_tau_saturation, delta_tau_limit, ...
            use_tau_abs_saturation, tau_limit, ...
            Qexec, Rexec, compute_settling_time_flag, settling_threshold, settling_hold_samples, mpcTracker);

        simOut.controller_name = cname;
        simOut.controller_label = controllerLabels{i};
        simOut.trajectory_name = trName;
        simOut.trajectory_label = trLabel;
        simOut.endpoint = reducedEndpointSeries(robot, reducedBody, simOut.q, q_fix_456);
        simOut.ref_endpoint = trStruct.ref_endpoint;

        O = computePlanningVsExecutionOptimality( ...
            q_ref, qd_ref, u_nom, ...
            simOut.q, simOut.qd, simOut.delta_tau_cmd, ...
            qg, qdg, t, Qgoal, Rgoal);
        simOut.optimality = O;

        trStruct.trackers.(cname) = simOut;

        M = simOut.metrics;
        fprintf('%-10s | RMS e = %.6e | Final e = %.6e | J_track = %.6e | J_goal = %.6e | Loss = %.3f%% | Eabs = %.6e\n', ...
            cname, M.rms_position_error, M.final_position_error, M.final_executed_cost, ...
            O.final_executed_goal_cost, 100*O.relative_optimality_loss, M.final_cumulative_abs_energy);
    end

    exec_results.trajectories.(trName) = trStruct;
    fprintf('\n');
end

%% Backward compatibility aliases
% Older animation/debug scripts can still use exec_results.opt / pid / lqr.
exec_results.opt = exec_results.trajectories.opt.trackers.pid_ff;
exec_results.pid = exec_results.trajectories.pid.trackers.pid_ff;
exec_results.lqr = exec_results.trajectories.lqr.trackers.pid_ff;

% Older repaired single-reference code expected exec_results.trackers for opt.
exec_results.trackers = exec_results.trajectories.opt.trackers;

%% Save results
save(outputFile, 'exec_results');

fprintf('====================================================================\n');
fprintf('Saved all trajectory-tracker execution results to %s\n', outputFile);
fprintf('Use exec_results.trajectories.<opt|pid|lqr>.trackers.<pid|pid_ff|lqr|ltv_mpc>\n');
fprintf('====================================================================\n\n');

%% =========================================================================
%% PUBLICATION-QUALITY TRACKING RESULTS, TABLES, AND PLOTS
%% =========================================================================
% This section only formats, compares, plots, and exports tracking results.
% It does not change the planner, controller objectives, tracking laws,
% nonlinear execution, or saved logging structure used by the animation script.
if enableSummaryPlots
    trackingResultsRoot = fullfile(pwd, 'publication_tracking_results_ur5');
    generateTrackingPublicationResults(exec_results, trackingResultsRoot, plotTrajectoryName);
end

%% Local helper functions

function generateTrackingPublicationResults(exec_results, resultsRoot, focusTrajectoryName)
% Publication-quality tracking result export for Chapter 8.
% Outputs:
%   resultsRoot/figures/*.png, *.pdf, *.fig
%   resultsRoot/tables/*.csv, *.xlsx, *.tex
%   resultsRoot/data/tracking_publication_summary.mat

assert(isfield(exec_results,'trajectories'), 'exec_results.trajectories is missing.');
ensureDirTracking(resultsRoot);
figDir   = fullfile(resultsRoot, 'figures');
tableDir = fullfile(resultsRoot, 'tables');
dataDir  = fullfile(resultsRoot, 'data');
ensureDirTracking(figDir); ensureDirTracking(tableDir); ensureDirTracking(dataDir);

% Publication defaults matched to the planner plotting style.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultFigureColor', 'w');

pub = struct();
pub.lineWidthMain = 2.15;
pub.lineWidthAlt  = 1.75;
pub.fontSizeAxes  = 13;
pub.fontSizeTitle = 15;
pub.markerSize    = 48;
pub.dpi           = 600;
pub.formats       = {'png','pdf','fig'};
pub.styles        = {'-','--','-.',':'};
pub.markers       = {'o','s','^','d'};
pub.jointLabels   = {'$q_1$','$q_2$','$q_3$'};
pub.jointErrLabels = {'$e_{q,1}$','$e_{q,2}$','$e_{q,3}$'};
pub.torqueLabels  = {'$\Delta\tau_1$','$\Delta\tau_2$','$\Delta\tau_3$'};

trajectoryNames = getFieldOrDefault(exec_results.common, 'trajectory_names', fieldnames(exec_results.trajectories));
controllerNames = getFieldOrDefault(exec_results.common, 'controller_names', {'pid','pid_ff','lqr','ltv_mpc'});
trajectoryLabels = getFieldOrDefault(exec_results.common, 'trajectory_labels', trajectoryNames);
controllerLabels = getFieldOrDefault(exec_results.common, 'controller_labels', controllerNames);
trajectoryLabelsLatex = cellfun(@escapeLatexTracking, trajectoryLabels, 'UniformOutput', false);
controllerLabelsLatex = prettyControllerLabelsLatex(controllerNames, controllerLabels);

t = exec_results.common.t(:);
if nargin < 3 || isempty(focusTrajectoryName) || ~isfield(exec_results.trajectories, focusTrajectoryName)
    focusTrajectoryName = trajectoryNames{1};
end

% Build quantitative tables once from saved results.
summaryTable = buildTrackingSummaryTable(exec_results, trajectoryNames, controllerNames, trajectoryLabels, controllerLabels);
exportTrackingTable(summaryTable, tableDir, 'tracking_comparative_summary_all');

costTable     = summaryTable(:, {'Trajectory','Controller','FinalTrackingCost','FinalGoalCost','PlannedGoalCost','AbsoluteOptimalityLoss','RelativeOptimalityLossPercent'});
rmsTable      = summaryTable(:, {'Trajectory','Controller','RMSJointTrackingError','RMSVelocityTrackingError','FinalJointTrackingError','PeakJointTrackingError'});
effortTable   = summaryTable(:, {'Trajectory','Controller','RMSDeltaTau','IntegratedAbsDeltaTau','RMSFeedbackCorrection','IntegratedFeedbackEnergy','PeakAbsTorque','FinalAbsMechanicalEnergy'});
settleTable   = summaryTable(:, {'Trajectory','Controller','SettlingTime','FinalJointTrackingError','PeakJointTrackingError'});
optimalityTable = summaryTable(:, {'Trajectory','Controller','PlannedGoalCost','FinalGoalCost','AbsoluteOptimalityLoss','RelativeOptimalityLossPercent','ControlDeviationEnergy'});
exportTrackingTable(costTable, tableDir, 'cost_comparison');
exportTrackingTable(rmsTable, tableDir, 'rms_tracking_error_comparison');
exportTrackingTable(effortTable, tableDir, 'control_effort_comparison');
exportTrackingTable(settleTable, tableDir, 'settling_time_comparison');
exportTrackingTable(optimalityTable, tableDir, 'optimality_loss_comparison');

% Matrices used for heatmaps/bar charts.
metrics = trackingMetricMatrices(summaryTable, trajectoryNames, controllerNames);
save(fullfile(dataDir, 'tracking_publication_summary.mat'), 'summaryTable', 'metrics', 'trajectoryNames', 'controllerNames', 'trajectoryLabels', 'controllerLabels');

% Mandatory same-reference comparisons.  These are generated for every
% planner reference so the thesis can compare PID, PID+FF, LQR, and LTV-MPC
% trackers on the exact same reference trajectory, and also across planners.
for itrajPlot = 1:numel(trajectoryNames)
    trPlot = trajectoryNames{itrajPlot};
    plotSameReferenceJointTracking(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotSameReferenceTrackingErrors(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotSameReferenceControlEffort(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotSameReferenceEndpointTracking(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotSameReferenceScalarComparisons(metrics, trajectoryNames, controllerNames, trajectoryLabelsLatex, controllerLabelsLatex, trPlot, pub, figDir);
    plotCostOptimalityTimeHistories(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotCorrectionAndEnergy(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
    plotCostBreakdownBars(exec_results, trPlot, controllerNames, controllerLabelsLatex, pub, figDir);
end

% Multiple trajectory/controller comparisons.
plotMetricHeatmaps(metrics, trajectoryLabelsLatex, controllerLabelsLatex, pub, figDir);
plotAllTrajectoryEndpointComparisons(exec_results, trajectoryNames, trajectoryLabelsLatex, controllerNames, controllerLabelsLatex, pub, figDir);
plotTrackerTrajectoryTradeoff(summaryTable, trajectoryNames, controllerNames, trajectoryLabelsLatex, controllerLabelsLatex, pub, figDir);

fprintf('\nPublication tracking results saved to:\n  %s\n', resultsRoot);
fprintf('Figures: %s\nTables : %s\nData   : %s\n\n', figDir, tableDir, dataDir);
end

function summaryTable = buildTrackingSummaryTable(exec_results, trajectoryNames, controllerNames, trajectoryLabels, controllerLabels)
rows = {};
for it = 1:numel(trajectoryNames)
    tr = trajectoryNames{it};
    if ~isfield(exec_results.trajectories, tr), continue; end
    for ic = 1:numel(controllerNames)
        cn = controllerNames{ic};
        if ~isfield(exec_results.trajectories.(tr).trackers, cn), continue; end
        s = exec_results.trajectories.(tr).trackers.(cn);
        M = s.metrics;
        O = s.optimality;
        ctrlDevEnergy = NaN;
        if isfield(O,'control_deviation_energy') && ~isempty(O.control_deviation_energy)
            ctrlDevEnergy = O.control_deviation_energy(end);
        end
        rows(end+1,:) = { ... %#ok<AGROW>
            tr, trajectoryLabels{it}, cn, controllerLabels{ic}, ...
            safeMetric(M,'rms_position_error'), safeMetric(M,'rms_velocity_error'), ...
            safeMetric(M,'final_position_error'), safeMetric(M,'peak_position_error'), ...
            safeMetric(M,'rms_delta_tau'), safeMetric(M,'integrated_abs_delta_tau'), ...
            safeMetric(M,'rms_feedback_component'), safeMetric(M,'integrated_feedback_energy'), ...
            safeMetric(M,'peak_abs_torque'), safeMetric(M,'final_executed_cost'), ...
            safeMetric(M,'peak_abs_mechanical_power'), safeMetric(M,'final_cumulative_abs_energy'), ...
            safeMetric(M,'settling_time'), ...
            O.final_planned_goal_cost, O.final_executed_goal_cost, ...
            O.absolute_optimality_loss, 100*O.relative_optimality_loss, ctrlDevEnergy};
    end
end
summaryTable = cell2table(rows, 'VariableNames', { ...
    'Trajectory','TrajectoryLabel','Controller','ControllerLabel', ...
    'RMSJointTrackingError','RMSVelocityTrackingError', ...
    'FinalJointTrackingError','PeakJointTrackingError', ...
    'RMSDeltaTau','IntegratedAbsDeltaTau', ...
    'RMSFeedbackCorrection','IntegratedFeedbackEnergy', ...
    'PeakAbsTorque','FinalTrackingCost', ...
    'PeakAbsMechanicalPower','FinalAbsMechanicalEnergy', ...
    'SettlingTime','PlannedGoalCost','FinalGoalCost', ...
    'AbsoluteOptimalityLoss','RelativeOptimalityLossPercent','ControlDeviationEnergy'});
end

function metrics = trackingMetricMatrices(T, trajectoryNames, controllerNames)
metricNames = {'RMSJointTrackingError','FinalTrackingCost','FinalGoalCost', ...
    'RMSDeltaTau','IntegratedAbsDeltaTau','SettlingTime', ...
    'AbsoluteOptimalityLoss','RelativeOptimalityLossPercent', ...
    'FinalAbsMechanicalEnergy','ControlDeviationEnergy'};
for im = 1:numel(metricNames)
    A = NaN(numel(trajectoryNames), numel(controllerNames));
    for it = 1:numel(trajectoryNames)
        for ic = 1:numel(controllerNames)
            mask = strcmp(T.Trajectory, trajectoryNames{it}) & strcmp(T.Controller, controllerNames{ic});
            if any(mask), A(it,ic) = T{find(mask,1), metricNames{im}}; end
        end
    end
    metrics.(metricNames{im}) = A;
end
end

function plotSameReferenceJointTracking(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
t = exec_results.common.t(:);
qref = exec_results.trajectories.(tr).q_ref;
trLabel = escapeLatexTracking(exec_results.trajectories.(tr).label);
fig = figure('Name',['Joint tracking - ', tr], 'Color','w', 'Position',[80 80 980 820]);
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for j = 1:3
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, qref(j,:).', 'k--', 'LineWidth', pub.lineWidthMain);
    for ic = 1:numel(controllerNames)
        cn = controllerNames{ic};
        q = exec_results.trajectories.(tr).trackers.(cn).q;
        plot(ax, t, q(j,:).', pub.styles{ic}, 'LineWidth', pub.lineWidthAlt);
    end
    ylabel(ax, sprintf('%s [rad]', pub.jointLabels{j}), 'FontSize', pub.fontSizeAxes);
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
    if j == 1
        title(ax, ['Joint-space tracking for ', trLabel], 'FontSize', pub.fontSizeTitle);
        legend(ax, [{'Reference'}, controllerLabelsLatex], 'Location','best', 'NumColumns',3);
    end
    if j == 3, xlabel(ax, '$t$ [s]', 'FontSize', pub.fontSizeAxes); end
end
title(tl, '$q(t)$ versus $q_{\mathrm{ref}}(t)$ for the same reference trajectory', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_joint_tracking_', tr], pub);
end

function plotSameReferenceTrackingErrors(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
t = exec_results.common.t(:);
trLabel = escapeLatexTracking(exec_results.trajectories.(tr).label);
fig = figure('Name',['Joint tracking errors - ', tr], 'Color','w', 'Position',[100 100 980 820]);
tl = tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
for j = 1:3
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for ic = 1:numel(controllerNames)
        cn = controllerNames{ic};
        e = exec_results.trajectories.(tr).trackers.(cn).err_q;
        plot(ax, t, e(j,:).', pub.styles{ic}, 'LineWidth', pub.lineWidthAlt);
    end
    ylabel(ax, sprintf('%s [rad]', pub.jointErrLabels{j}), 'FontSize', pub.fontSizeAxes);
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
    if j == 1
        title(ax, ['Joint tracking error for ', trLabel], 'FontSize', pub.fontSizeTitle);
        legend(ax, controllerLabelsLatex, 'Location','best', 'NumColumns',4);
    end
end
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    plot(ax, t, exec_results.trajectories.(tr).trackers.(cn).err_norm(:), pub.styles{ic}, 'LineWidth', pub.lineWidthMain);
end
ylabel(ax, '$\Vert e_q\Vert_2$ [rad]', 'FontSize', pub.fontSizeAxes);
xlabel(ax, '$t$ [s]', 'FontSize', pub.fontSizeAxes);
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
title(tl, 'Joint-space tracking error comparison', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_tracking_errors_', tr], pub);
end

function plotSameReferenceControlEffort(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
t = exec_results.common.t(:);
uNom = exec_results.trajectories.(tr).u_nom;
trLabel = escapeLatexTracking(exec_results.trajectories.(tr).label);
fig = figure('Name',['Control effort - ', tr], 'Color','w', 'Position',[120 120 980 820]);
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for j = 1:3
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    plot(ax, t, uNom(j,:).', 'k--', 'LineWidth', pub.lineWidthMain);
    for ic = 1:numel(controllerNames)
        cn = controllerNames{ic};
        u = exec_results.trajectories.(tr).trackers.(cn).delta_tau_cmd;
        plot(ax, t, u(j,:).', pub.styles{ic}, 'LineWidth', pub.lineWidthAlt);
    end
    ylabel(ax, sprintf('%s [N m]', pub.torqueLabels{j}), 'FontSize', pub.fontSizeAxes);
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
    if j == 1
        title(ax, ['Commanded control effort for ', trLabel], 'FontSize', pub.fontSizeTitle);
        legend(ax, [{'Nominal'}, controllerLabelsLatex], 'Location','best', 'NumColumns',3);
    end
    if j == 3, xlabel(ax, '$t$ [s]', 'FontSize', pub.fontSizeAxes); end
end
title(tl, '$\Delta\tau(t)$ comparison for the same reference trajectory', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_control_effort_', tr], pub);
end

function plotSameReferenceEndpointTracking(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
pref = exec_results.trajectories.(tr).ref_endpoint;
trLabel = escapeLatexTracking(exec_results.trajectories.(tr).label);
fig = figure('Name',['End-effector tracking - ', tr], 'Color','w', 'Position',[140 140 900 760]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
plot3(ax, pref(1,:), pref(2,:), pref(3,:), 'k--', 'LineWidth', pub.lineWidthMain);
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    p = exec_results.trajectories.(tr).trackers.(cn).endpoint;
    plot3(ax, p(1,:), p(2,:), p(3,:), pub.styles{ic}, 'LineWidth', pub.lineWidthAlt);
end
scatter3(ax, pref(1,1), pref(2,1), pref(3,1), pub.markerSize, 'filled', 'Marker','o');
scatter3(ax, pref(1,end), pref(2,end), pref(3,end), pub.markerSize, 'filled', 'Marker','s');
xlabel(ax, '$x$ [m]'); ylabel(ax, '$y$ [m]'); zlabel(ax, '$z$ [m]');
title(ax, ['End-effector executed path versus reference: ', trLabel], 'FontSize', pub.fontSizeTitle);
legend(ax, [{'Reference'}, controllerLabelsLatex, {'Start','Goal'}], 'Location','bestoutside');
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0); view(ax, 3);
savePubFigureTracking(fig, figDir, ['same_reference_endpoint_tracking_', tr], pub);
end

function plotSameReferenceScalarComparisons(metrics, trajectoryNames, controllerNames, trajectoryLabelsLatex, controllerLabelsLatex, tr, pub, figDir)
it = find(strcmp(trajectoryNames, tr), 1);
if isempty(it), it = 1; tr = trajectoryNames{1}; end
labels = categorical(controllerLabelsLatex); labels = reordercats(labels, controllerLabelsLatex);
fig = figure('Name',['Scalar comparison - ', tr], 'Color','w', 'Position',[100 100 1100 780]);
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
metricList = {'RMSJointTrackingError','FinalTrackingCost','RelativeOptimalityLossPercent','IntegratedAbsDeltaTau'};
ylabels = {'RMS $\Vert e_q\Vert_2$ [rad]', '$J_{\mathrm{track}}$', 'Optimality loss [\%]', '$\int \sum_i |\Delta\tau_i|\,dt$'};
titles = {'RMS tracking error','Cumulative tracking cost','Optimality loss','Integrated control effort'};
for k = 1:4
    ax = nexttile; grid(ax,'on'); box(ax,'on');
    bar(ax, labels, metrics.(metricList{k})(it,:));
    ylabel(ax, ylabels{k}, 'FontSize', pub.fontSizeAxes);
    title(ax, titles{k}, 'FontSize', pub.fontSizeTitle);
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
end
title(tl, ['Controller comparison on ', trajectoryLabelsLatex{it}], 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_scalar_comparison_', tr], pub);
end

function plotCostOptimalityTimeHistories(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
t = exec_results.common.t(:);
tCost = t(1:end-1);
fig = figure('Name',['Cumulative cost and optimality - ', tr], 'Color','w', 'Position',[130 130 1080 760]);
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    plot(ax, t, exec_results.trajectories.(tr).trackers.(cn).Jcum_exec(:), pub.styles{ic}, 'LineWidth', pub.lineWidthMain);
end
ylabel(ax, '$J_{\mathrm{track}}(0,t)$');
title(ax, 'Cumulative tracking cost');
legend(ax, controllerLabelsLatex, 'Location','best', 'NumColumns',4);
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    O = exec_results.trajectories.(tr).trackers.(cn).optimality;
    plot(ax, tCost, O.Jcum_exec_goal(:) - O.Jcum_plan(:), pub.styles{ic}, 'LineWidth', pub.lineWidthMain);
end
xlabel(ax, '$t$ [s]'); ylabel(ax, '$J_{\mathrm{exec}}-J_{\mathrm{plan}}$');
title(ax, 'Accumulated optimality gap relative to the planner objective');
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
title(tl, 'Cost and optimality-loss histories', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_cost_optimality_histories_', tr], pub);
end

function plotCorrectionAndEnergy(exec_results, tr, controllerNames, controllerLabelsLatex, pub, figDir)
t = exec_results.common.t(:);
fig = figure('Name',['Correction and energy - ', tr], 'Color','w', 'Position',[160 160 1080 760]);
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    fb = exec_results.trajectories.(tr).trackers.(cn).feedback_component;
    plot(ax, t, vecnorm(fb,2,1), pub.styles{ic}, 'LineWidth', pub.lineWidthMain);
end
ylabel(ax, '$\Vert\Delta\tau_{\mathrm{fb}}\Vert_2$ [N m]');
title(ax, 'Feedback correction magnitude');
legend(ax, controllerLabelsLatex, 'Location','best', 'NumColumns',4);
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    plot(ax, t, exec_results.trajectories.(tr).trackers.(cn).EmechAbs(:), pub.styles{ic}, 'LineWidth', pub.lineWidthMain);
end
xlabel(ax, '$t$ [s]'); ylabel(ax, '$E_{\mathrm{abs}}$ [J]');
title(ax, 'Cumulative absolute mechanical energy');
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
title(tl, 'Origin and energetic consequence of tracking corrections', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, ['same_reference_correction_energy_', tr], pub);
end

function plotMetricHeatmaps(metrics, trajectoryLabelsLatex, controllerLabelsLatex, pub, figDir)
items = { ...
    'RMSJointTrackingError', 'RMS $\Vert e_q\Vert_2$ [rad]', 'heatmap_rms_tracking_error'; ...
    'FinalTrackingCost', '$J_{\mathrm{track}}$', 'heatmap_cumulative_tracking_cost'; ...
    'RMSDeltaTau', 'RMS $\Delta\tau$ [N m]', 'heatmap_rms_control_effort'; ...
    'SettlingTime', '$t_s$ [s]', 'heatmap_settling_time'; ...
    'RelativeOptimalityLossPercent', 'Optimality loss [\%]', 'heatmap_relative_optimality_loss'; ...
    'FinalAbsMechanicalEnergy', '$E_{\mathrm{abs}}$ [J]', 'heatmap_abs_mechanical_energy'};
for k = 1:size(items,1)
    fig = figure('Name',items{k,3}, 'Color','w', 'Position',[180 180 900 620]);
    ax = axes(fig);
    imagesc(ax, metrics.(items{k,1}));
    colorbar(ax);
    xticks(ax, 1:numel(controllerLabelsLatex)); xticklabels(ax, controllerLabelsLatex);
    yticks(ax, 1:numel(trajectoryLabelsLatex)); yticklabels(ax, trajectoryLabelsLatex);
    xlabel(ax, 'Tracking controller'); ylabel(ax, 'Reference trajectory');
    title(ax, items{k,2}, 'FontSize', pub.fontSizeTitle);
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0, 'TickLabelInterpreter','latex');
    annotateMatrixTracking(ax, metrics.(items{k,1}));
    savePubFigureTracking(fig, figDir, items{k,3}, pub);
end
end

function plotAllTrajectoryEndpointComparisons(exec_results, trajectoryNames, trajectoryLabelsLatex, controllerNames, controllerLabelsLatex, pub, figDir)
fig = figure('Name','All trajectory endpoint comparisons', 'Color','w', 'Position',[80 80 1180 900]);
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for it = 1:numel(trajectoryNames)
    tr = trajectoryNames{it};
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
    pref = exec_results.trajectories.(tr).ref_endpoint;
    plot3(ax, pref(1,:), pref(2,:), pref(3,:), 'k--', 'LineWidth', pub.lineWidthMain);
    for ic = 1:numel(controllerNames)
        cn = controllerNames{ic};
        p = exec_results.trajectories.(tr).trackers.(cn).endpoint;
        plot3(ax, p(1,:), p(2,:), p(3,:), pub.styles{ic}, 'LineWidth', pub.lineWidthAlt);
    end
    title(ax, trajectoryLabelsLatex{it}, 'FontSize', pub.fontSizeTitle);
    xlabel(ax,'$x$ [m]'); ylabel(ax,'$y$ [m]'); zlabel(ax,'$z$ [m]');
    set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0); view(ax,3);
    if it == 1
        legend(ax, [{'Reference'}, controllerLabelsLatex], 'Location','best');
    end
end
title(tl, 'Executed end-effector paths for all references and controllers', 'FontSize', pub.fontSizeTitle);
savePubFigureTracking(fig, figDir, 'all_trajectories_endpoint_tracking_grid', pub);
end

function plotTrackerTrajectoryTradeoff(T, trajectoryNames, controllerNames, trajectoryLabelsLatex, controllerLabelsLatex, pub, figDir)
fig = figure('Name','Tracking error versus optimality loss tradeoff', 'Color','w', 'Position',[100 100 950 720]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for ic = 1:numel(controllerNames)
    xs = []; ys = [];
    for it = 1:numel(trajectoryNames)
        mask = strcmp(T.Trajectory, trajectoryNames{it}) & strcmp(T.Controller, controllerNames{ic});
        if any(mask)
            xs(end+1) = T.RMSJointTrackingError(find(mask,1)); %#ok<AGROW>
            ys(end+1) = T.RelativeOptimalityLossPercent(find(mask,1)); %#ok<AGROW>
        end
    end
    plot(ax, xs, ys, pub.styles{ic}, 'Marker', pub.markers{ic}, 'LineWidth', pub.lineWidthMain, 'MarkerSize', 8);
end
xlabel(ax, 'RMS joint tracking error [rad]');
ylabel(ax, 'Relative optimality loss [\%]');
title(ax, 'Tracking accuracy versus preservation of planner optimality', 'FontSize', pub.fontSizeTitle);
legend(ax, controllerLabelsLatex, 'Location','best');
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
for it = 1:numel(trajectoryNames)
    mask = strcmp(T.Trajectory, trajectoryNames{it}) & strcmp(T.Controller, controllerNames{end});
    if any(mask)
        text(ax, T.RMSJointTrackingError(find(mask,1)), T.RelativeOptimalityLossPercent(find(mask,1)), ['  ', trajectoryLabelsLatex{it}], 'FontSize', pub.fontSizeAxes-1);
    end
end
savePubFigureTracking(fig, figDir, 'tracking_error_vs_optimality_loss_tradeoff', pub);
end

function plotCostBreakdownBars(exec_results, focusTrajectoryName, controllerNames, controllerLabelsLatex, pub, figDir)
if ~isfield(exec_results,'common') || ~isfield(exec_results.common,'Qexec') || ~isfield(exec_results.common,'Rexec')
    return;
end
Qexec = exec_results.common.Qexec;
Rexec = exec_results.common.Rexec;
stateCost = NaN(1,numel(controllerNames));
controlCost = NaN(1,numel(controllerNames));
for ic = 1:numel(controllerNames)
    cn = controllerNames{ic};
    if ~isfield(exec_results.trajectories.(focusTrajectoryName).trackers, cn), continue; end
    s = exec_results.trajectories.(focusTrajectoryName).trackers.(cn);
    Nsamples = size(s.err_q,2);
    Js = zeros(1,Nsamples);
    Ju = zeros(1,Nsamples);
    for k = 1:Nsamples
        ex = [s.err_q(:,k); s.err_qd(:,k)];
        uk = s.delta_tau_cmd(:,k);
        Js(k) = ex.'*Qexec*ex;
        Ju(k) = uk.'*Rexec*uk;
    end
    stateCost(ic) = sum(Js);
    controlCost(ic) = sum(Ju);
end
labels = categorical(controllerLabelsLatex); labels = reordercats(labels, controllerLabelsLatex);
fig = figure('Name',['Cost breakdown - ', focusTrajectoryName], 'Color','w', 'Position',[100 100 900 650]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
bar(ax, labels, [stateCost(:) controlCost(:)], 'stacked');
ylabel(ax, '$J_{\mathrm{track}}$');
title(ax, 'Exact tracking-cost decomposition: state error versus control effort', 'FontSize', pub.fontSizeTitle);
legend(ax, {'State/error cost','Control-effort cost'}, 'Location','best');
set(ax,'FontSize',pub.fontSizeAxes,'LineWidth',1.0);
savePubFigureTracking(fig, figDir, ['cost_breakdown_exact_', focusTrajectoryName], pub);
end

function exportTrackingTable(T, tableDir, baseName)
writetable(T, fullfile(tableDir, [baseName '.csv']));
try
    writetable(T, fullfile(tableDir, [baseName '.xlsx']));
catch ME
    warning('Could not write Excel table %s: %s', baseName, ME.message);
end
writeLatexTableTracking(T, fullfile(tableDir, [baseName '.tex']));
end

function writeLatexTableTracking(T, filePath)
fid = fopen(filePath, 'w');
if fid < 0, warning('Could not open %s for writing.', filePath); return; end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
vars = T.Properties.VariableNames;
fprintf(fid, '%% Auto-generated by tracking publication results script.\n');
fprintf(fid, '\\begin{table}[!t]\n\\centering\n\\scriptsize\n');
fprintf(fid, '\\begin{tabular}{%s}\n\\hline\n', repmat('l',1,numel(vars)));
fprintf(fid, '%s \\\\ \\hline\n', strjoin(cellfun(@escapeLatexTracking, vars, 'UniformOutput', false), ' & '));
for i = 1:height(T)
    vals = cell(1,numel(vars));
    for j = 1:numel(vars)
        v = T{i,j};
        if iscell(v), v = v{1}; end
        if isnumeric(v)
            if isnan(v), vals{j} = '--'; else, vals{j} = sprintf('%.4g', v); end
        else
            vals{j} = escapeLatexTracking(char(string(v)));
        end
    end
    fprintf(fid, '%s \\\\ \n', strjoin(vals, ' & '));
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fprintf(fid, '\\caption{Auto-generated comparative tracking results.}\n');
fprintf(fid, '\\end{table}\n');
end

function savePubFigureTracking(fig, figDir, baseName, pub)
set(fig, 'Color','w');
for i = 1:numel(pub.formats)
    fmt = pub.formats{i};
    out = fullfile(figDir, [baseName '.' fmt]);
    try
        switch lower(fmt)
            case 'png'
                exportgraphics(fig, out, 'Resolution', pub.dpi);
            case 'pdf'
                exportgraphics(fig, out, 'ContentType','vector');
            case 'fig'
                savefig(fig, out);
        end
    catch
        if strcmpi(fmt,'png')
            print(fig, out, '-dpng', sprintf('-r%d', pub.dpi));
        elseif strcmpi(fmt,'pdf')
            print(fig, out, '-dpdf', '-painters');
        elseif strcmpi(fmt,'fig')
            saveas(fig, out);
        end
    end
end
end

function annotateMatrixTracking(ax, A)
[nr,nc] = size(A);
for r = 1:nr
    for c = 1:nc
        val = A(r,c);
        if isnan(val), label = '--'; else, label = sprintf('%.3g', val); end
        text(ax, c, r, label, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 10, 'Interpreter','latex');
    end
end
end

function labelsLatex = prettyControllerLabelsLatex(controllerNames, controllerLabels)
labelsLatex = cell(size(controllerNames));
for i = 1:numel(controllerNames)
    switch lower(controllerNames{i})
        case 'pid'
            labelsLatex{i} = 'PID';
        case 'pid_ff'
            labelsLatex{i} = 'PID+FF';
        case 'lqr'
            labelsLatex{i} = 'LQR';
        case 'ltv_mpc'
            labelsLatex{i} = 'LTV-MPC tracker';
        otherwise
            labelsLatex{i} = escapeLatexTracking(controllerLabels{i});
    end
end
end

function v = safeMetric(M, name)
if isfield(M,name), v = M.(name); else, v = NaN; end
end

function out = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    out = s.(fieldName);
else
    out = defaultValue;
end
end

function s = escapeLatexTracking(s)
s = char(string(s));
s = strrep(s, '\', '\textbackslash{}');
s = strrep(s, '_', '\_');
s = strrep(s, '%', '\%');
s = strrep(s, '&', '\&');
s = strrep(s, '#', '\#');
s = strrep(s, '{', '\{');
s = strrep(s, '}', '\}');
end

function ensureDirTracking(d)
if exist(d,'dir') ~= 7
    mkdir(d);
end
end



function simOut = simulateOptimalReferenceExecution( ...
    controllerName, robot, idx, q_fix_456, ...
    q0, qd0, q_ref, qd_ref, qdd_ref, u_nom, ...
    t, dt, Kp, Ki, Kd, ei_limit, lqrDesign, ...
    gravity_compensation_enabled, ...
    use_delta_tau_saturation, delta_tau_limit, ...
    use_tau_abs_saturation, tau_limit, ...
    Qexec, Rexec, compute_settling_time_flag, settling_threshold, settling_hold_samples, mpcTracker)

numSamples = numel(t);
q               = zeros(3, numSamples);
qd              = zeros(3, numSamples);
qdd             = zeros(3, numSamples);
delta_tau_ff    = zeros(3, numSamples);
delta_tau_fb    = zeros(3, numSamples);
delta_tau_cmd   = zeros(3, numSamples);
tau_abs         = zeros(3, numSamples);
err_q           = zeros(3, numSamples);
err_qd          = zeros(3, numSamples);
err_norm        = zeros(1, numSamples);
feedback_component = zeros(3, numSamples);
Jstage_exec     = zeros(1, numSamples);
Jcum_exec       = zeros(1, numSamples);
mpc_exitflag    = NaN(1, numSamples);
mpc_horizon_used = NaN(1, numSamples);
mpc_du0         = NaN(3, numSamples);

x = [q0(:); qd0(:)];
ei = zeros(3,1);

for k = 1:numSamples
    q(:,k)  = x(1:3);
    qd(:,k) = x(4:6);

    eq_ref_minus_exec  = q_ref(:,k)  - q(:,k);
    eqd_ref_minus_exec = qd_ref(:,k) - qd(:,k);
    err_q(:,k)  = eq_ref_minus_exec;
    err_qd(:,k) = eqd_ref_minus_exec;
    err_norm(k) = norm(eq_ref_minus_exec);

    % Feedforward/nominal component
    switch lower(controllerName)
        case 'pid'
            delta_tau_ff_k = zeros(3,1);
        case 'pid_ff'
            delta_tau_ff_k = reducedInverseDynamicsFeedforward(robot, idx, q_fix_456, q_ref(:,k), qd_ref(:,k), qdd_ref(:,k));
        case 'lqr'
            if isfield(lqrDesign,'use_nominal_feedforward') && lqrDesign.use_nominal_feedforward
                delta_tau_ff_k = u_nom(:,k);
            else
                delta_tau_ff_k = zeros(3,1);
            end
        case 'ltv_mpc'
            % The LTV MPC tracker optimizes the full commanded delta torque
            % through increments, so keep the plotted feedforward part as the
            % nominal/reference command and the feedback part as the MPC correction.
            delta_tau_ff_k = u_nom(:,k);
        otherwise
            error('Unknown controllerName: %s', controllerName);
    end

    % Feedback/correction component
    switch lower(controllerName)
        case {'pid','pid_ff'}
            ei = ei + eq_ref_minus_exec * dt;
            ei = min(max(ei, -ei_limit), ei_limit);
            correction_k = Kp*eq_ref_minus_exec + Ki*ei + Kd*eqd_ref_minus_exec;
        case 'lqr'
            x_tilde = [q(:,k) - q_ref(:,k); qd(:,k) - qd_ref(:,k)];
            correction_k = -lqrDesign.K(:,:,k) * x_tilde;
        case 'ltv_mpc'
            % Repaired MPC tracking: optimize only the correction c around
            % the nominal planner torque u_nom, not the full torque from zero.
            % This is the important change that makes the MPC tracker comparable
            % to LQR/PID+FF while preserving the Sec. 10.3 incremental-input idea.
            c_prev_for_mpc = zeros(3,1);
            if k > 1
                c_prev_for_mpc = delta_tau_fb(:,k-1); % previous correction = previous cmd - previous u_nom
            end
            [delta_tau_cmd_from_mpc, mpcStepLog] = solveLtvMpcTrackingStep( ...
                robot, idx, q_fix_456, x, c_prev_for_mpc, ...
                q_ref, qd_ref, u_nom, k, dt, gravity_compensation_enabled, ...
                mpcTracker, use_delta_tau_saturation, delta_tau_limit);
            correction_k = delta_tau_cmd_from_mpc - delta_tau_ff_k;
            mpc_exitflag(k) = mpcStepLog.exitflag;
            mpc_horizon_used(k) = mpcStepLog.Np;
            if isfield(mpcStepLog,'dc0')
                mpc_du0(:,k) = mpcStepLog.dc0;
            elseif isfield(mpcStepLog,'du0')
                mpc_du0(:,k) = mpcStepLog.du0;
            end
        otherwise
            error('Unknown controllerName: %s', controllerName);
    end

    delta_tau_cmd_k = delta_tau_ff_k + correction_k;

    if use_delta_tau_saturation
        delta_tau_cmd_k = min(max(delta_tau_cmd_k, -delta_tau_limit), delta_tau_limit);
    end

    tau_abs_k = composeAbsoluteTorque(robot, idx, q_fix_456, q(:,k), delta_tau_cmd_k, gravity_compensation_enabled);

    if use_tau_abs_saturation
        tau_abs_k = min(max(tau_abs_k, -tau_limit), tau_limit);
        if gravity_compensation_enabled
            G3 = reducedGravity(robot, idx, q_fix_456, q(:,k));
            delta_tau_cmd_k = tau_abs_k - G3;
        else
            delta_tau_cmd_k = tau_abs_k;
        end
    end

    delta_tau_ff(:,k) = delta_tau_ff_k;
    delta_tau_cmd(:,k) = delta_tau_cmd_k;
    delta_tau_fb(:,k) = delta_tau_cmd_k - delta_tau_ff_k;
    feedback_component(:,k) = delta_tau_cmd_k - delta_tau_ff_k;
    tau_abs(:,k) = tau_abs_k;

    xdot_now = reducedStateDerivative(robot, idx, q_fix_456, x, delta_tau_cmd_k, gravity_compensation_enabled);
    qdd(:,k) = xdot_now(4:6);

    ex = [eq_ref_minus_exec; eqd_ref_minus_exec];
    Jstage_exec(k) = ex.'*Qexec*ex + delta_tau_cmd_k.'*Rexec*delta_tau_cmd_k;
    if k == 1
        Jcum_exec(k) = Jstage_exec(k);
    else
        Jcum_exec(k) = Jcum_exec(k-1) + Jstage_exec(k);
    end

    if k < numSamples
        x = rk4StepReducedPlant(robot, idx, q_fix_456, x, delta_tau_cmd_k, dt, gravity_compensation_enabled);
    end
end

PjointSigned = tau_abs .* qd;
PjointAbs    = abs(PjointSigned);
PmechSigned  = sum(PjointSigned,1);
PmechAbs     = sum(PjointAbs,1);
EmechSigned  = cumtrapz(t, PmechSigned);
EmechAbs     = cumtrapz(t, PmechAbs);
correction_energy_inst = sum(feedback_component.^2,1);
correction_energy_cum  = cumtrapz(t, correction_energy_inst);

metrics = struct();
metrics.rms_position_error             = sqrt(mean(sum(err_q.^2,1)));
metrics.final_position_error           = err_norm(end);
metrics.peak_position_error            = max(err_norm);
metrics.rms_velocity_error             = sqrt(mean(sum(err_qd.^2,1)));
metrics.rms_delta_tau                  = sqrt(mean(sum(delta_tau_cmd.^2,1)));
metrics.integrated_abs_delta_tau       = trapz(t, sum(abs(delta_tau_cmd),1));
metrics.rms_feedback_component         = sqrt(mean(sum(feedback_component.^2,1)));
metrics.integrated_feedback_energy     = correction_energy_cum(end);
metrics.peak_abs_torque                = max(abs(tau_abs), [], 'all');
metrics.final_executed_cost            = Jcum_exec(end);
metrics.peak_abs_mechanical_power      = max(PmechAbs);
metrics.final_cumulative_abs_energy    = EmechAbs(end);
if compute_settling_time_flag
    metrics.settling_time = computeSettlingTime(err_norm, t, settling_threshold, settling_hold_samples);
else
    metrics.settling_time = NaN;
end

simOut = struct();
simOut.q = q;
simOut.qd = qd;
simOut.qdd = qdd;
simOut.delta_tau_ff = delta_tau_ff;
simOut.delta_tau_fb = delta_tau_fb;
simOut.delta_tau_cmd = delta_tau_cmd;
simOut.feedback_component = feedback_component;
simOut.tau_abs = tau_abs;
simOut.err_q = err_q;
simOut.err_qd = err_qd;
simOut.err_norm = err_norm;
simOut.Jstage_exec = Jstage_exec;
simOut.Jcum_exec = Jcum_exec;
simOut.PjointSigned = PjointSigned;
simOut.PjointAbs = PjointAbs;
simOut.PmechSigned = PmechSigned;
simOut.PmechAbs = PmechAbs;
simOut.EmechSigned = EmechSigned;
simOut.EmechAbs = EmechAbs;
simOut.correction_energy_inst = correction_energy_inst;
simOut.correction_energy_cum = correction_energy_cum;
simOut.metrics = metrics;
simOut.q_ref = q_ref;
simOut.qd_ref = qd_ref;
simOut.qdd_ref = qdd_ref;
simOut.u_nom = u_nom;
simOut.mpc.exitflag = mpc_exitflag;
simOut.mpc.horizon_used = mpc_horizon_used;
simOut.mpc.dc0 = mpc_du0;
simOut.mpc.du0 = mpc_du0; % backward-compatible alias; here this is Delta correction, not full Delta torque
simOut.mpc.description = 'For ltv_mpc only: online QP diagnostics. dc0 is the first receding-horizon correction increment, with delta_tau_cmd = u_nom + correction.';
end


function mpcTracker = designLtvMpcTracker(Q, Qf, Rdu, Np, Nc, useQuadprog, useInputBounds, useIncrementBounds, deltaUMax, eps_x, eps_u)
% Design/settings container for the online LTV MPC tracking controller.
% The actual LTV model is rebuilt online because the nonlinear plant is
% linearized along the future reference segment at every execution sample.
mpcTracker = struct();
mpcTracker.type = 'online_receding_horizon_ltv_mpc_tracking';
mpcTracker.Q = Q;
mpcTracker.Qf = Qf;
mpcTracker.Rdu = Rdu;
mpcTracker.Np = Np;
mpcTracker.Nc = Nc;
mpcTracker.useQuadprog = useQuadprog;
mpcTracker.useInputBounds = useInputBounds;
mpcTracker.useIncrementBounds = useIncrementBounds;
mpcTracker.deltaUMax = deltaUMax(:);
mpcTracker.eps_x = eps_x;
mpcTracker.eps_u = eps_u;
mpcTracker.outputMatrix = [eye(6), zeros(6,3)];
mpcTracker.description = ['Book Sec. 10.3 style tracking MPC with incremental input: ', ...
    'u(k)=u(k-1)+Delta u(k). The augmented state is zeta=[z;u(k-1)], ', ...
    'z=[q-qref;qd-qdref]. A finite-horizon QP minimizes predicted tracking ', ...
    'error and Delta u effort, and only the first Delta u is applied.'];
end

function [u_cmd, info] = solveLtvMpcTrackingStep( ...
    robot, idx, q_fix_456, x_now, c_prev, ...
    q_ref, qd_ref, u_nom, k_now, dt, gravity_compensation_enabled, ...
    mpcTracker, use_delta_tau_saturation, delta_tau_limit)
% Repaired receding-horizon online LTV MPC tracker.
%
% This version is centered on the planner-generated nominal trajectory:
%   delta_tau_cmd(k) = u_nom(k) + c(k)
% where c(k) is a correction.  The optimized variable follows the tracking
% MPC idea from Sec. 10.3 of the note/book chapter, but applied to correction:
%   c(k) = c(k-1) + Delta c(k)
%
% Error dynamics are built along the future MPC reference:
%   e_{i+1} = A_i e_i + B_i c_i + g_i
%   g_i = x_nom_next - x_ref_next
% where x_nom_next is obtained by propagating the nonlinear model from
% x_ref_i using u_nom_i.  This defect term is important whenever the saved
% reference and saved nominal input are not exactly consistent after sampling,
% numerical differentiation, endpoint clipping, or planner post-processing.
%
% Augmented correction model:
%   zeta_i = [e_i; c_{i-1}]
%   zeta_{i+1} = [A_i B_i; 0 I] zeta_i + [B_i; I] Delta c_i + [g_i; 0]
% The QP minimizes predicted tracking error and correction increments.  Only
% the first correction increment is applied, then the problem is rebuilt at
% the next sample from the updated nonlinear state.

numSamples = size(q_ref,2);
Np = min(mpcTracker.Np, numSamples - k_now);
if Np < 1
    u_cmd = u_nom(:,k_now) + c_prev(:);
    if use_delta_tau_saturation
        u_cmd = min(max(u_cmd, -delta_tau_limit(:)), delta_tau_limit(:));
    end
    info = struct('exitflag',0,'Np',0,'message','terminal sample, reused previous correction');
    return;
end

nx = 6; nu = 3; nxa = nx + nu;
e0 = x_now(:) - [q_ref(:,k_now); qd_ref(:,k_now)];
zeta0 = [e0; c_prev(:)];

Aseq = cell(1,Np);
Bseq = cell(1,Np);
gseq = cell(1,Np);
for i = 1:Np
    kk = min(k_now + i - 1, numSamples);
    kkNext = min(k_now + i, numSamples);

    x_ref_i = [q_ref(:,kk); qd_ref(:,kk)];
    x_ref_next = [q_ref(:,kkNext); qd_ref(:,kkNext)];
    u_nom_i = u_nom(:,kk);

    [Ad_i, Bd_i, x_nom_next] = discreteLinearizeAroundNominalLocal( ...
        robot, idx, q_fix_456, x_ref_i, u_nom_i, dt, ...
        gravity_compensation_enabled, mpcTracker.eps_x, mpcTracker.eps_u);

    g_i = x_nom_next - x_ref_next;

    Aseq{i} = Ad_i;
    Bseq{i} = Bd_i;
    gseq{i} = g_i;
end

% Build augmented affine model in the correction increment Delta c.
AaugSeq = cell(1,Np);
BaugSeq = cell(1,Np);
caugSeq = cell(1,Np);
for i = 1:Np
    AaugSeq{i} = [Aseq{i}, Bseq{i}; zeros(nu,nx), eye(nu)];
    BaugSeq{i} = [Bseq{i}; eye(nu)];
    caugSeq{i} = [gseq{i}; zeros(nu,1)];
end

[Abar, Bbar, cbar] = buildStackedPredictionMatricesLTVAffineLocal(AaugSeq, BaugSeq, caugSeq, Np);
Ctilde = mpcTracker.outputMatrix; % selects e from zeta=[e;c]
Qbar = zeros(nxa*Np);
for i = 1:Np
    rows = (i-1)*nxa + (1:nxa);
    if i == Np
        Qbar(rows,rows) = Ctilde.'*mpcTracker.Qf*Ctilde;
    else
        Qbar(rows,rows) = Ctilde.'*mpcTracker.Q*Ctilde;
    end
end
Rbar = kron(eye(Np), mpcTracker.Rdu);

H = Bbar.'*Qbar*Bbar + Rbar;
H = 0.5*(H + H.') + 1e-9*eye(size(H));
fqp = Bbar.'*Qbar*(Abar*zeta0 + cbar);

Aineq = [];
bineq = [];
LB = [];
UB = [];
if mpcTracker.useIncrementBounds
    LB = repmat(-mpcTracker.deltaUMax, Np, 1);
    UB = repmat( mpcTracker.deltaUMax, Np, 1);
end

if mpcTracker.useInputBounds && use_delta_tau_saturation
    % Constrain the predicted total command:
    %   u_total_i = u_nom_i + c_prev + sum_{j=0}^i Delta c_j
    % inside the same delta_tau limits used by the nonlinear execution.
    S = kron(tril(ones(Np)), eye(nu));
    umaxStack = repmat(delta_tau_limit(:), Np, 1);
    uminStack = -umaxStack;
    cPrevStack = repmat(c_prev(:), Np, 1);
    uNomStack = zeros(nu*Np,1);
    for i = 1:Np
        kk = min(k_now + i - 1, numSamples);
        uNomStack((i-1)*nu + (1:nu)) = u_nom(:,kk);
    end
    Aineq = [ S; -S ];
    bineq = [ umaxStack - uNomStack - cPrevStack; ...
             -uminStack + uNomStack + cPrevStack ];
end

quadprog_available = (exist('quadprog','file') == 2 || exist('quadprog','builtin') == 5);
if mpcTracker.useQuadprog && quadprog_available
    options = optimoptions('quadprog', 'Display','off', ...
        'Algorithm','interior-point-convex', ...
        'OptimalityTolerance',1e-9, 'ConstraintTolerance',1e-9, 'MaxIterations',200);
    [dC,~,exitflag,output] = quadprog(H, fqp, Aineq, bineq, [], [], LB, UB, zeros(nu*Np,1), options);
    if exitflag <= 0 || isempty(dC)
        % Fall back to the unconstrained FONC solution. This keeps the script
        % executable even if quadprog is temporarily infeasible near saturation.
        dC = -H\fqp;
        exitflag = -100;
        output = struct('message','quadprog failed; used regularized unconstrained correction solution');
    end
else
    dC = -H\fqp;
    exitflag = 999;
    output = struct('message','closed-form unconstrained correction solution used');
end

dc0 = dC(1:nu);
c_cmd = c_prev(:) + dc0;
u_cmd = u_nom(:,k_now) + c_cmd;
if use_delta_tau_saturation
    u_cmd = min(max(u_cmd, -delta_tau_limit(:)), delta_tau_limit(:));
    c_cmd = u_cmd - u_nom(:,k_now);
end

info = struct();
info.exitflag = exitflag;
info.Np = Np;
info.dc0 = dc0;
info.du0 = dc0; % alias for older plot/debug code
info.c_prev = c_prev;
info.c_cmd = c_cmd;
info.u_nom_k = u_nom(:,k_now);
info.u_cmd = u_cmd;
info.g0 = gseq{1};
info.predicted_initial_error_norm = norm(e0(1:3));
if isfield(output,'message')
    info.message = output.message;
else
    info.message = '';
end
end

function [Ad, Bd, x_nom_next] = discreteLinearizeAroundNominalLocal( ...
    robot, idx, q_fix_456, x_ref, u_ref, dt, gravity_compensation_enabled, eps_x, eps_u)
% Discrete-time numerical linearization of the RK4 plant map
%   x_next = F_d(x,u)
% around (x_ref,u_ref).  Returning x_nom_next lets the MPC include the
% reference/nominal mismatch defect g_i = x_nom_next - x_ref_next.
nx = numel(x_ref);
nu = numel(u_ref);
Ad = zeros(nx,nx);
Bd = zeros(nx,nu);
x_nom_next = rk4StepReducedPlant(robot, idx, q_fix_456, x_ref, u_ref, dt, gravity_compensation_enabled);
for ii = 1:nx
    dx = zeros(nx,1);
    dx(ii) = eps_x;
    xp = rk4StepReducedPlant(robot, idx, q_fix_456, x_ref + dx, u_ref, dt, gravity_compensation_enabled);
    xm = rk4StepReducedPlant(robot, idx, q_fix_456, x_ref - dx, u_ref, dt, gravity_compensation_enabled);
    Ad(:,ii) = (xp - xm)/(2*eps_x);
end
for jj = 1:nu
    du = zeros(nu,1);
    du(jj) = eps_u;
    xp = rk4StepReducedPlant(robot, idx, q_fix_456, x_ref, u_ref + du, dt, gravity_compensation_enabled);
    xm = rk4StepReducedPlant(robot, idx, q_fix_456, x_ref, u_ref - du, dt, gravity_compensation_enabled);
    Bd(:,jj) = (xp - xm)/(2*eps_u);
end
end

function [Ad, Bd] = zohDiscretizeLocal(Ac, Bc, dt)
nx = size(Ac,1);
nu = size(Bc,2);
M = expm([Ac, Bc; zeros(nu,nx+nu)]*dt);
Ad = M(1:nx,1:nx);
Bd = M(1:nx,nx+1:nx+nu);
end

function [Abar, Bbar, cbar] = buildStackedPredictionMatricesLTVAffineLocal(Aseq, Bseq, cseq, Np)
nx = size(Aseq{1},1);
nu = size(Bseq{1},2);
Abar = zeros(nx*Np, nx);
Bbar = zeros(nx*Np, nu*Np);
cbar = zeros(nx*Np, 1);
Phi = eye(nx);
cprev = zeros(nx,1);
for i = 1:Np
    Ai = Aseq{i};
    Bi = Bseq{i};
    ci = cseq{i};
    Phi = Ai*Phi;
    cprev = Ai*cprev + ci;
    rows = (i-1)*nx + (1:nx);
    Abar(rows,:) = Phi;
    cbar(rows) = cprev;
    for j = 1:i
        G = eye(nx);
        for ell = i:-1:(j+1)
            G = G*Aseq{ell};
        end
        cols = (j-1)*nu + (1:nu);
        Bbar(rows,cols) = G*Bseq{j};
    end
end
end

function lqrDesign = designMidpointLTILQR( ...
    robot, idx, q_fix_456, q_ref, qd_ref, u_nom, t, ...
    Q, R, gravity_compensation_enabled, eps_x, eps_u, use_nominal_feedforward)

% Midpoint LTI-LQR design.
% The nonlinear reduced plant is numerically linearized once at the middle of
% the MPC optimal trajectory. The resulting constant linear model is:
%
%       x_tilde_dot = A_mid*x_tilde + B_mid*u_tilde
%
% with x_tilde = [q-q_ref; qd-qd_ref]. The infinite-horizon LTI-LQR gain is
% obtained from the continuous algebraic Riccati equation:
%
%       A'*P + P*A - P*B*inv(R)*B'*P + Q = 0
%       K = inv(R)*B'*P
%
% In simulation, the tracking command is:
%
%       delta_tau_cmd = u_nom - K*x_tilde
%
% when use_nominal_feedforward=true. This is the standard stabilizing LQR
% tracking form around a moving nominal trajectory. Setting it false gives the
% pure feedback-only form delta_tau_cmd=-K*x_tilde.

numSamples = numel(t);
midIndex = max(1, min(numSamples, round((numSamples+1)/2)));
x_mid = [q_ref(:,midIndex); qd_ref(:,midIndex)];
if use_nominal_feedforward
    u_mid = u_nom(:,midIndex);
else
    u_mid = zeros(3,1);
end

[A_mid, B_mid] = numericalLinearizationReducedPlant( ...
    robot, idx, q_fix_456, x_mid, u_mid, gravity_compensation_enabled, eps_x, eps_u);

% Stabilizability/conditioning diagnostics. These are saved for thesis/debugging.
ctrb_rank = rank(ctrbLocal(A_mid, B_mid), 1e-9);
eig_A_open_loop = eig(A_mid);

% Prefer MATLAB's lqr if available; otherwise use a Hamiltonian CARE solver.
if exist('lqr','file') == 2 || exist('lqr','builtin') == 5
    [K_mid, P_mid, eig_Acl] = lqr(A_mid, B_mid, Q, R);
else
    P_mid = solveCareHamiltonian(A_mid, B_mid, Q, R);
    K_mid = R \ (B_mid.'*P_mid);
    eig_Acl = eig(A_mid - B_mid*K_mid);
end

% If the gain is unrealistically large, soften it automatically. This prevents
% the LQR branch from destroying the nonlinear simulation if the midpoint model
% is poorly conditioned for a particular MPC trajectory.
maxAllowedGainNorm = 120;
gainNorm = norm(K_mid, 'fro');
gainScale = 1.0;
if gainNorm > maxAllowedGainNorm
    gainScale = maxAllowedGainNorm/gainNorm;
    K_mid = gainScale*K_mid;
    eig_Acl = eig(A_mid - B_mid*K_mid);
end

K = repmat(K_mid, 1, 1, numSamples);
A = repmat(A_mid, 1, 1, numSamples);
B = repmat(B_mid, 1, 1, numSamples);
P = repmat(P_mid, 1, 1, numSamples);

lqrDesign = struct();
lqrDesign.type = 'midpoint_lti_lqr';
lqrDesign.midIndex = midIndex;
lqrDesign.midTime = t(midIndex);
lqrDesign.x_mid = x_mid;
lqrDesign.u_mid = u_mid;
lqrDesign.A_mid = A_mid;
lqrDesign.B_mid = B_mid;
lqrDesign.P_mid = P_mid;
lqrDesign.K_mid = K_mid;
lqrDesign.K = K;
lqrDesign.A = A;
lqrDesign.B = B;
lqrDesign.P = P;
lqrDesign.Q = Q;
lqrDesign.R = R;
lqrDesign.use_nominal_feedforward = use_nominal_feedforward;
lqrDesign.ctrb_rank = ctrb_rank;
lqrDesign.eig_A_open_loop = eig_A_open_loop;
lqrDesign.eig_A_closed_loop = eig_Acl;
lqrDesign.gain_fro_norm_before_scaling = gainNorm;
lqrDesign.gain_scale_applied = gainScale;
lqrDesign.description = 'Constant midpoint-linearized continuous-time LTI-LQR tracking gain. K is replicated over time for animation/log compatibility.';
end

function Ctrb = ctrbLocal(A,B)
nx = size(A,1);
Ctrb = zeros(nx, nx*size(B,2));
block = B;
for i = 1:nx
    cols = (i-1)*size(B,2) + (1:size(B,2));
    Ctrb(:,cols) = block;
    block = A*block;
end
end

function P = solveCareHamiltonian(A,B,Q,R)
% Fallback continuous-time CARE solver for:
% A'*P + P*A - P*B*inv(R)*B'*P + Q = 0
nx = size(A,1);
G = B*(R\B.');
H = [A, -G; -Q, -A.'];
[V,D] = eig(H);
lambda = diag(D);
stableIdx = find(real(lambda) < -1e-8);
if numel(stableIdx) < nx
    % If numerical eigenvalues are close to the imaginary axis, select the nx
    % most stable eigenvalues.
    [~, order] = sort(real(lambda), 'ascend');
    stableIdx = order(1:nx);
else
    [~, order] = sort(real(lambda(stableIdx)), 'ascend');
    stableIdx = stableIdx(order(1:nx));
end
Vstable = V(:, stableIdx);
V1 = Vstable(1:nx,:);
V2 = Vstable(nx+1:end,:);
P = real(V2 / V1);
P = 0.5*(P + P.');
end

function [A,B] = numericalLinearizationReducedPlant(robot, idx, q_fix_456, x0, u0, gravity_compensation_enabled, eps_x, eps_u)
nx = numel(x0);
nu = numel(u0);
A = zeros(nx,nx);
B = zeros(nx,nu);
for i = 1:nx
    dx = zeros(nx,1);
    dx(i) = eps_x;
    fp = reducedStateDerivative(robot, idx, q_fix_456, x0 + dx, u0, gravity_compensation_enabled);
    fm = reducedStateDerivative(robot, idx, q_fix_456, x0 - dx, u0, gravity_compensation_enabled);
    A(:,i) = (fp - fm)/(2*eps_x);
end
for j = 1:nu
    du = zeros(nu,1);
    du(j) = eps_u;
    fp = reducedStateDerivative(robot, idx, q_fix_456, x0, u0 + du, gravity_compensation_enabled);
    fm = reducedStateDerivative(robot, idx, q_fix_456, x0, u0 - du, gravity_compensation_enabled);
    B(:,j) = (fp - fm)/(2*eps_u);
end
end

function xnext = rk4StepReducedPlant(robot, idx, q_fix_456, x, delta_tau_cmd, dt, gravity_compensation_enabled)
k1 = reducedStateDerivative(robot, idx, q_fix_456, x,              delta_tau_cmd, gravity_compensation_enabled);
k2 = reducedStateDerivative(robot, idx, q_fix_456, x + 0.5*dt*k1, delta_tau_cmd, gravity_compensation_enabled);
k3 = reducedStateDerivative(robot, idx, q_fix_456, x + 0.5*dt*k2, delta_tau_cmd, gravity_compensation_enabled);
k4 = reducedStateDerivative(robot, idx, q_fix_456, x + dt*k3,     delta_tau_cmd, gravity_compensation_enabled);
xnext = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function xdot = reducedStateDerivative(robot, idx, q_fix_456, x, delta_tau_cmd, gravity_compensation_enabled)
q  = x(1:3);
qd = x(4:6);
q_full  = [q; q_fix_456];
qd_full = [qd; zeros(3,1)];
D3 = reducedMassMatrix(robot, q_full, idx);
C3 = reducedVelocityProduct(robot, q_full, qd_full, idx);
if gravity_compensation_enabled
    qdd = D3 \ (delta_tau_cmd - C3);
else
    G3 = reducedGravity(robot, idx, q_fix_456, q);
    qdd = D3 \ (delta_tau_cmd - C3 - G3);
end
xdot = [qd; qdd];
end

function tau_abs_k = composeAbsoluteTorque(robot, idx, q_fix_456, q, delta_tau_k, gravity_compensation_enabled)
if gravity_compensation_enabled
    G3 = reducedGravity(robot, idx, q_fix_456, q);
    tau_abs_k = G3 + delta_tau_k;
else
    tau_abs_k = delta_tau_k;
end
end

function delta_tau_ff = reducedInverseDynamicsFeedforward(robot, idx, q_fix_456, q_ref, qd_ref, qdd_ref)
q_full_ref  = [q_ref(:);  q_fix_456(:)];
qd_full_ref = [qd_ref(:); zeros(3,1)];
D3 = reducedMassMatrix(robot, q_full_ref, idx);
C3 = reducedVelocityProduct(robot, q_full_ref, qd_full_ref, idx);
delta_tau_ff = D3*qdd_ref(:) + C3;
end

function G3 = reducedGravity(robot, idx, q_fix_456, q)
q_full = [q(:); q_fix_456(:)];
G6 = gravityTorque(robot, q_full);
G3 = G6(idx);
end

function D3 = reducedMassMatrix(robot, q_full, idx)
D6 = massMatrix(robot, q_full);
D3 = D6(idx, idx);
end

function C3 = reducedVelocityProduct(robot, q_full, qd_full, idx)
if exist('velocityProduct','file') == 2 || exist('velocityProduct','builtin') == 5
    C6 = velocityProduct(robot, q_full, qd_full);
    C3 = C6(idx);
    return;
end
gravityBackup = robot.Gravity;
cleanupObj = onCleanup(@() restoreRobotGravity(robot, gravityBackup)); %#ok<NASGU>
robot.Gravity = [0 0 0];
tauVel6 = inverseDynamics(robot, q_full, qd_full, zeros(size(qd_full)));
C3 = tauVel6(idx);
end

function restoreRobotGravity(robot, gravityVector)
robot.Gravity = gravityVector;
end

function endpointSeries = reducedEndpointSeries(robot, bodyName, q_traj, q_fix_456)
numSamples = size(q_traj,2);
endpointSeries = zeros(3, numSamples);
for k = 1:numSamples
    endpointSeries(:,k) = reducedFKPosition(robot, bodyName, q_traj(:,k), q_fix_456);
end
end

function p = reducedFKPosition(robot, bodyName, q123, q_fix_456)
q_full = [q123(:); q_fix_456(:)];
T = getTransform(robot, q_full, bodyName);
p = T(1:3,4);
end

function settlingTime = computeSettlingTime(err_norm, t, threshold, holdSamples)
settlingTime = NaN;
numSamples = numel(err_norm);
for k = 1:(numSamples - holdSamples + 1)
    window = err_norm(k:k+holdSamples-1);
    if all(window <= threshold)
        settlingTime = t(k);
        return;
    end
end
end

function checkTrajectoryDimensions(q_ref, qd_ref, qdd_ref, expectedSamples, methodName)
assert(isequal(size(q_ref),   [3, expectedSamples]), '%s: q_ref must be 3-by-%d.', methodName, expectedSamples);
assert(isequal(size(qd_ref),  [3, expectedSamples]), '%s: qd_ref must be 3-by-%d.', methodName, expectedSamples);
assert(isequal(size(qdd_ref), [3, expectedSamples]), '%s: qdd_ref must be 3-by-%d.', methodName, expectedSamples);
assert(all(isfinite(q_ref), 'all'),   '%s: q_ref contains non-finite values.', methodName);
assert(all(isfinite(qd_ref), 'all'),  '%s: qd_ref contains non-finite values.', methodName);
assert(all(isfinite(qdd_ref), 'all'), '%s: qdd_ref contains non-finite values.', methodName);
end

function O = computePlanningVsExecutionOptimality(q_plan, qd_plan, u_plan, q_exec, qd_exec, u_exec, qg, qdg, t, Qgoal, Rgoal)
N = numel(t) - 1;
Jstage_plan = zeros(1,N);
Jstage_exec_goal = zeros(1,N);
traj_deviation_norm = zeros(1,N+1);
state_deviation_norm = zeros(1,N+1);
for k = 1:N+1
    dq = q_exec(:,k) - q_plan(:,k);
    dqd = qd_exec(:,k) - qd_plan(:,k);
    traj_deviation_norm(k) = norm(dq);
    state_deviation_norm(k) = norm([dq; dqd]);
end
for k = 1:N
    x_plan = [q_plan(:,k) - qg(:); qd_plan(:,k) - qdg(:)];
    x_exec = [q_exec(:,k) - qg(:); qd_exec(:,k) - qdg(:)];
    uk_plan = u_plan(:,k);
    uk_exec = u_exec(:,k);
    Jstage_plan(k) = 0.5*(x_plan.'*Qgoal*x_plan) + 0.5*(uk_plan.'*Rgoal*uk_plan);
    Jstage_exec_goal(k) = 0.5*(x_exec.'*Qgoal*x_exec) + 0.5*(uk_exec.'*Rgoal*uk_exec);
end
Jcum_plan = cumsum(Jstage_plan);
Jcum_exec_goal = cumsum(Jstage_exec_goal);
final_planned_goal_cost = Jcum_plan(end);
final_executed_goal_cost = Jcum_exec_goal(end);
absolute_optimality_loss = final_executed_goal_cost - final_planned_goal_cost;
relative_optimality_loss = absolute_optimality_loss / max(abs(final_planned_goal_cost), eps);
control_deviation = u_exec - u_plan;
control_deviation_norm = vecnorm(control_deviation,2,1);
control_deviation_energy = cumtrapz(t, sum(control_deviation.^2,1));
O = struct();
O.Jstage_plan = Jstage_plan;
O.Jstage_exec_goal = Jstage_exec_goal;
O.Jcum_plan = Jcum_plan;
O.Jcum_exec_goal = Jcum_exec_goal;
O.final_planned_goal_cost = final_planned_goal_cost;
O.final_executed_goal_cost = final_executed_goal_cost;
O.absolute_optimality_loss = absolute_optimality_loss;
O.relative_optimality_loss = relative_optimality_loss;
O.traj_deviation_norm = traj_deviation_norm;
O.state_deviation_norm = state_deviation_norm;
O.control_deviation = control_deviation;
O.control_deviation_norm = control_deviation_norm;
O.control_deviation_energy = control_deviation_energy;
end