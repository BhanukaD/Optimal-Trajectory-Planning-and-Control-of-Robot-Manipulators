%% Reduced UR5 Trajectory Generation and Controller Comparison
% This script generates a reduced 3-DoF UR5 point-to-point trajectory using
% online receding-horizon LTV MPC and compares it with nonlinear PID and LQR
% goal-regulation baselines. The generated data are saved for tracking,
% visualization, and performance analysis.

clear; clc; close all;

%% Simulation and Controller Settings
% Sampling time and total simulation duration.
dt = 0.02;                 % [s]
Tf = 2.00;                 % [s]
t  = (0:dt:Tf).';          % [s], N+1 samples
N  = numel(t) - 1;         % number of control intervals

% Active joint subset and fixed passive joint values.
idx = [1 2 3];
q_fix_456 = [0; 0; 0];

% Initial and goal configurations for the point-to-point motion.
q0   = [-1.35; -2.05;  2.35];   % [rad] lower/back corner posture
qg   = [ 1.15; -0.65;  0.75];   % [rad] upper/front swept posture

qd0  = zeros(3,1);              % [rad/s]
qdg  = zeros(3,1);              % [rad/s]
qdd0 = zeros(3,1);              % [rad/s^2]
qddg = zeros(3,1);              % [rad/s^2]

% State and input weighting matrices for the LTV MPC problem.
% State vector: x = [q - qg; qd - qdg].
Qq = 150*eye(3);
Qd =  80*eye(3);
Qx = blkdiag(Qq, Qd);

% Input weighting matrix for the MPC quadratic program.
Ru = 35*eye(3);

% Terminal state weighting matrix.
Qf = 50*Qx;

% Prediction and control horizon settings.
Np_max = 50;                       % prediction horizon in samples
mpc_control_horizon = 1;           % receding horizon: apply first action only

% Number of local LTV refinement iterations per MPC sample.
mpc_ltv_refinement_iterations = 1;

% Input bounds for the torque increment command.
use_delta_tau_bounds = true;
delta_tau_min = -60*ones(3,1);     % [Nm]
delta_tau_max =  60*ones(3,1);     % [Nm]

% Terminal equality option for the final local MPC horizon.
enforce_terminal_equality_at_final_horizon = true;
terminal_equality_tolerance_for_acceptance = 5e-3;

% Numerical perturbations used for finite-difference linearization.
eps_x = 1e-6;
eps_u = 1e-6;

% Trapezoidal profile design parameter.
trap_accel_fraction = 0.20;

% Plot formatting parameters.
lineWidthMain = 1.8;
fontSizeAxes  = 11;

%% Load Reduced UR5 Model
assert(exist('loadrobot','file') == 2 || exist('loadrobot','builtin') == 5, ...
    'Robotics System Toolbox is required: loadrobot not found.');

robot = loadrobot("universalUR5", ...
    "DataFormat","column", ...
    "Gravity",[0 0 -9.81]);

reducedBody = pickReducedUR5EndpointBody(robot);

%% Input Validation
validateattributes(q0,       {'double'},{'column','numel',3,'finite','real'});
validateattributes(qg,       {'double'},{'column','numel',3,'finite','real'});
validateattributes(qd0,      {'double'},{'column','numel',3,'finite','real'});
validateattributes(qdg,      {'double'},{'column','numel',3,'finite','real'});
validateattributes(q_fix_456,{'double'},{'column','numel',3,'finite','real'});
assert(N >= 1, 'The horizon must contain at least one interval.');
assert(Np_max >= 1, 'Np_max must be at least 1.');

start_endpoint = reducedFKPosition(robot, reducedBody, q0, q_fix_456);
goal_endpoint  = reducedFKPosition(robot, reducedBody, qg, q_fix_456);

fprintf('------------------------------------------------------------\n');
fprintf('Online LTV MPC reduced-UR5 planner setup\n');
fprintf('Reduced endpoint body : %s\n', reducedBody);
fprintf('dt = %.4f s, Tf = %.4f s, N = %d, Np_max = %d\n', dt, Tf, N, Np_max);
fprintf('q0 = [% .4f % .4f % .4f]^T rad\n', q0(1), q0(2), q0(3));
fprintf('qg = [% .4f % .4f % .4f]^T rad\n', qg(1), qg(2), qg(3));
fprintf('Start endpoint = [% .4f % .4f % .4f]^T m\n', start_endpoint(1), start_endpoint(2), start_endpoint(3));
fprintf('Goal endpoint  = [% .4f % .4f % .4f]^T m\n', goal_endpoint(1), goal_endpoint(2), goal_endpoint(3));
fprintf('------------------------------------------------------------\n');

%% Nonlinear Error Dynamics
f_error = @(x,u) reducedErrorDynamics(robot, idx, q_fix_456, qg, qdg, x, u);

% % ONLINE RECEDING-HORIZON LTV MPC PLANNER
x0 = [q0 - qg; qd0 - qdg];
x_now = x0;

q_ref_opt   = zeros(3, N+1);
qd_ref_opt  = zeros(3, N+1);
x_ref_opt   = zeros(6, N+1);
u_ref_opt   = zeros(3, N);
tau_abs_opt = zeros(3, N);

stage_cost_opt      = zeros(N,1);
cumulative_cost_opt = zeros(N,1);
mpc_exitflag        = zeros(N,1);
mpc_horizon_used    = zeros(N,1);
mpc_terminal_eq_used = false(N,1);

% Warm-start / nominal input sequence. It is shifted after every MPC solve.
U_warm = zeros(3*Np_max,1);

quadprog_available = exist('quadprog','file') == 2 || exist('quadprog','builtin') == 5;
if quadprog_available
    qp_options = optimoptions('quadprog', ...
        'Display','off', ...
        'Algorithm','interior-point-convex', ...
        'OptimalityTolerance',1e-8, ...
        'ConstraintTolerance',1e-8, ...
        'MaxIterations',200);
else
    qp_options = [];
    warning('quadprog not found. Falling back to unconstrained U = -H\\f.');
end

for k = 1:N
% Log current state before applying control at interval k.
    x_ref_opt(:,k) = x_now;
    q_ref_opt(:,k) = x_now(1:3) + qg;
    qd_ref_opt(:,k) = x_now(4:6) + qdg;

% Shrink horizon near the end, exactly as a finite-duration MPC problem.
    Np = min(Np_max, N-k+1);
    mpc_horizon_used(k) = Np;

% Warm-start controls for this local horizon.
    U_guess = U_warm(1:3*Np);

% Re-linearize / re-solve, optionally more than once.
    for refine = 1:mpc_ltv_refinement_iterations
% Build nominal predicted rollout from current measured/planned state.
        [x_nom, u_nom] = rolloutNominalTrajectory(f_error, x_now, U_guess, dt, Np);

% Linearize the nonlinear plant at every nominal predicted stage.
        Aseq = cell(1,Np);
        Bseq = cell(1,Np);
        cseq = cell(1,Np);
        for i = 1:Np
            xi = x_nom(:,i);
            ui = u_nom(:,i);
            [Ac_i, Bc_i] = numericalJacobians(f_error, xi, ui, eps_x, eps_u);
            [Ad_i, Bd_i] = zohDiscretize(Ac_i, Bc_i, dt);
            Aseq{i} = Ad_i;
            Bseq{i} = Bd_i;
            cseq{i} = x_nom(:,i+1) - Ad_i*x_nom(:,i) - Bd_i*ui;
        end

% Stacked LTV affine prediction:
% X = Abar*x_now + Bbar*U + cbar
        [Abar, Bbar, cbar] = buildStackedPredictionMatricesLTVAffine(Aseq, Bseq, cseq, Np);
        Qbar = blkdiag(kron(eye(max(Np-1,0)), Qx), Qf);
        Rbar = kron(eye(Np), Ru);

        H = Bbar.'*Qbar*Bbar + Rbar;
        H = 0.5*(H + H.') + 1e-9*eye(size(H));
        fqp = Bbar.'*Qbar*(Abar*x_now + cbar);

% Bounds for delta_tau over prediction horizon.
        if use_delta_tau_bounds
            LB = repmat(delta_tau_min, Np, 1);
            UB = repmat(delta_tau_max, Np, 1);
        else
            LB = [];
            UB = [];
        end

% Optional final-horizon terminal equality in the local LTV prediction.
        Aeq = [];
        beq = [];
        useTerminalEqThisSolve = false;
        if enforce_terminal_equality_at_final_horizon && (k + Np - 1 == N)
            n = 6;
            Aeq = Bbar((Np-1)*n+1:Np*n, :);
            beq = -(Abar((Np-1)*n+1:Np*n,:)*x_now + cbar((Np-1)*n+1:Np*n));
            useTerminalEqThisSolve = true;
        end

        if quadprog_available
            [U_sol,~,exitflag] = quadprog(H, fqp, [], [], Aeq, beq, LB, UB, U_guess, qp_options);
            if exitflag <= 0 && useTerminalEqThisSolve
% is too restrictive or numerically infeasible.
                [U_sol,~,exitflag] = quadprog(H, fqp, [], [], [], [], LB, UB, U_guess, qp_options);
                useTerminalEqThisSolve = false;
            end
            if exitflag <= 0 || isempty(U_sol)
                warning('quadprog failed at k=%d with exitflag=%d. Using unconstrained regularized solution.', k, exitflag);
                U_sol = -H\fqp;
            end
        else
            U_sol = -H\fqp;
            exitflag = 999;
        end

        U_guess = U_sol;
    end

    mpc_exitflag(k) = exitflag;
    mpc_terminal_eq_used(k) = useTerminalEqThisSolve;

% Receding-horizon action: apply only the first control signal.
    delta_tau = U_guess(1:3);
    if use_delta_tau_bounds
        delta_tau = min(max(delta_tau, delta_tau_min), delta_tau_max);
    end
    u_ref_opt(:,k) = delta_tau;

% Absolute torque log, using the gravity-compensated input convention.
    q_now = x_now(1:3) + qg;
    q_full = [q_now; q_fix_456];
    G6 = gravityTorque(robot, q_full);
    tau_abs_opt(:,k) = G6(idx) + delta_tau;

% Stage cost in the same quadratic form as the planner.
    stage_cost_opt(k) = 0.5*(x_now.'*Qx*x_now) + 0.5*(delta_tau.'*Ru*delta_tau);
    if k == 1
        cumulative_cost_opt(k) = stage_cost_opt(k);
    else
        cumulative_cost_opt(k) = cumulative_cost_opt(k-1) + stage_cost_opt(k);
    end

% Propagate the actual nonlinear reduced plant by one sample.
    x_now = rk4Step(f_error, x_now, delta_tau, dt);

% Shift warm-start sequence for the next receding horizon solve.
    U_warm = shiftWarmStart(U_guess, Np, Np_max, 3);
end

% Log terminal state sample.
x_ref_opt(:,N+1) = x_now;
q_ref_opt(:,N+1) = x_now(1:3) + qg;
qd_ref_opt(:,N+1) = x_now(4:6) + qdg;
qdd_ref_opt = estimateAccelerationFromVelocity(qd_ref_opt, dt);
J_plan_opt = cumulative_cost_opt(end);
terminal_error_mpc = x_ref_opt(:,end);
terminal_error_norm_mpc = norm(terminal_error_mpc);

% reference endpoint exactly qg,qdg for clean downstream plotting/tracking.
if terminal_error_norm_mpc < terminal_equality_tolerance_for_acceptance
    q_ref_opt(:,end) = qg;
    qd_ref_opt(:,end) = qdg;
    x_ref_opt(:,end) = zeros(6,1);
    qdd_ref_opt = estimateAccelerationFromVelocity(qd_ref_opt, dt);
end

fprintf('Online LTV MPC completed. Final ||x_N|| = %.6e, J = %.8g\n', ...
    terminal_error_norm_mpc, J_plan_opt);

% % COMPARISON PLANNERS: PID / LQR GOAL REGULATION
% These are nonlinear closed-loop executions from the same start state to the
% same goal state. They are not kinematic baselines and they do not track the
% MPC trajectory. They directly regulate q -> qg and qd -> qdg using the same
% reduced nonlinear UR5 dynamics and the same delta_tau convention as MPC.

% because the plant is propagated online through the nonlinear reduced model.
pid_Kp = diag([140 140 140]);
pid_Ki = diag([  8   8   8]);
pid_Kd = diag([ 35  35  35]);
pid_integral_limit = 0.75*ones(3,1);  % anti-windup clamp [rad*s]

pid_use_delta_tau_bounds = use_delta_tau_bounds;
lqr_use_delta_tau_bounds = use_delta_tau_bounds;

[q_ref_pid, qd_ref_pid, u_ref_pid, pid_info] = generatePIDGoalRegulationTrajectory( ...
    f_error, robot, idx, q_fix_456, q0, qd0, qg, qdg, t, ...
    pid_Kp, pid_Ki, pid_Kd, pid_integral_limit, ...
    pid_use_delta_tau_bounds, delta_tau_min, delta_tau_max);
qdd_ref_pid = estimateAccelerationFromVelocity(qd_ref_pid, dt);

% LQR gain designed about the goal equilibrium x = 0, delta_tau = 0.
% This gives a true goal-regulation LQR baseline rather than a tracker of MPC.
lqr_Q = blkdiag(250*eye(3), 25*eye(3));
lqr_R = 0.35*eye(3);
[K_lqr_goal, lqr_design_info] = designGoalLQRGain(f_error, dt, eps_x, eps_u, lqr_Q, lqr_R);

[q_ref_lqr, qd_ref_lqr, u_ref_lqr, lqr_info] = generateLQRGoalRegulationTrajectory( ...
    f_error, robot, idx, q_fix_456, q0, qd0, qg, qdg, t, K_lqr_goal, ...
    lqr_use_delta_tau_bounds, delta_tau_min, delta_tau_max);
qdd_ref_lqr = estimateAccelerationFromVelocity(qd_ref_lqr, dt);

fprintf('PID goal regulation completed. Final ||x_N|| = %.6e\n', pid_info.final_error_norm);
fprintf('LQR goal regulation completed. Final ||x_N|| = %.6e\n', lqr_info.final_error_norm);

%% Consistency Checks
assert(isequal(size(q_ref_opt),     [3, N+1]));
assert(isequal(size(qd_ref_opt),    [3, N+1]));
assert(isequal(size(qdd_ref_opt),   [3, N+1]));
assert(isequal(size(u_ref_opt),     [3, N]));

assert(isequal(size(q_ref_pid),     [3, N+1]));
assert(isequal(size(qd_ref_pid),    [3, N+1]));
assert(isequal(size(qdd_ref_pid),   [3, N+1]));
assert(isequal(size(u_ref_pid),     [3, N]));

assert(isequal(size(q_ref_lqr),     [3, N+1]));
assert(isequal(size(qd_ref_lqr),    [3, N+1]));
assert(isequal(size(qdd_ref_lqr),   [3, N+1]));
assert(isequal(size(u_ref_lqr),     [3, N]));

%% Endpoint Path Computation
path_opt = reducedPathFromTrajectory(robot, reducedBody, q_ref_opt, q_fix_456);
path_pid = reducedPathFromTrajectory(robot, reducedBody, q_ref_pid, q_fix_456);
path_lqr = reducedPathFromTrajectory(robot, reducedBody, q_ref_lqr, q_fix_456);

%% Output Data Structure
planner_data = struct();

planner_data.common.t              = t;
planner_data.common.dt             = dt;
planner_data.common.Tf             = Tf;
planner_data.common.N              = N;
planner_data.common.q0             = q0;
planner_data.common.qg             = qg;
planner_data.common.qd0            = qd0;
planner_data.common.qdg            = qdg;
planner_data.common.idx            = idx;
planner_data.common.q_fix_456      = q_fix_456;
planner_data.common.reducedBody    = reducedBody;
planner_data.common.start_endpoint = start_endpoint;
planner_data.common.goal_endpoint  = goal_endpoint;

planner_data.opt.q_ref             = q_ref_opt;
planner_data.opt.qd_ref            = qd_ref_opt;
planner_data.opt.qdd_ref           = qdd_ref_opt;
planner_data.opt.u_ref             = u_ref_opt;
planner_data.opt.x_ref             = x_ref_opt;
planner_data.opt.tau_abs_ref       = tau_abs_opt;
planner_data.opt.Qx                = Qx;
planner_data.opt.Ru                = Ru;
planner_data.opt.Qf                = Qf;
planner_data.opt.J_plan            = J_plan_opt;
planner_data.opt.stage_cost        = stage_cost_opt;
planner_data.opt.cumulative_cost   = cumulative_cost_opt;
planner_data.opt.endpoint_path     = path_opt;

planner_data.opt.algorithm.type = 'online_receding_horizon_LTV_MPC';
planner_data.opt.algorithm.description = ...
    'At each sample, linearize the reduced nonlinear UR5 plant along the MPC predicted nominal trajectory, solve an affine LTV QP, apply the first control, and repeat.';
planner_data.opt.algorithm.Np_max = Np_max;
planner_data.opt.algorithm.control_horizon = mpc_control_horizon;
planner_data.opt.algorithm.ltv_refinement_iterations = mpc_ltv_refinement_iterations;
planner_data.opt.algorithm.uses_quadprog = quadprog_available;
planner_data.opt.algorithm.use_delta_tau_bounds = use_delta_tau_bounds;
planner_data.opt.algorithm.delta_tau_min = delta_tau_min;
planner_data.opt.algorithm.delta_tau_max = delta_tau_max;
planner_data.opt.algorithm.enforce_terminal_equality_at_final_horizon = enforce_terminal_equality_at_final_horizon;
planner_data.opt.algorithm.exitflag = mpc_exitflag;
planner_data.opt.algorithm.horizon_used = mpc_horizon_used;
planner_data.opt.algorithm.terminal_eq_used = mpc_terminal_eq_used;
planner_data.opt.algorithm.final_x = terminal_error_mpc;
planner_data.opt.algorithm.final_x_norm = terminal_error_norm_mpc;
planner_data.opt.algorithm.input_convention = 'tau_abs = G(q) + delta_tau; D(q)qdd = delta_tau - C(q,qd)';

planner_data.pid.q_ref           = q_ref_pid;
planner_data.pid.qd_ref          = qd_ref_pid;
planner_data.pid.qdd_ref         = qdd_ref_pid;
planner_data.pid.u_ref           = u_ref_pid;
planner_data.pid.endpoint_path   = path_pid;
planner_data.pid.Kp              = pid_Kp;
planner_data.pid.Ki              = pid_Ki;
planner_data.pid.Kd              = pid_Kd;
planner_data.pid.integral_limit  = pid_integral_limit;
planner_data.pid.info            = pid_info;
planner_data.pid.algorithm.type  = 'nonlinear_PID_goal_regulation';
planner_data.pid.algorithm.description = ...
    'Closed-loop PID directly regulates the reduced nonlinear UR5 from q0,qd0 to qg,qdg using delta_tau input; it does not track the MPC trajectory.';
planner_data.pid.algorithm.use_delta_tau_bounds = pid_use_delta_tau_bounds;

planner_data.lqr.q_ref           = q_ref_lqr;
planner_data.lqr.qd_ref          = qd_ref_lqr;
planner_data.lqr.qdd_ref         = qdd_ref_lqr;
planner_data.lqr.u_ref           = u_ref_lqr;
planner_data.lqr.endpoint_path   = path_lqr;
planner_data.lqr.K               = K_lqr_goal;
planner_data.lqr.Q               = lqr_Q;
planner_data.lqr.R               = lqr_R;
planner_data.lqr.design_info     = lqr_design_info;
planner_data.lqr.info            = lqr_info;
planner_data.lqr.algorithm.type  = 'discrete_LQR_goal_regulation';
planner_data.lqr.algorithm.description = ...
    'LQR gain is designed from the linearized reduced nonlinear UR5 dynamics at the goal equilibrium and directly regulates q0,qd0 to qg,qdg.';
planner_data.lqr.algorithm.use_delta_tau_bounds = lqr_use_delta_tau_bounds;

% % Add MPC/PID/LQR comparison metrics for cost/control/power/energy plots

plannerNames = {'Online LTV MPC','PID','LQR'};
numPlanners  = 3;

q_all   = {q_ref_opt,   q_ref_pid,   q_ref_lqr};
qd_all  = {qd_ref_opt,  qd_ref_pid,  qd_ref_lqr};
qdd_all = {qdd_ref_opt, qdd_ref_pid, qdd_ref_lqr};

x_all = cell(numPlanners,1);
u_all = cell(numPlanners,1);
state_cost_all = cell(numPlanners,1);
control_cost_all = cell(numPlanners,1);
stage_cost_all = cell(numPlanners,1);
cumulative_cost_all = cell(numPlanners,1);
cost_to_go_all = cell(numPlanners,1);
joint_power_all = cell(numPlanners,1);
total_power_all = cell(numPlanners,1);
abs_total_power_all = cell(numPlanners,1);
cum_energy_all = cell(numPlanners,1);

total_cost_all = zeros(numPlanners,1);
final_energy_all = zeros(numPlanners,1);
control_rms_all = zeros(numPlanners,3);
peak_power_all = zeros(numPlanners,1);

tau_abs_all = cell(numPlanners,1);
joint_power_actual_all = cell(numPlanners,1);
total_power_actual_all = cell(numPlanners,1);
abs_total_power_actual_all = cell(numPlanners,1);
cum_energy_actual_all = cell(numPlanners,1);
final_energy_actual_all = zeros(numPlanners,1);
peak_power_actual_all = zeros(numPlanners,1);
endpoint_error_all = zeros(numPlanners,1);
path_length_all = zeros(numPlanners,1);

for p = 1:numPlanners
    q_now   = q_all{p};
    qd_now  = qd_all{p};
    qdd_now = qdd_all{p};

    x_now = [q_now - qg; qd_now - qdg];
    x_all{p} = x_now;

    if p == 1
        u_now = u_ref_opt;
    elseif p == 2
        u_now = u_ref_pid;
    else
        u_now = u_ref_lqr;
    end
    u_all{p} = u_now;

    state_cost_k = zeros(N,1);
    control_cost_k = zeros(N,1);
    stage_cost_k = zeros(N,1);
    cumulative_cost_k = zeros(N,1);
    for k = 1:N
        xk = x_now(:,k+1);
        uk = u_now(:,k);
        state_cost_k(k)   = 0.5*(xk.'*Qx*xk);
        control_cost_k(k) = 0.5*(uk.'*Ru*uk);
        stage_cost_k(k)   = state_cost_k(k) + control_cost_k(k);
        if k == 1
            cumulative_cost_k(k) = 0.5*(x_now(:,1).'*Qx*x_now(:,1)) + stage_cost_k(k);
        else
            cumulative_cost_k(k) = cumulative_cost_k(k-1) + stage_cost_k(k);
        end
    end

    state_cost_all{p} = state_cost_k;
    control_cost_all{p} = control_cost_k;
    stage_cost_all{p} = stage_cost_k;
    cumulative_cost_all{p} = cumulative_cost_k;

    cost_to_go_k = zeros(N,1);
    for k = 1:N
        if k == 1
            cost_to_go_k(k) = cumulative_cost_k(end);
        else
            cost_to_go_k(k) = cumulative_cost_k(end) - cumulative_cost_k(k-1);
        end
    end
    cost_to_go_all{p} = cost_to_go_k;
    total_cost_all(p) = cumulative_cost_k(end);

    joint_power_k = u_now .* qd_now(:,1:N);
    total_power_k = sum(joint_power_k, 1).';
    abs_total_power_k = sum(abs(joint_power_k), 1).';
    joint_power_all{p} = joint_power_k;
    total_power_all{p} = total_power_k;
    abs_total_power_all{p} = abs_total_power_k;
    cum_energy_k = cumtrapz(t(1:N), abs_total_power_k);
    cum_energy_all{p} = cum_energy_k;
    final_energy_all(p) = cum_energy_k(end);
    control_rms_all(p,:) = sqrt(mean(u_now.^2, 2)).';
    peak_power_all(p) = max(abs_total_power_k);

    tau_abs_now = zeros(6,N);
    joint_power_actual_now = zeros(6,N);
    total_power_actual_now = zeros(N,1);
    abs_total_power_actual_now = zeros(N,1);
    for k = 1:N
        q_full   = [q_now(:,k);   q_fix_456];
        qd_full  = [qd_now(:,k);  zeros(3,1)];
        qdd_full = [qdd_now(:,k); zeros(3,1)];
        tau_abs_k = inverseDynamics(robot, q_full, qd_full, qdd_full);
        joint_power_actual_k = tau_abs_k .* qd_full;
        tau_abs_now(:,k) = tau_abs_k;
        joint_power_actual_now(:,k) = joint_power_actual_k;
        total_power_actual_now(k) = sum(joint_power_actual_k);
        abs_total_power_actual_now(k) = sum(abs(joint_power_actual_k));
    end
    cum_energy_actual_k = cumtrapz(t(1:N), abs_total_power_actual_now);
    tau_abs_all{p} = tau_abs_now;
    joint_power_actual_all{p} = joint_power_actual_now;
    total_power_actual_all{p} = total_power_actual_now;
    abs_total_power_actual_all{p} = abs_total_power_actual_now;
    cum_energy_actual_all{p} = cum_energy_actual_k;
    final_energy_actual_all(p) = cum_energy_actual_k(end);
    peak_power_actual_all(p) = max(abs_total_power_actual_now);

    endpoint_error_all(p) = norm(q_now(:,end)-qg) + norm(qd_now(:,end)-qdg);
    path_now = reducedPathFromTrajectory(robot, reducedBody, q_now, q_fix_456);
    path_length_all(p) = sum(sqrt(sum(diff(path_now,1,2).^2,1)));
end

planner_data.compare.names = plannerNames;
planner_data.compare.total_cost = total_cost_all;
planner_data.compare.final_energy = final_energy_all;
planner_data.compare.peak_power = peak_power_all;
planner_data.compare.control_rms = control_rms_all;
planner_data.compare.actual_final_energy = final_energy_actual_all;
planner_data.compare.actual_peak_power = peak_power_actual_all;
planner_data.compare.endpoint_error = endpoint_error_all;
planner_data.compare.endpoint_path_length = path_length_all;
planner_data.compare.stage_cost = stage_cost_all;
planner_data.compare.cumulative_cost = cumulative_cost_all;
planner_data.compare.cost_to_go = cost_to_go_all;
planner_data.compare.delta_tau = u_all;
planner_data.compare.tau_abs = tau_abs_all;
planner_data.compare.total_power = total_power_all;
planner_data.compare.abs_total_power = abs_total_power_all;
planner_data.compare.cumulative_energy = cum_energy_all;
planner_data.compare.total_power_actual = total_power_actual_all;
planner_data.compare.abs_total_power_actual = abs_total_power_actual_all;
planner_data.compare.cumulative_energy_actual = cum_energy_actual_all;

planner_data.opt.delta_tau = u_all{1};
planner_data.opt.state_cost = state_cost_all{1};
planner_data.opt.control_cost = control_cost_all{1};
planner_data.opt.stage_cost = stage_cost_all{1};
planner_data.opt.cumulative_cost = cumulative_cost_all{1};
planner_data.opt.cost_to_go = cost_to_go_all{1};
planner_data.opt.joint_power = joint_power_all{1};
planner_data.opt.total_power = total_power_all{1};
planner_data.opt.abs_total_power = abs_total_power_all{1};
planner_data.opt.cumulative_energy = cum_energy_all{1};
planner_data.opt.tau_abs_actual = tau_abs_all{1};
planner_data.opt.actual_joint_power = joint_power_actual_all{1};
planner_data.opt.actual_total_power = total_power_actual_all{1};
planner_data.opt.actual_abs_total_power = abs_total_power_actual_all{1};
planner_data.opt.actual_cumulative_energy = cum_energy_actual_all{1};

planner_data.pid.delta_tau = u_all{2};
planner_data.pid.state_cost = state_cost_all{2};
planner_data.pid.control_cost = control_cost_all{2};
planner_data.pid.stage_cost = stage_cost_all{2};
planner_data.pid.cumulative_cost = cumulative_cost_all{2};
planner_data.pid.cost_to_go = cost_to_go_all{2};
planner_data.pid.joint_power = joint_power_all{2};
planner_data.pid.total_power = total_power_all{2};
planner_data.pid.abs_total_power = abs_total_power_all{2};
planner_data.pid.cumulative_energy = cum_energy_all{2};
planner_data.pid.tau_abs_actual = tau_abs_all{2};
planner_data.pid.actual_joint_power = joint_power_actual_all{2};
planner_data.pid.actual_total_power = total_power_actual_all{2};
planner_data.pid.actual_abs_total_power = abs_total_power_actual_all{2};
planner_data.pid.actual_cumulative_energy = cum_energy_actual_all{2};

planner_data.lqr.delta_tau = u_all{3};
planner_data.lqr.state_cost = state_cost_all{3};
planner_data.lqr.control_cost = control_cost_all{3};
planner_data.lqr.stage_cost = stage_cost_all{3};
planner_data.lqr.cumulative_cost = cumulative_cost_all{3};
planner_data.lqr.cost_to_go = cost_to_go_all{3};
planner_data.lqr.joint_power = joint_power_all{3};
planner_data.lqr.total_power = total_power_all{3};
planner_data.lqr.abs_total_power = abs_total_power_all{3};
planner_data.lqr.cumulative_energy = cum_energy_all{3};
planner_data.lqr.tau_abs_actual = tau_abs_all{3};
planner_data.lqr.actual_joint_power = joint_power_actual_all{3};
planner_data.lqr.actual_total_power = total_power_actual_all{3};
planner_data.lqr.actual_abs_total_power = abs_total_power_actual_all{3};
planner_data.lqr.actual_cumulative_energy = cum_energy_actual_all{3};

%% Save Results
save('planner_data.mat', 'planner_data');

% % PUBLICATION-QUALITY RESULTS, TABLES, AND PLOTS
% constraints, or trajectory-generation method.

resultsRoot = fullfile(pwd, 'publication_results_mpc_pid_lqr');
figDir      = fullfile(resultsRoot, 'figures');
tableDir    = fullfile(resultsRoot, 'tables');
dataDir     = fullfile(resultsRoot, 'data');
ensureDir(resultsRoot); ensureDir(figDir); ensureDir(tableDir); ensureDir(dataDir);

save(fullfile(dataDir, 'planner_data_raw.mat'), 'planner_data');

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
pub.markerSize    = 58;
pub.dpi           = 600;
pub.formats       = {'png','pdf','fig'};

styles      = {'-','--','-.'};
labelsLatex = {'LTV-MPC','PID','LQR'};
jointLabels = {'$q_1$','$q_2$','$q_3$'};
paths_all   = {path_opt, path_pid, path_lqr};

%% Torque decomposition
% Torque components are computed from the already generated trajectories only.
% They are used for reporting/plotting and do not feed back into the planners.
tau_components_all = cell(numPlanners,1);
for p = 1:numPlanners
    tau_components_all{p} = computeTorqueComponentsReduced(robot, idx, q_fix_456, q_all{p}, qd_all{p}, qdd_all{p}, N);
end
planner_data.opt.torque_components = tau_components_all{1};
planner_data.pid.torque_components = tau_components_all{2};
planner_data.lqr.torque_components = tau_components_all{3};

%% Summary metrics table
max_abs_acc_all          = zeros(numPlanners,3);
mean_abs_jerk_all        = zeros(numPlanners,3);
rms_delta_tau_all        = zeros(numPlanners,3);
peak_delta_tau_all       = zeros(numPlanners,3);
total_variation_tau_all  = zeros(numPlanners,1);
final_q_error_all        = zeros(numPlanners,1);
final_qd_error_all       = zeros(numPlanners,1);
max_overshoot_all        = zeros(numPlanners,1);
settling_time_all        = zeros(numPlanners,1);
constraint_violation_all = zeros(numPlanners,1);
endpoint_straightness_all = zeros(numPlanners,1);
endpoint_terminal_error_all = zeros(numPlanners,1);

for p = 1:numPlanners
    q_now   = q_all{p};
    qd_now  = qd_all{p};
    qdd_now = qdd_all{p};
    u_now   = u_all{p};

    jerk_now = diff(qdd_now,1,2)/dt;
    max_abs_acc_all(p,:)   = max(abs(qdd_now),[],2).';
    mean_abs_jerk_all(p,:) = mean(abs(jerk_now),2).';
    rms_delta_tau_all(p,:) = sqrt(mean(u_now.^2,2)).';
    peak_delta_tau_all(p,:) = max(abs(u_now),[],2).';
    total_variation_tau_all(p) = sum(sum(abs(diff(u_now,1,2))));

    final_q_error_all(p)  = norm(q_now(:,end)-qg,2);
    final_qd_error_all(p) = norm(qd_now(:,end)-qdg,2);
    max_overshoot_all(p)  = computeMaxJointOvershoot(q_now, q0, qg);
    settling_time_all(p)  = computeSettlingTime(t, q_now, qd_now, qg, qdg, 0.02, 0.05);

    if use_delta_tau_bounds
        lbViol = max(delta_tau_min - u_now, [], 'all');
        ubViol = max(u_now - delta_tau_max, [], 'all');
        constraint_violation_all(p) = max([0, lbViol, ubViol]);
    end

    endpoint_terminal_error_all(p) = norm(paths_all{p}(:,end) - goal_endpoint, 2);
    straightDistance = norm(goal_endpoint - start_endpoint, 2);
    endpoint_straightness_all(p) = path_length_all(p) / max(straightDistance, eps);
end

summaryTable = table(plannerNames(:), total_cost_all, final_energy_actual_all, peak_power_actual_all, ...
    final_q_error_all, final_qd_error_all, endpoint_terminal_error_all, path_length_all, endpoint_straightness_all, ...
    max_overshoot_all, settling_time_all, constraint_violation_all, total_variation_tau_all, ...
    'VariableNames', {'Planner','TotalCost','ActualMechanicalEnergy_J','PeakActualPower_W', ...
    'TerminalPositionError_rad','TerminalVelocityError_radps','TerminalEndpointError_m','EndpointPathLength_m', ...
    'EndpointStraightnessRatio','MaxJointOvershoot_rad','SettlingTime_s','MaxInputBoundViolation_Nm','TotalDeltaTauVariation_Nm'});

jointMetricTable = table(plannerNames(:), ...
    max_abs_acc_all(:,1), max_abs_acc_all(:,2), max_abs_acc_all(:,3), ...
    mean_abs_jerk_all(:,1), mean_abs_jerk_all(:,2), mean_abs_jerk_all(:,3), ...
    rms_delta_tau_all(:,1), rms_delta_tau_all(:,2), rms_delta_tau_all(:,3), ...
    peak_delta_tau_all(:,1), peak_delta_tau_all(:,2), peak_delta_tau_all(:,3), ...
    'VariableNames', {'Planner','PeakAbsQdd_J1','PeakAbsQdd_J2','PeakAbsQdd_J3', ...
    'MeanAbsJerk_J1','MeanAbsJerk_J2','MeanAbsJerk_J3', ...
    'RMSDeltaTau_J1','RMSDeltaTau_J2','RMSDeltaTau_J3', ...
    'PeakAbsDeltaTau_J1','PeakAbsDeltaTau_J2','PeakAbsDeltaTau_J3'});

writetable(summaryTable, fullfile(tableDir, 'mpc_pid_lqr_comparison_summary.csv'));
writetable(jointMetricTable, fullfile(tableDir, 'mpc_pid_lqr_joint_smoothness_control_metrics.csv'));
writeTextFile(fullfile(tableDir, 'mpc_pid_lqr_comparison_summary.tex'), tableToLatex(summaryTable, 'MPC--PID--LQR trajectory comparison summary'));
writeTextFile(fullfile(tableDir, 'mpc_pid_lqr_joint_smoothness_control_metrics.tex'), tableToLatex(jointMetricTable, 'Joint smoothness and control-effort metrics'));
save(fullfile(dataDir, 'planner_data_publication.mat'), 'planner_data', 'summaryTable', 'jointMetricTable');
save('planner_data.mat', 'planner_data');

fig = pubFigure('Joint Position Comparison');
plotJointFamily(t, q_all, labelsLatex, styles, '$q_j(t)$ [rad]', ...
    'Reduced UR5 joint position trajectories', jointLabels, pub, false);
savePublicationFigure(fig, figDir, 'comparison_joint_positions_overlay', pub);

fig = pubFigure('Joint Velocity Comparison');
plotJointFamily(t, qd_all, labelsLatex, styles, '$\dot{q}_j(t)$ [rad/s]', ...
    'Reduced UR5 joint velocity trajectories', jointLabels, pub, false);
savePublicationFigure(fig, figDir, 'comparison_joint_velocities_overlay', pub);

fig = pubFigure('Joint Acceleration Comparison');
plotJointFamily(t, qdd_all, labelsLatex, styles, '$\ddot{q}_j(t)$ [rad/s$^2$]', ...
    'Reduced UR5 joint acceleration trajectories', jointLabels, pub, false);
savePublicationFigure(fig, figDir, 'comparison_joint_accelerations_overlay', pub);

fig = pubFigure('3D Endpoint Path Comparison');
hold on; grid on; axis equal; box on;
for p = 1:numPlanners
    plot3(paths_all{p}(1,:), paths_all{p}(2,:), paths_all{p}(3,:), styles{p}, 'LineWidth', pub.lineWidthMain);
end
scatter3(start_endpoint(1), start_endpoint(2), start_endpoint(3), pub.markerSize, 'filled', 'MarkerEdgeColor','k');
scatter3(goal_endpoint(1),  goal_endpoint(2),  goal_endpoint(3),  pub.markerSize, 'filled', 'MarkerEdgeColor','k');
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]');
title(sprintf('Reduced UR5 end-effector paths (%s)', reducedBody));
legend([labelsLatex, {'Start','Goal'}], 'Location','bestoutside');
pubAxes(gca, pub); view(35,25);
savePublicationFigure(fig, figDir, 'comparison_end_effector_path_3d', pub);

fig = pubFigure('Endpoint Path Projection Comparison');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
projPairs = [1 2; 1 3; 2 3];
projNames = {'$x$--$y$ projection','$x$--$z$ projection','$y$--$z$ projection'};
axisNames = {'$x$ [m]','$y$ [m]','$z$ [m]'};
for ii = 1:3
    nexttile; hold on; grid on; box on; axis equal;
    a = projPairs(ii,1); b = projPairs(ii,2);
    for p = 1:numPlanners
        plot(paths_all{p}(a,:), paths_all{p}(b,:), styles{p}, 'LineWidth', pub.lineWidthMain);
    end
    scatter(start_endpoint(a), start_endpoint(b), pub.markerSize, 'filled', 'MarkerEdgeColor','k');
    scatter(goal_endpoint(a), goal_endpoint(b), pub.markerSize, 'filled', 'MarkerEdgeColor','k');
    xlabel(axisNames{a}); ylabel(axisNames{b}); title(projNames{ii}); pubAxes(gca,pub);
end
legend([labelsLatex, {'Start','Goal'}], 'Location','bestoutside');
savePublicationFigure(fig, figDir, 'comparison_end_effector_path_projections', pub);

for p = 1:numPlanners
    cleanName = matlab.lang.makeValidName(plannerNames{p});

    fig = pubFigure(sprintf('%s Joint Position Profiles', plannerNames{p}));
    plotSinglePlannerJointProfiles(t, q_all{p}, '$q_j(t)$ [rad]', ...
        sprintf('%s joint position profiles', labelsLatex{p}), jointLabels, pub, false);
    savePublicationFigure(fig, figDir, sprintf('%s_joint_position_profiles', cleanName), pub);

    fig = pubFigure(sprintf('%s Joint Velocity Profiles', plannerNames{p}));
    plotSinglePlannerJointProfiles(t, qd_all{p}, '$\dot{q}_j(t)$ [rad/s]', ...
        sprintf('%s joint velocity profiles', labelsLatex{p}), jointLabels, pub, false);
    savePublicationFigure(fig, figDir, sprintf('%s_joint_velocity_profiles', cleanName), pub);

    fig = pubFigure(sprintf('%s Joint Acceleration Profiles', plannerNames{p}));
    plotSinglePlannerJointProfiles(t, qdd_all{p}, '$\ddot{q}_j(t)$ [rad/s$^2$]', ...
        sprintf('%s joint acceleration profiles', labelsLatex{p}), jointLabels, pub, false);
    savePublicationFigure(fig, figDir, sprintf('%s_joint_acceleration_profiles', cleanName), pub);

    fig = pubFigure(sprintf('%s Delta Torque Profiles', plannerNames{p}));
    plotSinglePlannerJointProfiles(t(1:N), u_all{p}, '$\Delta\tau_j$ [N m]', ...
        sprintf('%s gravity-compensated torque input', labelsLatex{p}), jointLabels, pub, true);
    savePublicationFigure(fig, figDir, sprintf('%s_delta_torque_profiles', cleanName), pub);

    fig = pubFigure(sprintf('%s Planned Torque Decomposition', plannerNames{p}));
    plotTorqueDecomposition4x1(t(1:N), tau_components_all{p}, pub, sprintf('%s planned torque decomposition', labelsLatex{p}));
    savePublicationFigure(fig, figDir, sprintf('%s_planned_torque_decomposition_4x1', cleanName), pub);
end

fig = pubFigure('Overlayed Torque Inputs');
plotJointFamily(t(1:N), u_all, labelsLatex, styles, '$\Delta\tau_j(t)$ [N m]', ...
    'Overlayed gravity-compensated torque inputs', jointLabels, pub, true);
savePublicationFigure(fig, figDir, 'comparison_overlayed_delta_torque_inputs', pub);

fig = pubFigure('Overlayed Absolute Torque Inputs');
plotJointFamily(t(1:N), {tau_components_all{1}.absolute, tau_components_all{2}.absolute, tau_components_all{3}.absolute}, ...
    labelsLatex, styles, '$\tau_j(t)$ [N m]', 'Overlayed absolute inverse-dynamics torques', jointLabels, pub, false);
savePublicationFigure(fig, figDir, 'comparison_overlayed_absolute_torque_inputs', pub);

fig = pubFigure('State Convergence Comparison');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t, vecnorm(q_all{p}-qg,2,1), styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('$\|q-q_g\|_2$ [rad]'); title('Position-state convergence'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t, vecnorm(qd_all{p}-qdg,2,1), styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('$\|\dot{q}-\dot{q}_g\|_2$ [rad/s]'); xlabel('Time [s]'); title('Velocity-state convergence'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_state_convergence', pub);

fig = pubFigure('Stage Cost Comparison');
hold on; grid on; box on;
for p = 1:numPlanners
    stairs(t(1:N), stage_cost_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
xlabel('Time [s]'); ylabel('Stage cost'); title('Stage cost comparison'); legend(labelsLatex,'Location','best'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_stage_cost', pub);

fig = pubFigure('Cumulative Cost Comparison');
hold on; grid on; box on;
for p = 1:numPlanners
    stairs(t(1:N), cumulative_cost_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
xlabel('Time [s]'); ylabel('Cumulative cost'); title('Cumulative cost comparison'); legend(labelsLatex,'Location','best'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_cumulative_cost', pub);

fig = pubFigure('Cost-to-Go Comparison');
hold on; grid on; box on;
for p = 1:numPlanners
    stairs(t(1:N), cost_to_go_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
xlabel('Time [s]'); ylabel('Cost-to-go'); title('Cost-to-go comparison'); legend(labelsLatex,'Location','best'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_cost_to_go', pub);

fig = pubFigure('RMS Control Effort Comparison');
bar(rms_delta_tau_all); grid on; box on;
setPlannerTicks(gca, labelsLatex);
xlabel('Planner/controller'); ylabel('RMS $\Delta\tau_j$ [N m]');
title('RMS control effort per joint'); legend({'Joint 1','Joint 2','Joint 3'},'Location','bestoutside'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_rms_control_effort', pub);

fig = pubFigure('Planned Power Comparison');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), total_power_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('Signed power [W]'); title('Signed commanded power'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), abs_total_power_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('Absolute power [W]'); xlabel('Time [s]'); title('Absolute commanded power'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_commanded_power', pub);

fig = pubFigure('Planned Energy Comparison');
hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), cum_energy_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
xlabel('Time [s]'); ylabel('Cumulative commanded energy [J]'); title('Commanded energy comparison'); legend(labelsLatex,'Location','best'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_commanded_energy', pub);

fig = pubFigure('Actual Mechanical Power Comparison');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), total_power_actual_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('Signed power [W]'); title('Signed actual mechanical power'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
nexttile; hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), abs_total_power_actual_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
ylabel('Absolute power [W]'); xlabel('Time [s]'); title('Absolute actual mechanical power'); legend(labelsLatex,'Location','bestoutside'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_actual_mechanical_power', pub);

fig = pubFigure('Actual Mechanical Energy Comparison');
hold on; grid on; box on;
for p = 1:numPlanners
    plot(t(1:N), cum_energy_actual_all{p}, styles{p}, 'LineWidth', pub.lineWidthMain);
end
xlabel('Time [s]'); ylabel('Cumulative actual mechanical energy [J]'); title('Actual mechanical energy comparison'); legend(labelsLatex,'Location','best'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_actual_mechanical_energy', pub);

fig = pubFigure('Energy and Peak Power Summary');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; bar(final_energy_actual_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Actual energy [J]'); title('Final actual mechanical energy'); pubAxes(gca,pub);
nexttile; bar(peak_power_actual_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Peak absolute power [W]'); title('Peak actual mechanical power'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_energy_peak_power_summary', pub);

fig = pubFigure('Endpoint Error and Path Length Summary');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; bar(endpoint_terminal_error_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Terminal endpoint error [m]'); title('End-effector terminal accuracy'); pubAxes(gca,pub);
nexttile; bar(path_length_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Path length [m]'); title('End-effector path length'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_endpoint_error_path_length', pub);

fig = pubFigure('Trajectory Quality Summary');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
nexttile; bar(sum(mean_abs_jerk_all,2)); grid on; box on; title('Smoothness'); ylabel('$\sum_j$ mean $|\dddot q_j|$'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
nexttile; bar(max(peak_delta_tau_all,[],2)); grid on; box on; title('Aggressiveness'); ylabel('$\max_j |\Delta\tau_j|$ [N m]'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
nexttile; bar(final_q_error_all); grid on; box on; title('Terminal accuracy'); ylabel('$\|q_N-q_g\|_2$ [rad]'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
nexttile; bar(total_variation_tau_all); grid on; box on; title('Control variation'); ylabel('$\sum |\Delta\tau_k-\Delta\tau_{k-1}|$'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
nexttile; bar(constraint_violation_all); grid on; box on; title('Constraint handling'); ylabel('Max violation [N m]'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
nexttile; bar(endpoint_straightness_all); grid on; box on; title('Motion quality'); ylabel('Path / straight-line length'); setPlannerTicks(gca,labelsLatex); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_quality_summary_2x3', pub);

fig = pubFigure('Settling and Overshoot Summary');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; bar(settling_time_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Settling time [s]'); title('Settling behavior'); pubAxes(gca,pub);
nexttile; bar(max_overshoot_all); grid on; box on; setPlannerTicks(gca,labelsLatex);
ylabel('Max overshoot [rad]'); title('Joint-space overshoot'); pubAxes(gca,pub);
savePublicationFigure(fig, figDir, 'comparison_settling_overshoot_summary', pub);

%% Console Summary
fprintf('\nReduced 3-DoF UR5 planner generation completed successfully.\n');
fprintf('Reduced endpoint body               : %s\n', reducedBody);
fprintf('Horizon intervals N                 : %d\n', N);
fprintf('Sampling time dt                    : %.4f s\n', dt);
fprintf('Total motion time Tf                : %.4f s\n', Tf);
fprintf('Optimal planner type                : online receding-horizon LTV MPC\n');
fprintf('Prediction horizon max Np           : %d\n', Np_max);
fprintf('quadprog available/used             : %d\n', quadprog_available);
fprintf('MPC final nonlinear error norm      : %.6e\n', terminal_error_norm_mpc);
fprintf('PID final nonlinear error norm      : %.6e\n', pid_info.final_error_norm);
fprintf('LQR final nonlinear error norm      : %.6e\n', lqr_info.final_error_norm);
fprintf('MPC planned cost J                  : %.8g\n', J_plan_opt);
fprintf('Output saved to                     : planner_data.mat\n');
fprintf('Publication figures saved to        : %s\n', figDir);
fprintf('Publication tables saved to         : %s\n', tableDir);
fprintf('Publication data saved to           : %s\n\n', dataDir);

fprintf('------------------------------------------------------------\n');
fprintf('Planner/controller comparison summary (with actual mechanical power):\n');
for p = 1:numPlanners
    fprintf(['%-15s | Total cost: %12.6g | Actual energy: %12.6g J | Actual peak power: %12.6g W' ...
             ' | Path length: %12.6g m | Terminal endpoint error: %12.6g m | Terminal q error: %12.6g rad\n'], ...
        plannerNames{p}, total_cost_all(p), final_energy_actual_all(p), peak_power_actual_all(p), ...
        path_length_all(p), endpoint_terminal_error_all(p), final_q_error_all(p));
end
fprintf('------------------------------------------------------------\n');

% % LOCAL HELPER FUNCTIONS

function ensureDir(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end

function fig = pubFigure(figName)
fig = figure('Name', figName, 'Color', 'w', 'Units', 'pixels', 'Position', [80 80 1220 840]);
end

function pubAxes(ax, pub)
set(ax, 'FontSize', pub.fontSizeAxes, 'LineWidth', 1.0, 'Box', 'on', ...
    'TickLabelInterpreter', 'latex', 'XMinorGrid', 'on', 'YMinorGrid', 'on');
grid(ax, 'on');
end

function savePublicationFigure(fig, figDir, baseName, pub)
for iFmt = 1:numel(pub.formats)
    fmt = pub.formats{iFmt};
    outPath = fullfile(figDir, [baseName '.' fmt]);
    switch lower(fmt)
        case 'fig'
            savefig(fig, outPath);
        case 'pdf'
            exportgraphics(fig, outPath, 'ContentType', 'vector');
        case 'png'
            exportgraphics(fig, outPath, 'Resolution', pub.dpi);
        otherwise
            saveas(fig, outPath);
    end
end
end

function plotJointFamily(tPlot, dataCells, labelsLatex, styles, yLabelText, plotTitle, jointLabels, pub, useStairs)
if nargin < 9
    useStairs = false;
end
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
numCurves = numel(dataCells);
for j = 1:3
    nexttile; hold on; grid on; box on;
    for p = 1:numCurves
        y = dataCells{p}(j,:).';
        if useStairs
            stairs(tPlot(:), y, styles{p}, 'LineWidth', pub.lineWidthMain);
        else
            plot(tPlot(:), y, styles{p}, 'LineWidth', pub.lineWidthMain);
        end
    end
    ylabel(yLabelText);
    title(sprintf('%s, %s', plotTitle, jointLabels{j}));
    if j == 3, xlabel('Time [s]'); end
    pubAxes(gca, pub);
end
legend(labelsLatex, 'Location', 'bestoutside');
end

function plotSinglePlannerJointProfiles(tPlot, dataMat, yLabelText, plotTitle, jointLabels, pub, useStairs)
if nargin < 7
    useStairs = false;
end
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for j = 1:3
    nexttile; hold on; grid on; box on;
    if useStairs
        stairs(tPlot(:), dataMat(j,:).', 'LineWidth', pub.lineWidthMain);
    else
        plot(tPlot(:), dataMat(j,:).', 'LineWidth', pub.lineWidthMain);
    end
    ylabel(yLabelText);
    title(sprintf('%s, %s', plotTitle, jointLabels{j}));
    if j == 3, xlabel('Time [s]'); end
    pubAxes(gca, pub);
end
end

function plotTorqueDecomposition4x1(tPlot, tauComp, pub, figTitle)
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
fields = {'absolute','inertial','coriolis','gravity'};
titles = {'Absolute torque $\tau$', 'Inertial torque $M(q)\ddot{q}$', ...
    'Coriolis/centrifugal torque $C(q,\dot{q})$', 'Gravity compensation $G(q)$'};
for ii = 1:4
    nexttile; hold on; grid on; box on;
    dataNow = tauComp.(fields{ii});
    for j = 1:3
        plot(tPlot(:), dataNow(j,:).', 'LineWidth', pub.lineWidthMain);
    end
    ylabel('[N m]'); title(sprintf('%s: %s', figTitle, titles{ii}));
    if ii == 4, xlabel('Time [s]'); end
    pubAxes(gca,pub);
end
legend({'Joint 1','Joint 2','Joint 3'}, 'Location','bestoutside');
end

function setPlannerTicks(ax, labelsLatex)
set(ax, 'XTick', 1:numel(labelsLatex), 'XTickLabel', labelsLatex);
xtickangle(ax, 20);
end

function writeTextFile(filePath, txt)
fid = fopen(filePath, 'w');
assert(fid > 0, 'Could not open file for writing: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
end

function tex = tableToLatex(T, captionText)
varNames = T.Properties.VariableNames;
tex = sprintf('%% Auto-generated by MATLAB publication export section\n');
tex = [tex, sprintf('\\begin{table}[!t]\n\\centering\n')];
tex = [tex, sprintf('\\caption{%s}\n', captionText)];
tex = [tex, sprintf('\\begin{tabular}{%s}\n\\hline\n', repmat('c',1,numel(varNames)))];
tex = [tex, strjoin(strrep(varNames,'_','\_'), ' & '), sprintf(' \\\\ \n\\hline\n')];
for r = 1:height(T)
    rowCells = cell(1,numel(varNames));
    for c = 1:numel(varNames)
        value = T{r,c};
        if iscell(value), value = value{1}; end
        if isstring(value), value = char(value); end
        if isnumeric(value)
            rowCells{c} = sprintf('%.4g', value);
        else
            rowCells{c} = char(value);
        end
        rowCells{c} = strrep(rowCells{c}, '_', '\_');
    end
    tex = [tex, strjoin(rowCells, ' & '), sprintf(' \\\\ \n')];
end
tex = [tex, sprintf('\\hline\n\\end{tabular}\n\\end{table}\n')];
end

function tauComp = computeTorqueComponentsReduced(robot, idx, q_fix_456, q, qd, qdd, N)
tauComp.absolute = zeros(3,N);
tauComp.inertial = zeros(3,N);
tauComp.coriolis = zeros(3,N);
tauComp.gravity  = zeros(3,N);
for k = 1:N
    q_full   = [q(:,k);   q_fix_456(:)];
    qd_full  = [qd(:,k);  zeros(3,1)];
    qdd_full = [qdd(:,k); zeros(3,1)];
    M6 = massMatrix(robot, q_full);
    G6 = gravityTorque(robot, q_full);
    if exist('velocityProduct','file') == 2 || exist('velocityProduct','builtin') == 5
        C6 = velocityProduct(robot, q_full, qd_full);
    else
        gBackup = robot.Gravity;
        robot.Gravity = [0 0 0];
        C6 = inverseDynamics(robot, q_full, qd_full, zeros(6,1));
        robot.Gravity = gBackup;
    end
    tauAbs6 = inverseDynamics(robot, q_full, qd_full, qdd_full);
    tauComp.absolute(:,k) = tauAbs6(idx);
    tauComp.inertial(:,k) = M6(idx,:)*qdd_full;
    tauComp.coriolis(:,k) = C6(idx);
    tauComp.gravity(:,k)  = G6(idx);
end
end

function maxOvershoot = computeMaxJointOvershoot(q, q0, qg)
maxOvershoot = 0;
for j = 1:3
    lo = min(q0(j), qg(j));
    hi = max(q0(j), qg(j));
    over = max([0, max(q(j,:) - hi), max(lo - q(j,:))]);
    maxOvershoot = max(maxOvershoot, over);
end
end

function ts = computeSettlingTime(t, q, qd, qg, qdg, posTolFrac, velTolAbs)
travel = max(abs(qg(:) - q(:,1)), 1e-6);
posTol = posTolFrac*travel;
posErr = abs(q - qg(:));
velErr = abs(qd - qdg(:));
inside = all(posErr <= posTol, 1) & all(velErr <= velTolAbs, 1);
ts = NaN;
for k = 1:numel(t)
    if all(inside(k:end))
        ts = t(k);
        return;
    end
end
end

function [q, qd, u, info] = generatePIDGoalRegulationTrajectory(f_error, robot, idx, q_fix_456, q0, qd0, qg, qdg, t, Kp, Ki, Kd, integralLimit, useBounds, uMin, uMax)
% Closed-loop nonlinear PID goal-regulation trajectory using delta_tau input.
dt = t(2) - t(1);
N = numel(t) - 1;
x = [q0(:)-qg(:); qd0(:)-qdg(:)];
q = zeros(3,N+1);
qd = zeros(3,N+1);
u = zeros(3,N);
ei = zeros(3,1);
for k = 1:N
    q(:,k) = x(1:3) + qg;
    qd(:,k) = x(4:6) + qdg;

    eq = qg - q(:,k);
    eqd = qdg - qd(:,k);
    ei = ei + eq*dt;
    ei = min(max(ei, -integralLimit(:)), integralLimit(:));

    delta_tau = Kp*eq + Ki*ei + Kd*eqd;
    if useBounds
        delta_tau = saturateVector(delta_tau, uMin, uMax);
    end
    u(:,k) = delta_tau;
    x = rk4Step(f_error, x, delta_tau, dt);
end
q(:,N+1) = x(1:3) + qg;
qd(:,N+1) = x(4:6) + qdg;

info.final_error = [q(:,end)-qg; qd(:,end)-qdg];
info.final_error_norm = norm(info.final_error);
info.integral_state_final = ei;
info.input_convention = 'tau_abs = G(q) + delta_tau; nonlinear propagation uses delta_tau';
info.note = 'robot, idx, and q_fix_456 are included for interface symmetry and downstream logging compatibility.';
end

function [K, info] = designGoalLQRGain(f_error, dt, eps_x, eps_u, Qlqr, Rlqr)
% Linearize and discretize the reduced nonlinear error dynamics at the goal.
x_eq = zeros(6,1);
u_eq = zeros(3,1);
[Ac, Bc] = numericalJacobians(f_error, x_eq, u_eq, eps_x, eps_u);
[Ad, Bd] = zohDiscretize(Ac, Bc, dt);
if exist('dlqr','file') == 2 || exist('dlqr','builtin') == 5
    [K,~,eigCL] = dlqr(Ad, Bd, Qlqr, Rlqr);
else
    warning('dlqr not found. Falling back to finite-iteration discrete Riccati recursion. Control System Toolbox is recommended.');
    P = Qlqr;
    for iter = 1:1000
        Pnext = Ad.'*P*Ad - Ad.'*P*Bd*((Rlqr + Bd.'*P*Bd)\(Bd.'*P*Ad)) + Qlqr;
        if norm(Pnext-P,'fro') < 1e-10
            P = Pnext;
            break;
        end
        P = Pnext;
    end
    K = (Rlqr + Bd.'*P*Bd)\(Bd.'*P*Ad);
    eigCL = eig(Ad - Bd*K);
end
info.Ac = Ac;
info.Bc = Bc;
info.Ad = Ad;
info.Bd = Bd;
info.closed_loop_eigs = eigCL;
info.linearization_point = 'x = zeros(6,1), delta_tau = zeros(3,1) at q = qg, qd = qdg';
end

function [q, qd, u, info] = generateLQRGoalRegulationTrajectory(f_error, robot, idx, q_fix_456, q0, qd0, qg, qdg, t, K, useBounds, uMin, uMax)
% Closed-loop nonlinear LQR goal-regulation trajectory using delta_tau = -K*x.
dt = t(2) - t(1);
N = numel(t) - 1;
x = [q0(:)-qg(:); qd0(:)-qdg(:)];
q = zeros(3,N+1);
qd = zeros(3,N+1);
u = zeros(3,N);
for k = 1:N
    q(:,k) = x(1:3) + qg;
    qd(:,k) = x(4:6) + qdg;

    delta_tau = -K*x;
    if useBounds
        delta_tau = saturateVector(delta_tau, uMin, uMax);
    end
    u(:,k) = delta_tau;
    x = rk4Step(f_error, x, delta_tau, dt);
end
q(:,N+1) = x(1:3) + qg;
qd(:,N+1) = x(4:6) + qdg;

info.final_error = [q(:,end)-qg; qd(:,end)-qdg];
info.final_error_norm = norm(info.final_error);
info.input_convention = 'tau_abs = G(q) + delta_tau; nonlinear propagation uses delta_tau';
info.note = 'robot, idx, and q_fix_456 are included for interface symmetry and downstream logging compatibility.';
end

function y = saturateVector(x, xmin, xmax)
y = min(max(x(:), xmin(:)), xmax(:));
end

function bodyName = pickReducedUR5EndpointBody(robot)
% Select a body/frame whose pose is determined by joints 1-3 only.
candidates = {'wrist_1_link','forearm_link'};
bodyName = '';
for i = 1:numel(candidates)
    if any(strcmp(robot.BodyNames, candidates{i}))
        bodyName = candidates{i};
        return;
    end
end
error('No suitable reduced UR5 endpoint body found. Expected one of: %s', ...
    strjoin(candidates, ', '));
end

function p = reducedFKPosition(robot, bodyName, q123, q_fix_456)
% Forward kinematics of the reduced 3-DoF UR5 endpoint position.
q_full = [q123(:); q_fix_456(:)];
T = getTransform(robot, q_full, bodyName);
p = T(1:3,4);
end

function path = reducedPathFromTrajectory(robot, bodyName, q_traj, q_fix_456)
% Convert a 3xM joint trajectory into a 3xM endpoint position path.
numSamples = size(q_traj,2);
path = zeros(3, numSamples);
for k = 1:numSamples
    path(:,k) = reducedFKPosition(robot, bodyName, q_traj(:,k), q_fix_456);
end
end

function xdot = reducedErrorDynamics(robot, idx, q_fix_456, qg, qdg, x, delta_tau)
% Reduced nonlinear error dynamics for active UR5 joints.
% Full dynamics:
% D(q)qdd + C(q,qd) + G(q) = tau_abs
% Gravity-compensated input convention:
% tau_abs = G(q) + delta_tau
% Therefore:
% D(q)qdd = delta_tau - C(q,qd)
% State:
% x = [q-qg; qd-qdg]

dq  = x(1:3);
dqd = x(4:6);

q  = qg  + dq;
qd = qdg + dqd;

q_full  = [q;  q_fix_456(:)];
qd_full = [qd; zeros(3,1)];

D3 = reducedMassMatrix(robot, q_full, idx);
C3 = reducedVelocityProduct(robot, q_full, qd_full, idx);

qdd = D3 \ (delta_tau(:) - C3);
xdot = [dqd; qdd];
end

function D3 = reducedMassMatrix(robot, q_full, idx)
D6 = massMatrix(robot, q_full);
D3 = D6(idx, idx);
end

function C3 = reducedVelocityProduct(robot, q_full, qd_full, idx)
% Reduced velocity-dependent generalized torque term.
if exist('velocityProduct','file') == 2 || exist('velocityProduct','builtin') == 5
    C6 = velocityProduct(robot, q_full, qd_full);
    C3 = C6(idx);
    return;
end

gravityBackup = robot.Gravity;
cleanupObj = onCleanup(@() restoreRobotGravity(robot, gravityBackup));
robot.Gravity = [0 0 0];
tau_vel_6 = inverseDynamics(robot, q_full, qd_full, zeros(size(qd_full)));
C3 = tau_vel_6(idx);
end

function restoreRobotGravity(robot, gravityVector)
robot.Gravity = gravityVector;
end

function [Ac, Bc] = numericalJacobians(f, x0, u0, eps_x, eps_u)
% Central-difference numerical Jacobians of continuous dynamics.
n = numel(x0);
m = numel(u0);
Ac = zeros(n,n);
Bc = zeros(n,m);
for i = 1:n
    dx = zeros(n,1);
    dx(i) = eps_x;
    Ac(:,i) = (f(x0 + dx, u0) - f(x0 - dx, u0)) / (2*eps_x);
end
for j = 1:m
    du = zeros(m,1);
    du(j) = eps_u;
    Bc(:,j) = (f(x0, u0 + du) - f(x0, u0 - du)) / (2*eps_u);
end
end

function [Ad, Bd] = zohDiscretize(Ac, Bc, dt)
% Zero-order-hold discretization using matrix exponential.
n = size(Ac,1);
m = size(Bc,2);
M = [Ac, Bc;
     zeros(m,n+m)];
Md = expm(M*dt);
Ad = Md(1:n, 1:n);
Bd = Md(1:n, n+1:n+m);
end

% Advance the nonlinear dynamics by one fourth-order Runge-Kutta step.
function xnext = rk4Step(f, x, u, dt)
% One Runge-Kutta 4 integration step.
k1 = f(x, u);
k2 = f(x + 0.5*dt*k1, u);
k3 = f(x + 0.5*dt*k2, u);
k4 = f(x + dt*k3, u);
xnext = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function [x_nom, u_nom] = rolloutNominalTrajectory(f, x0, U, dt, Np)
% Roll out a nonlinear nominal trajectory using the current candidate U.
m = numel(U)/Np;
x_nom = zeros(numel(x0), Np+1);
u_nom = reshape(U, m, Np);
x_nom(:,1) = x0;
for i = 1:Np
    x_nom(:,i+1) = rk4Step(f, x_nom(:,i), u_nom(:,i), dt);
end
end

function U_warm_next = shiftWarmStart(U_current, Np_current, Np_max, m)
% Shift optimal predicted controls by one sample and pad with the final value.
Umat = reshape(U_current, m, Np_current);
if Np_current > 1
    lastCol = Umat(:,end);
    Ushift = [Umat(:,2:end), lastCol];
else
    Ushift = Umat;
end
if size(Ushift,2) < Np_max
    Ushift = [Ushift, repmat(Ushift(:,end), 1, Np_max-size(Ushift,2))];
elseif size(Ushift,2) > Np_max
    Ushift = Ushift(:,1:Np_max);
end
U_warm_next = Ushift(:);
end

function [Abar, Bbar, cbar] = buildStackedPredictionMatricesLTVAffine(Aseq, Bseq, cseq, Np)
% Build stacked affine LTV prediction matrices:
% x_{i+1} = A_i x_i + B_i u_i + c_i
% X = Abar*x0 + Bbar*U + cbar
% where X = [x_1; x_2; ...; x_Np].
n = size(Aseq{1},1);
m = size(Bseq{1},2);
Abar = zeros(n*Np, n);
Bbar = zeros(n*Np, m*Np);
cbar = zeros(n*Np, 1);

Aprev = eye(n);
Bprev = zeros(n, m*Np);
cprev = zeros(n,1);
for i = 1:Np
    Ai = Aseq{i};
    Bi = Bseq{i};
    ci = cseq{i};

    Anew = Ai*Aprev;
    Bnew = Ai*Bprev;
    Bnew(:, (i-1)*m+1:i*m) = Bi;
    cnew = Ai*cprev + ci;

    Abar((i-1)*n+1:i*n, :) = Anew;
    Bbar((i-1)*n+1:i*n, :) = Bnew;
    cbar((i-1)*n+1:i*n) = cnew;

    Aprev = Anew;
    Bprev = Bnew;
    cprev = cnew;
end
end

function qdd = estimateAccelerationFromVelocity(qd, dt)
% Numerically estimate acceleration from velocity samples.
qdd = zeros(size(qd));
if size(qd,2) == 1
    return;
end
qdd(:,1) = (qd(:,2) - qd(:,1)) / dt;
for k = 2:size(qd,2)-1
    qdd(:,k) = (qd(:,k+1) - qd(:,k-1)) / (2*dt);
end
qdd(:,end) = (qd(:,end) - qd(:,end-1)) / dt;
end

function [q, qd, qdd] = generateCubicTrajectory(q0, qd0, qf, qdf, t)
% Joint-wise cubic polynomial trajectory with endpoint position/velocity.
t = t(:);
T = t(end) - t(1);
tau = t - t(1);
numJoints = numel(q0);
numSamples = numel(t);
q = zeros(numJoints, numSamples);
qd = zeros(numJoints, numSamples);
qdd = zeros(numJoints, numSamples);
for j = 1:numJoints
    a0 = q0(j);
    a1 = qd0(j);
    M = [T^2,   T^3;
         2*T, 3*T^2];
    b = [qf(j) - a0 - a1*T;
         qdf(j) - a1];
    a23 = M\b;
    a2 = a23(1);
    a3 = a23(2);
    q(j,:)   = a0 + a1*tau + a2*tau.^2 + a3*tau.^3;
    qd(j,:)  = a1 + 2*a2*tau + 3*a3*tau.^2;
    qdd(j,:) = 2*a2 + 6*a3*tau;
end
end

function [q, qd, qdd] = generateQuinticTrajectory(q0, qd0, qdd0, qf, qdf, qddf, t)
% Joint-wise quintic polynomial trajectory with pos/vel/acc endpoints.
t = t(:);
T = t(end) - t(1);
tau = t - t(1);
numJoints = numel(q0);
numSamples = numel(t);
q = zeros(numJoints, numSamples);
qd = zeros(numJoints, numSamples);
qdd = zeros(numJoints, numSamples);
for j = 1:numJoints
    a0 = q0(j);
    a1 = qd0(j);
    a2 = 0.5*qdd0(j);
    M = [T^3,    T^4,     T^5;
         3*T^2,  4*T^3,   5*T^4;
         6*T,   12*T^2,  20*T^3];
    b = [qf(j)   - (a0 + a1*T + a2*T^2);
         qdf(j)  - (a1 + 2*a2*T);
         qddf(j) - (2*a2)];
    a345 = M\b;
    a3 = a345(1);
    a4 = a345(2);
    a5 = a345(3);
    q(j,:) = a0 + a1*tau + a2*tau.^2 + a3*tau.^3 + a4*tau.^4 + a5*tau.^5;
    qd(j,:) = a1 + 2*a2*tau + 3*a3*tau.^2 + 4*a4*tau.^3 + 5*a5*tau.^4;
    qdd(j,:) = 2*a2 + 6*a3*tau + 12*a4*tau.^2 + 20*a5*tau.^3;
end
end

function [q, qd, qdd] = generateSynchronizedTrapezoidalTrajectory(q0, qf, t, accelFraction)
% Joint-wise synchronized trapezoidal velocity trajectory.
% All joints share the same normalized scalar progress s(t).
t = t(:);
T = t(end) - t(1);
tau = t - t(1);
numSamples = numel(t);
numJoints = numel(q0);

accelFraction = max(min(accelFraction, 0.49), 1e-3);
ta = accelFraction*T;
tc = T - 2*ta;

% Unit-distance trapezoid parameters.
vmax = 1/(T - ta);
acc = vmax/ta;

s = zeros(1,numSamples);
sd = zeros(1,numSamples);
sdd = zeros(1,numSamples);
for k = 1:numSamples
    tk = tau(k);
    if tk <= ta
        s(k) = 0.5*acc*tk^2;
        sd(k) = acc*tk;
        sdd(k) = acc;
    elseif tk <= ta + tc
        s(k) = 0.5*acc*ta^2 + vmax*(tk - ta);
        sd(k) = vmax;
        sdd(k) = 0;
    else
        td = tk - ta - tc;
        s(k) = 0.5*acc*ta^2 + vmax*tc + vmax*td - 0.5*acc*td^2;
        sd(k) = vmax - acc*td;
        sdd(k) = -acc;
    end
end

% Ensure exact endpoint values.
s(1) = 0; s(end) = 1;
sd(1) = 0; sd(end) = 0;

Delta = qf(:) - q0(:);
q = q0(:) + Delta*s;
qd = Delta*sd;
qdd = Delta*sdd;
assert(isequal(size(q), [numJoints, numSamples]));
end

function u = reconstructReducedDeltaTau(robot, idx, q_fix_456, q, qd, qdd, N)
u = zeros(3,N);
for k = 1:N
    q_full  = [q(:,k);  q_fix_456];
    qd_full = [qd(:,k); zeros(3,1)];

    D3 = reducedMassMatrix(robot, q_full, idx);
    C3 = reducedVelocityProduct(robot, q_full, qd_full, idx);

    u(:,k) = D3*qdd(:,k) + C3;
end
end

function stage = computeStageCost(q, qd, u, qg, qdg, Qx, Ru)
N = size(u,2);
stage = zeros(N,1);

for k = 1:N
    xk = [q(:,k)-qg; qd(:,k)-qdg];
    uk = u(:,k);
    stage(k) = 0.5*(xk.'*Qx*xk) + 0.5*(uk.'*Ru*uk);
end
end
