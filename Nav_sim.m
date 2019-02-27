tic
% clear all
% for n_test = 0:7
% r = 1;
% for r = 1:100
%% INITIALIZATION
for i = 1:1
    
clear all
addpath(genpath('gnc'));
rng(123908570);

n_test = 0;
n       = 16;
n_err   = 15;

%Setting global constants
Constants.r2d   = 180/pi;
Constants.d2r   = 1/Constants.r2d;
Constants.g_z   = 9.81;
Constants.g     = [0; 0; Constants.g_z];
Constants.h     = 0.01;
Constants.piInc = pi;
% Constants.r     = r;

loadStruct(Constants);

%Time constants
t0      = 0;
tend    = 20;
t       = t0:h:tend;
sim_len = length(t);
h_corr  = 10;
if mod(100,h_corr) ~= 0
    error('Indivisable correction frequency.')
end

Constants.sim_len = sim_len;

%Utility matrices
O3  = zeros(3);
O31 = zeros(3,1);
O13 = zeros(1,3);
I3  = eye(3);
I31 = ones(3,1);
I13 = ones(1,3);

end
%% SWITCHES
for i = 1:1

Plotting.switch = 1;

Plotting.pos    = 1;
Plotting.vel    = 0;
Plotting.ang    = 0;

Plotting.p_hat          = 0;
Plotting.v_hat          = 0;
Plotting.bacc_hat       = 0;
Plotting.ang_hat        = 0;
Plotting.ang_hat_eul    = 0;
Plotting.bars_hat       = 0;
Plotting.acc_hat        = 0;

Switches.Estimate.Bacc          = 1;
Switches.Estimate.Bars          = 1;
Switches.Estimate.Correction    = 1;
Switches.Estimate.Prediction    = 1;

Switches.Measurements.vel       = mod(n_test,2);
Switches.Measurements.mag       = mod(floor(n_test/2),2);
Switches.Measurements.com       = 0;
Switches.Measurements.vel_b     = mod(floor(n_test/4),2);
Switches.Measurements.acc       = 0;

Switches.Path                   = 2;

end
%% MULITCOPTER STRUCT
for i = 1:1
    
Heli.p = O31;
Heli.v = O31;
Heli.theta = O31;
Heli.omega = O31;

Heli.m = 1;
Heli.I = diag([0.2; 0.2; 0.2]); 

end
%% SENSOR DATA
for i = 1:1

%IMU - ADIS16465 at +-500deg dynamic range.
%Accelerometer
Sensor.acc.data         = zeros(3,1);
Sensor.acc.noise        = sqrt(2*10^(-4));
Sensor.acc.bias.data    = zeros(3,1);
Sensor.acc.bias.noise   = sqrt(3.6*10^(-9));

%Angular Rate Sensor
Sensor.ars.data         = zeros(3,1);
Sensor.ars.noise        = sqrt(4.36*10^(-5));
Sensor.ars.bias.data    = zeros(3,1);
Sensor.ars.bias.noise   = sqrt(1.21*10^(-6));


%Aiding sensors
%Position
Sensor.pos.data         = zeros(3,1);
Sensor.pos.noise        = [0.05; 0.05; 0.13];

%Velocity
Sensor.vel.data         = zeros(3,1);
Sensor.vel.noise        = [0.01; 0.01; 0.03];

Sensor.mag.vec_n        = [13.5758; 0.8132; 50.0765];
Sensor.mag.data         = zeros(3,1);
Sensor.mag.noise        = [0.138; 0.089; 0.165];

Sensor.com.vec_n        = [1.00; 0.00; 0.00];
Sensor.com.data         = zeros(3,1);
Sensor.com.noise        = [0.01; 0.01; 0];

Sensor.vel_b.data       = zeros(3,1);
Sensor.vel_b.noise      = [0.2; 0.2; 0.2];

Sensor.fsn.data         = zeros(3,1);
Sensor.fsn.noise        = [0.2; 0.2; 0.1];

Sensor.True_a           = -g;

end
%% ESTIMATOR DATA
for i = 1:1

%Process and measurement noises
w  = [I31*0; I31*Sensor.acc.noise; I31*Sensor.acc.bias.noise*10; I31*Sensor.ars.noise; I31*Sensor.ars.bias.noise*5];
v  = Sensor.pos.noise;

if Switches.Measurements.vel == 1
    v = [v; Sensor.vel.noise];
end
if Switches.Measurements.mag == 1
    v = [v; Sensor.mag.noise];
elseif Switches.Measurements.com == 1
    v = [v; Sensor.com.noise];
end
if Switches.Measurements.vel_b == 1
    v = [v; Sensor.vel_b.noise];
end
if Switches.Measurements.acc == 1
    v = [v; Sensor.fsn.noise];
end

%Initial error
e0 = [I31*0.1;   I31*0.01;   I31*0.001;   I31*0.01*d2r;   I31*0.01*d2r;];
       %p-m       %v-m/s    %bacc-m/s^2     %ang-rad       %bars-rad/s

%Estimator
T_t             = Sensor.acc.bias.noise;
T_r             = Sensor.ars.bias.noise;

%Set initial quaternion
Estimate.x_hat  = zeros(n,1); Estimate.x_hat(10:13) = euler2q(Heli.theta(1), Heli.theta(2), Heli.theta(3));
Estimate.dx     = zeros(n_err,1);
Estimate.F      = @(q_hat, f_hat, w_hat) ...
                   [O3  I3       O3                O3               O3   ;
                    O3  O3    -R(q_hat)  -R(q_hat)*skew(f_hat)      O3   ;
                    O3  O3     -T_t*I3             O3               O3   ;
                    O3  O3       O3           -skew(w_hat)         -I3   ;
                    O3  O3       O3                O3            -T_r*I3];
if Switches.Estimate.Bacc == 0 && Switches.Estimate.Bars == 1
    Estimate.F      = @(q_hat, f_hat, w_hat) ...
                        [O3  I3       O3                O3               O3   ;
                         O3  O3       O3      -R(q_hat)*skew(f_hat)      O3   ;
                         O3  O3       O3                O3               O3   ;
                         O3  O3       O3           -skew(w_hat)         -I3   ;
                         O3  O3       O3                O3            -T_r*I3];
elseif Switches.Estimate.Bacc == 1 && Switches.Estimate.Bars == 0
    Estimate.F      = @(q_hat, f_hat, w_hat) ...
                        [O3  I3       O3                O3               O3   ;
                         O3  O3    -R(q_hat)  -R(q_hat)*skew(f_hat)      O3   ;
                         O3  O3     -T_t*I3             O3               O3   ;
                         O3  O3       O3           -skew(w_hat)          O3   ;
                         O3  O3       O3                O3               O3   ];
elseif Switches.Estimate.Bacc == 0 && Switches.Estimate.Bars == 0
    Estimate.F      = @(q_hat, f_hat, w_hat) ...
                        [O3  I3       O3                O3               O3   ;
                         O3  O3       O3      -R(q_hat)*skew(f_hat)      O3   ;
                         O3  O3       O3                O3               O3   ;
                         O3  O3       O3           -skew(w_hat)          O3   ;
                         O3  O3       O3                O3               O3   ];
end

Estimate.Q      = diag(w.^2);
Estimate.R      = diag(v.^2);
Estimate.P      = diag(e0.^2);
Estimate.H      = zeros(3,n_err);
Estimate.theta  = zeros(3,1);
clear w v e0

end
%% ACCELERATION ESTIMATOR
for i = 1:1
TMO.x_hat = Estimate.x_hat(1:9);
TMO.filtered_acc = zeros(3,20);
T_k = 0.97;
if Switches.Measurements.vel == 1
    TMO.y = zeros(6,1);
else
    TMO.y = zeros(3,1);
end
    
w  = [I31*0; I31*0; I31*0.15];
v  = Sensor.pos.noise;
if Switches.Measurements.vel == 1
    v = [v; Sensor.vel.noise];
end
%Initial error
e0 = [I31*0.1;   I31*0.01;   I31*0.1];
       %p-m       %v-m/s    %acc-m/s^2 

TMO.F = [O3 I3 O3;
         O3 O3 I3;
         O3 O3 O3];

       
TMO.Q = diag(w.^2);
TMO.R = diag(v.^2);
TMO.P = diag(e0.^2);
TMO.F = Van_Loan_discretization(TMO.F, TMO.Q, eye(9), Constants.h);
TMO.H = [I3 O3 O3];
if Switches.Measurements.vel == 1
    TMO.H = [TMO.H; [O3 I3 O3]];
end
clear w v e0
end
%% LOG
for i = 1:1
    
log = Logger();
log.init(t,h);

%States
log.add('p', 3);
log.add('v', 3);
log.add('theta', 3);
log.add('omega', 3);
log.add('bacc', 3);
log.add('bars', 3);

%Desired states
log.add('p_d', 3);
log.add('v_d', 3);
log.add('theta_d', 3);

%Estimated states
log.add('x_hat', n);
log.add('theta_hat', 3);
log.add('theta_imu', 3);
log.add('dx', n_err);
log.add('P', n_err);

%Measurements
log.add('true_a', 3);
log.add('fsn', 3);
log.add('fsn_lp', 3);
log.add('P_TMO', 9);
end
%% GENERATE PATH
for i = 1:1
    
if Switches.Path == 0
    Path_p          = [zeros(1,sim_len); zeros(1,sim_len); zeros(1,sim_len)]*10;
elseif Switches.Path == 1
    Path_p          = [cos((2*pi/5)*t); sin((2*pi/5)*t); sin((2*pi/5)*t) + cos((2*pi/5)*t)];
    Path_p          = Path_p - [cos((2*pi/10)*t); sin((2*pi/10)*t); sin((2*pi/10)*t) + cos((2*pi/10)*t)]*0.5;
elseif Switches.Path == 2
    Path_p          = [cos((2*pi/5)*t); sin((2*pi/5)*t); sin((2*pi/5)*t) + cos((2*pi/5)*t)]*2;
    Path_p          = Path_p - [cos((2*pi/10)*t); sin((2*pi/10)*t); sin((2*pi/10)*t) + cos((2*pi/10)*t)];      
end

Path_v          = [zeros(3,1), diff(Path_p, 1, 2)]/h;
Desired.p       = Heli.p;
Desired.v       = Heli.v;
Desired.theta   = Heli.theta;

end
%% SIMULATION
for k = 1:sim_len
    
    %Breakpoint at specific iteration
    if mod(k,(sim_len-1)/2) == 0
%         Switches.Measurements.vel = 0;
%         Switches.Measurements.mag = 0;
%         Estimate.R = diag([Sensor.pos.noise].^2);
%         if mod(floor(n_test/4),2) == 1
%             Estimate.R = diag([Sensor.pos.noise; Sensor.vel_b.noise].^2);
%         end
    end
    
    %Storing of variables
    log.store('p', Heli.p, k);
    log.store('v', Heli.v, k);
    log.store('theta', Heli.theta, k);
    log.store('omega', Heli.omega, k);
    log.store('bacc', Sensor.acc.bias.data, k);
    log.store('bars', Sensor.ars.bias.data, k);
    
    log.store('true_a', Sensor.True_a, k);
    log.store('fsn', TMO.x_hat(7:9)-g, k);
    log.store('fsn_lp', TMO.filtered_acc(:,end)-g, k);
    log.store('P_TMO', diag(TMO.P), k);
    
    log.store('p_d', Desired.p, k);
    log.store('v_d', Desired.v, k);
    log.store('theta_d', Desired.theta, k);
    
    log.store('x_hat', Estimate.x_hat, k);
    log.store('theta_hat', q2euler(Estimate.x_hat(10:13)), k);
    log.store('theta_imu', Estimate.theta, k);
    log.store('P', diag(Estimate.P), k);
    log.store('dx', Estimate.dx, k)
    
    %Generating desired states
    Desired.p   = Desired.p + 0.5*h*(Path_p(:,k) - Desired.p);
    Desired.v   = Desired.v + 0.5*h*(Path_v(:,k) - Desired.v);
    Desired.theta(3) = sin(t(k)*2*pi/10)*pi/4;
    
    [Heli, Sensor, Desired] = state_update(Heli, Sensor, Desired, h, Constants, Switches);
    
    if Switches.Estimate.Prediction == 1
        Estimate = navigation_prediction(Estimate, Sensor, Constants);
        if Switches.Measurements.acc == 1
            TMO = TMO_prediction(TMO, Constants);
        end
    end
        
    if Switches.Estimate.Correction == 1 && mod(k,100/h_corr) == 0
        Sensor = generate_measurements(Heli, Estimate, Sensor, Switches);
        if Switches.Measurements.acc == 1
            if Switches.Measurements.vel == 1
                TMO.y = [Sensor.pos.data; Sensor.vel.data];
            else
                TMO.y = Sensor.pos.data;
            end
            TMO = TMO_correction(TMO, Constants, Switches);
            TMO.filtered_acc = one_step_lp_filtering(TMO.x_hat(7:9), TMO.filtered_acc, k, T_k);
            
            Sensor.y(end-2:end) = TMO.filtered_acc(:,end);
            Estimate.R(end-2:end,end-2:end) = TMO.P(7:9,7:9);
        end
        Estimate.H = measurement_function(Estimate, Sensor, Switches);
        Estimate = extended_kalman_filter(Estimate, Sensor, Constants, Switches);

        Estimate.x_hat(1:9) = Estimate.x_hat(1:9) + Estimate.dx(1:9);
        Estimate.x_hat(14:16) = Estimate.x_hat(14:16) + Estimate.dx(13:15);

        a_g = Estimate.dx(10:12);
        dq  = (1/sqrt(4+norm(a_g)^2))*[2; a_g];
        Estimate.x_hat(10:13) = qprod(Estimate.x_hat(10:13), dq);
    end
end
if Plotting.switch == 1
    print(log,t,Constants, Plotting);
end
% end
% load('results')
% results.final_results.pos(:,n_test+1) = mean(results.pos,2);
% results.final_results.vel(:,n_test+1) = mean(results.vel,2);
% results.final_results.angles(:,n_test+1) = mean(results.angles,2);
% save('results','results');
% end
% results.final_results.pos'
% results.final_results.vel'
% results.final_results.angles'
% toc
%% FUNCTIONS

function print(log, t, Constants, Plotting)
    %Init
    for k = 1:1
        fig = 1;
        loadStruct(Constants);
        p           = log.get('p', ':');
        v           = log.get('v', ':');
        theta       = log.get('theta', ':');
        bacc        = log.get('bacc', ':');
        bars        = log.get('bars', ':');
        
        p_d         = log.get('p_d', ':');
        v_d         = log.get('v_d', ':');
        theta_d     = log.get('theta_d', ':');
        
        x_hat       = log.get('x_hat', ':');
        p_hat       = x_hat(1:3,:);
        v_hat       = x_hat(4:6,:);
        bacc_hat    = x_hat(7:9,:);
        q_hat       = x_hat(10:13,:);
        theta_hat   = log.get('theta_hat', ':');
        theta_e     = log.get('theta_imu', ':');
        bars_hat    = x_hat(14:16,:);
        P           = sqrt(log.get('P', ':'));
        
        fsn         = log.get('fsn_lp', ':');
        true_a      = log.get('true_a', ':');
        P_TMO       = log.get('P_TMO',':');
        
%         acceleration = [var((fsn-true_a)')', mean((fsn-true_a)')']
%         halfway = (sim_len-1)/2;
% %         ang_err = r2d*(theta-theta_hat);
% %         angles = sqrt([var(ang_err(:,1:halfway),0,2), var(ang_err(:,halfway+1:end),0,2)])
% 
%         load('results');
%         halfway = 2500;
%         ang_err = r2d*(theta'-theta_hat');
%         ang_err = sqrt([var(ang_err(1:halfway,:))', var(ang_err(halfway+1:end,:))']);
%         results.angles(1:3,r) = ang_err(:,1);
%         results.angles(4:6,r) = ang_err(:,2);
%         
%         pos_err = (p'-p_hat');
%         pos_err = sqrt([var(pos_err(1:halfway,:))', var(pos_err(halfway+1:end,:))']);
%         results.pos(1:3,r) = pos_err(:,1);
%         results.pos(4:6,r) = pos_err(:,2);
%         
%         vel_err = (v'-v_hat');
%         vel_err = sqrt([var(vel_err(1:halfway,:))', var(vel_err(halfway+1:end,:))']);
%         results.vel(1:3,r) = vel_err(:,1);
%         results.vel(4:6,r) = vel_err(:,2);
%         
%         save('results', 'results');
        
    end
    
    %Position
    if Plotting.pos == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,1,1)
            plot(t,p(1,:),'r');
            hold on
            plot(t,p_d(1,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Position')
            legend({'$p_x$', '$p_{d_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,2)
            plot(t,p(2,:),'r');
            hold on
            plot(t,p_d(2,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$p_y$', '$p_{d_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,3)
            plot(t,p(3,:),'r');
            hold on
            plot(t,p_d(3,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$p_z$', '$p_{d_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %Velocity
    if Plotting.vel == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;

            subplot(3,1,1)
            plot(t,v(1,:),'r');
            hold on
            plot(t,v_d(1,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Velocity')
            legend({'$v_x$', '$v_{d_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,2)
            plot(t,v(2,:),'r');
            hold on
            plot(t,v_d(2,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$v_y$', '$v_{d_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,3)
            plot(t,v(3,:),'r');
            hold on
            plot(t,v_d(3,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$v_z$', '$v_{d_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %Angles
    if Plotting.ang == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;

            subplot(3,1,1)
            plot(t,r2d*theta(1,:),'r');
            hold on
            plot(t,r2d*theta_d(1,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Angles')
            legend({'$\phi$', '$\phi_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,2)
            plot(t,r2d*theta(2,:),'r');
            hold on
            plot(t,r2d*theta_d(2,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\theta$', '$\theta_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,3)
            plot(t,r2d*theta(3,:),'r');
            hold on
            plot(t,r2d*theta_d(3,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\psi$', '$\psi_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %State Estimates: Position
    if Plotting.p_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,2,1)
            plot(t,p_hat(1,:),'r');
            hold on
            plot(t,p(1,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Position Estimates')
            legend({'$\hat{p}_x$', '$p_x$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,2)
            plot(t,p(1,:)-p_hat(1,:), 'r')
            hold on
            plot(t,-P(1,:), 'k--')
            plot(t,P(1,:), 'k--')
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Position Estimation Errors')
            legend({'$\tilde{p}_x$', '$\sigma_{p_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,p_hat(2,:),'r');
            hold on
            plot(t,p(2,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{p}_y$', '$p_y$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,4)
            plot(t,p(2,:)-p_hat(2,:), 'r')
            hold on
            plot(t,-P(2,:), 'k--')
            plot(t,P(2,:), 'k--')
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{p}_y$', '$\sigma_{p_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,p_hat(3,:),'r');
            hold on
            plot(t,p(3,:),'k--');
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{p}_z$', '$p_z$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,p(3,:)-p_hat(3,:), 'r')
            hold on
            plot(t,-P(3,:), 'k--')
            plot(t,P(3,:), 'k--')
            hold off
            grid on
            ylabel('$\{m\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{p}_z$', '$\sigma_{p_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %State Estimate: Velocity
    if Plotting.v_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,2,1)
            plot(t,v_hat(1,:),'r');
            hold on
            plot(t,v(1,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Velocity Estimates')
            legend({'$\hat{v}_x$', '$v_x$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,2)
            plot(t,v(1,:)-v_hat(1,:), 'r')
            hold on
            plot(t,-P(4,:), 'k--')
            plot(t,P(4,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Velocity Estimation Errors')
            legend({'$\tilde{v}_x$', '$\sigma_{v_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,v_hat(2,:),'r');
            hold on
            plot(t,v(2,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{v}_y$', '$v_y$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,4)
            plot(t,v(2,:)-v_hat(2,:), 'r')
            hold on
            plot(t,-P(5,:), 'k--')
            plot(t,P(5,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{v}_y$', '$\sigma_{v_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,v_hat(3,:),'r');
            hold on
            plot(t,v(3,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{v}_z$', '$v_z$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,v(3,:)-v_hat(3,:), 'r')
            hold on
            plot(t,-P(6,:), 'k--')
            plot(t,P(6,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{v}_z$', '$\sigma_{v_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %State Estimate: Acceleration bias
    if Plotting.bacc_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,2,1)
            plot(t,bacc_hat(1,:),'r');
            hold on
            plot(t,bacc(1,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Accelerometer Bias Estimates')
            legend({'$\hat{b}_{acc_x}$', '$b_{acc_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,2)
            plot(t,bacc(1,:)-bacc_hat(1,:), 'r')
            hold on
            plot(t,-P(7,:), 'k--')
            plot(t,P(7,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Acc Bias Estimation Errors')
            legend({'$\tilde{b}_{acc_x}$', '$\sigma_{b_{acc}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,bacc_hat(2,:),'r');
            hold on
            plot(t,bacc(2,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{b}_{acc_y}$', '$b_{acc_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,4)
            plot(t,bacc(2,:)-bacc_hat(2,:), 'r')
            hold on
            plot(t,-P(8,:), 'k--')
            plot(t,P(8,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{b}_{acc_y}$', '$\sigma_{b_{acc}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,bacc_hat(3,:),'r');
            hold on
            plot(t,bacc(3,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{b}_{acc_z}$', '$b_{acc_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,bacc(3,:)-bacc_hat(3,:), 'r')
            hold on
            plot(t,-P(9,:), 'k--')
            plot(t,P(9,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{b}_{acc_z}$', '$\sigma_{b_{acc}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %State Estimate: Angles from quaternions
    if Plotting.ang_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;

            subplot(3,2,1)
            plot(t,r2d*theta_hat(1,:),'r');
            hold on
            plot(t,r2d*theta(1,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Angle Estimates')
            legend({'$\hat{\phi}$', '$\phi$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,2)
            plot(t,r2d*(theta(1,:)-theta_hat(1,:)), 'r')
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Angle Estimation Errors')
            legend({'$\tilde{\phi}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,r2d*theta_hat(2,:),'r');
            hold on
            plot(t,r2d*theta(2,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{\theta}$', '$\theta$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,4)
            plot(t,r2d*(theta(2,:)-theta_hat(2,:)), 'r')
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{\theta}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,r2d*theta_hat(3,:),'r');
            hold on
            plot(t,r2d*theta(3,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{\psi}$', '$\psi$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,r2d*(theta(3,:)-theta_hat(3,:)), 'r')
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{\psi}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

        end
    end
    
    %State Estimate: Angles from euler
    if Plotting.ang_hat_eul == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;

            subplot(3,1,1)
            plot(t,r2d*theta_e(1,:),'r');
            hold on
            plot(t,r2d*theta(1,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Euler Angle Estimates')
            legend({'$\phi$', '$\phi_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,2)
            plot(t,r2d*theta_e(2,:),'r');
            hold on
            plot(t,r2d*theta(2,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\theta$', '$\theta_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,1,3)
            plot(t,r2d*theta_e(3,:),'r');
            hold on
            plot(t,r2d*theta(3,:),'k--');
            hold off
            grid on
            ylabel('$\{deg\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\psi$', '$\psi_d$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %State Estimate: Gyro bias
    if Plotting.bars_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,2,1)
            plot(t,bars_hat(1,:),'r');
            hold on
            plot(t,bars(1,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{rad}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Gyrometer Bias Estimates')
            legend({'$\hat{b}_{ars_x}$', '$b_{ars_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,2)
            plot(t,bars(1,:)-bars_hat(1,:), 'r')
            hold on
            plot(t,-P(13,:), 'k--')
            plot(t,P(13,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{rad}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Gyro Bias Estimation Error')
            legend({'$\tilde{b}_{ars_x}$', '$\sigma_{b_{ars}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,bars_hat(2,:),'r');
            hold on
            plot(t,bars(2,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{rad}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{b}_{ars_y}$', '$b_{ars_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,4)
            plot(t,bars(2,:)-bars_hat(2,:), 'r')
            hold on
            plot(t,-P(14,:), 'k--')
            plot(t,P(14,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{rad}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{b}_{ars_y}$', '$\sigma_{b_{ars}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,bars_hat(3,:),'r');
            hold on
            plot(t,bars(3,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{b}_{ars_z}$', '$b_{ars_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,bars(3,:)-bars_hat(3,:), 'r')
            hold on
            plot(t,-P(15,:), 'k--')
            plot(t,P(15,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{rad}{s}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{b}_{ars_z}$', '$\sigma_{b_{ars}}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
    
    %Estimated Acceleration
    if Plotting.acc_hat == 1
        for k = 1:1
            figure(fig)
            fig = fig + 1;
            subplot(3,2,1)
            plot(t,fsn(1,:),'r');
            hold on
            plot(t,true_a(1,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Acceleration Estimates')
            legend({'$\hat{a}_x$', '$a_x$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,2)
            plot(t,fsn(1,:)-true_a(1,:), 'r')
            hold on
            plot(t,-P_TMO(7,:), 'k--')
            plot(t,P_TMO(7,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            title('Acceleration Estimation Errors')
            legend({'$\tilde{a}_x$', '$\sigma_{a_x}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,3)
            plot(t,fsn(2,:),'r');
            hold on
            plot(t,true_a(2,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{a}_y$', '$a_y$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,4)
            plot(t,fsn(2,:)-true_a(2,:), 'r')
            hold on
            plot(t,-P_TMO(8,:), 'k--')
            plot(t,P_TMO(8,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{a}_y$', '$\sigma_{a_y}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')

            subplot(3,2,5)
            plot(t,fsn(3,:),'r');
            hold on
            plot(t,true_a(3,:),'k--');
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\hat{a}_z$', '$a_z$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
            
            subplot(3,2,6)
            plot(t,fsn(3,:)-true_a(3,:), 'r')
            hold on
            plot(t,-P_TMO(9,:), 'k--')
            plot(t,P_TMO(9,:), 'k--')
            hold off
            grid on
            ylabel('$\{\frac{m}{s^2}\}$', 'Interpreter', 'latex', 'FontSize', 20)
            xlabel('Time $\{s\}$', 'Interpreter', 'latex', 'FontSize', 20)
            legend({'$\tilde{a}_z$', '$\sigma_{a_z}$'}, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northeastoutside')
        end
    end
end

function q = qprod(q1, q2)
    if any(~isreal(q1(:)))
        error('First input elements are not real.');
    end
    if any(size(q1)~=[4,1])
        error('First input dimension is not 4-by-1.');
    end
    if any(~isreal(q2(:)))
        error('Second input elements are not real.');
    end
    if any(size(q2)~=[4,1])
        error('Second input dimension is not 4-by-1.');
    end

    eta1 = q1(1);
    eps1 = q1(2:4);
    eta2 = q2(1);
    eps2 = q2(2:4);

    eta = eta1*eta2 - eps1'*eps2;
    eps = eta1*eps2 + eta2*eps1 + skew(eps1)*eps2;

    q = [eta; eps];
end

function S = skew(vec)
    S = [      0, -vec(3),  vec(2);
          vec(3),       0, -vec(1);
         -vec(2),  vec(1),       0];
end

function Rmat = R(att)
    if length(att) == 3
        Rmat = Rzyx(att(1), att(2), att(3));
    elseif length(att) == 4
        Rmat = Rquat(att);
    else
        error('Invalid length of input argument');
    end
end

function filtered_sig = one_step_lp_filtering(sig, filtered_sig, k, T_k)
    if T_k > 1 && T_k < 0
        error('Invalid T_k');
    end
    order = size(filtered_sig,2);
    filtered_sig(:,1) = filtered_sig(:,1) + T_k*(sig - filtered_sig(:,1));
    for i = 2:order
        filtered_sig(:,i) = filtered_sig(:,i) + T_k*(filtered_sig(:,i-1) - filtered_sig(:,i));
    end
end

function Sensor = generate_measurements(Heli, Estimate, Sensor, Switches)
    %GPosition Measurements
    Sensor.pos.data = Heli.p + randn(3,1).*Sensor.pos.noise;
    Sensor.y = [Sensor.pos.data];
    
    %Velocity ;easurements
    if Switches.Measurements.vel == 1
        Sensor.vel.data = Heli.v + randn(3,1).*Sensor.vel.noise;
        Sensor.y = [Sensor.y; Sensor.vel.data];
    end
    %Magnetometer/Compass Measurements
    if Switches.Measurements.mag == 1
        Sensor.mag.data = R(Heli.theta)'*(Sensor.mag.vec_n);% + randn(3,1).*Sensor.mag.noise);
        Sensor.y = [Sensor.y; Sensor.mag.vec_n];
    elseif Switches.Measurements.com == 1
        Sensor.com.data = Rzyx(0,0,Heli.theta(3))'*(Sensor.com.vec_n + randn(3,1).*Sensor.com.noise);
        Sensor.y = [Sensor.y; Sensor.com.data];
    end
    
    %Optical flow
    if Switches.Measurements.vel_b == 1
        Sensor.vel_b.data = Sensor.True_v_b + randn(3,1) .* Sensor.vel_b.noise;
        Sensor.y = [Sensor.y; Estimate.x_hat(4:6)];
    end
    
    %Specific Force Measurements
    if Switches.Measurements.acc == 1
        Sensor.fsn.data = Sensor.True_a + randn(3,1).*Sensor.fsn.noise;
        Sensor.y = [Sensor.y; Sensor.fsn.data];
    end

end

function H = measurement_function(Estimate, Sensor, Switches)
    I3 = eye(3); O3 = zeros(3); O31 = zeros(3,1); Rmat = R(Estimate.x_hat(10:13));
    %Position Measurements
    H  = [I3 O3 O3 O3 O3];
       
    %Velocity Measurements
    if Switches.Measurements.vel == 1
        H1 = [O3 I3 O3 O3 O3];
        H  = [H; H1];
    end
    %Magnetometer Measurements
    if Switches.Measurements.mag == 1
        H2 = [O3 O3 O3 -Rmat*skew(Sensor.mag.data) O3];
        H  = [H; H2];
    elseif Switches.Measurements.com == 1
        H2 = [O3 O3 O3 [O31 O31 -skew(Sensor.com.data)*[0;0;1]] O3];
        H  = [H; H2];
    end
    
    %Optical Flow
    if Switches.Measurements.vel_b == 1
        H3 = [O3 O3 O3 -skew(Sensor.vel_b.data) O3];
        H = [H; H3];
    end
    
    %Specific Force Measurements
    if Switches.Measurements.acc
        H4 = [O3 O3 O3 -Rmat*skew(Sensor.acc.data - Estimate.x_hat(7:9)) O3];
        H = [H; H4];
    end
end

function Estimate = navigation_prediction(Estimate, Sensor, Constants)
    loadStruct(Constants);

    x                       = Estimate.x_hat;
    u                       = [Sensor.acc.data - x(7:9); Sensor.ars.data - x(14:16)];
    [~, R_mat, T]           = quatern(x(10:13));
    F                       = Estimate.F(x(10:13), u(1:3), u(4:6));
    [Fd, Qd]                = Van_Loan_discretization(F, Estimate.Q, eye(15), h);
    ang                     = Estimate.theta;
    [~, ~, T_eul]   = eulerang(ang(1), ang(2), ang(3));
    
    f = @(x,u)(                 [      x(4:6);
                                  R_mat*u(1:3) + g;
                                      zeros(3,1);
                                       T*u(4:6);
                                      zeros(3,1);]);
    Estimate.theta          = ang + h*T_eul*u(4:6);
    %Perform state prediction 
    Estimate.x_hat(1:9)     = [eye(9), zeros(9,7)]*(x + h*f(x,u));
    Estimate.x_hat(14:16)   = [zeros(3,13), eye(3)]*(x + h*f(x,u));
    
    if norm(u(4:6)) > 0
        dTheta                  = u(4:6)*h;
        dquatern                = [cos(norm(dTheta)/2); dTheta/norm(dTheta)*sin(norm(dTheta)/2)];
        Estimate.x_hat(10:13)   = qprod(x(10:13), dquatern);
        Estimate.x_hat(10:13)   = Estimate.x_hat(10:13)/norm(Estimate.x_hat(10:13));
    end
    
    Estimate.P              = Fd*Estimate.P*Fd' + Qd;
end

function TMO = TMO_prediction(TMO, Constants)
    loadStruct(Constants);
    
    x = TMO.x_hat; P = TMO.P; F = TMO.F; Q = TMO.Q;
    
    x = F*x;
    P = F*P*F' + Q;
    TMO.x_hat = x; TMO.P = P;
end

function Estimate = extended_kalman_filter(Estimate, Sensor, Constants, Switches)
    loadStruct(Constants);
    R_mat = R(Estimate.x_hat(10:13));
    if Switches.Measurements.com == 1
        Rz = atan2(R_mat(2,1), R_mat(1,1));
    end
    P = Estimate.P; H = Estimate.H; Rd = Estimate.R;
    %Measurements
    y = Sensor.y;
    y_hat = [Estimate.x_hat(1:3)];
    
    if Switches.Measurements.vel == 1
        y_hat = [y_hat; Estimate.x_hat(4:6)];
    end
    
    if Switches.Measurements.mag == 1
        y_hat = [y_hat; R_mat*Sensor.mag.data];
    elseif Switches.Measurements.com == 1
        y_hat = [y_hat; Rz*Sensor.com.data];
    end
    if Switches.Measurements.vel_b == 1
        y_hat = [y_hat; R_mat*Sensor.vel_b.data];
    end
    if Switches.Measurements.acc == 1
        y_hat = [y_hat; R_mat*(Sensor.acc.data-Estimate.x_hat(7:9))];
    end
    y_tilde = y-y_hat;
    
    %Correction
    K = P*H'/(H*P*H' + Rd);
    dx = K*y_tilde;
    P = (eye(15) - K*H)*P*(eye(15) - K*H)' + K*Rd*K';
    
    Estimate.P = P_check(P, 3);
    Estimate.dx = dx;
end

function TMO = TMO_correction(TMO, Constants, Switches)
    loadStruct(Constants);

    x = TMO.x_hat; P = TMO.P; H = TMO.H; Rd = TMO.R;
    
    y = TMO.y;
    
    y_hat = x(1:3);
    if Switches.Measurements.vel == 1
        y_hat = [y_hat; x(4:6)];
    end
    y_tilde = y-y_hat;
    K = P*H'/(H*P*H' + Rd);
    dx = K*y_tilde;
    
    P = (eye(9) - K*H)*P*(eye(9) - K*H)' + K*Rd*K';
    x = x + dx;
    
    TMO.P = P_check(P, 3);  TMO.x_hat = x;
end

function [u_b, M, Desired] = controller(Heli, Desired, p_hat, v_hat, theta_hat, Constants)
    loadStruct(Constants);
    
    %Transitional Regulator Parameters
    omega_n_t           = 4;
    zeta_t              = 1.7;
    k_p_t               = Heli.m * omega_n_t^2;
    k_d_t               = Heli.m * omega_n_t * 2 * zeta_t;

    K_p_t               = diag([k_p_t, k_p_t, k_p_t]);
    K_d_t               = diag([k_d_t, k_d_t, k_d_t]);
    
    %Angular Regulator Parameters
    omega_n_r           = 4*omega_n_t;
    zeta_r              = 1;

    k_p_r               = omega_n_r^2;
    k_d_r               = omega_n_r * 2 * zeta_r;

    K_p_r               = diag([k_p_r, k_p_r, k_p_r]);
    K_d_r               = diag([k_d_r, k_d_r, k_d_r]);
    

    [~, ~, T]           = eulerang(Heli.theta(1), Heli.theta(2), Heli.theta(3));

    % Transitional PD-controller
    p_tilde             = (p_hat-Desired.p);
    v_tilde             = (v_hat-Desired.v);
    u_d                 = -Heli.m * (K_d_t*v_tilde + K_p_t*p_tilde);
    
    % Limitations to speed and desired acceleration
    sat                 = @(x, max_a) min(max(x, -max_a), max_a);
    alpha_t             = sat(u_d, abs(Heli.m * g(3)/1.4));
    for i = 1:3
        if Heli.v(i) > 10
            alpha_t(i) = min(alpha_t(i), 0);
        elseif Heli.v(i) < -10
            alpha_t(i) = max(alpha_t(i), 0);
        end
    end
    alpha_t     = -alpha_t + Heli.m * g;
    
    % Angular PD-controller
    u_b                = -alpha_t(3) / (cos(Heli.theta(1)) * cos(Heli.theta(2)));
    
    
    alpha_t            = Rzyx(0, 0, theta_hat(3))'*alpha_t;
    Desired.theta(2)   = atan2(alpha_t(1), alpha_t(3));
    Desired.theta(1)   = -atan2(alpha_t(2), norm([alpha_t(1), alpha_t(3)]) );
    
    alpha_r             = -K_d_r*Heli.omega - T'*K_p_r*(theta_hat-Desired.theta);
    M                   = -skew(Heli.I*Heli.omega)*Heli.omega + alpha_r;
end

function [Heli, Sensor, Desired] = state_update(Heli, Sensor, Desired, h, Constants, Switches)
    loadStruct(Constants);

    [u_b, M, Desired] = controller(Heli, Desired, Heli.p, Heli.v, Heli.theta, Constants);
    
    %Diff equations
    [~, R_mat, T]   = eulerang(Heli.theta(1), Heli.theta(2), Heli.theta(3));
    p_dot       = Heli.v;
    v_dot       = (R_mat * [0;0;u_b])/Heli.m + g;
    theta_dot   = T * Heli.omega;
    omega_dot   = Heli.I \ (skew(Heli.I*Heli.omega) * Heli.omega + M);
    
    %Euler Integration
    Heli.p      = Heli.p + p_dot * h;
    Heli.v      = Heli.v + v_dot * h;
    Heli.theta  = Heli.theta + theta_dot * h;
    Heli.omega  = Heli.omega + omega_dot * h;
    
    %IMU measurements
    
    Sensor.acc.bias.data = Sensor.acc.bias.data + Sensor.acc.bias.noise * randn(3,1);
                              
    Sensor.ars.bias.data = Sensor.ars.bias.data + Sensor.ars.bias.noise * randn(3,1);                                                         
    
    Sensor.acc.data         = R_mat' * (v_dot - g + Sensor.acc.bias.data) + randn(3,1) * Sensor.acc.noise;
    Sensor.ars.data         = Heli.omega + Sensor.ars.bias.data + randn(3,1) * Sensor.ars.noise;
    Sensor.True_a           = v_dot - g;
    Sensor.True_a_b         = R_mat' * Sensor.True_a;
    Sensor.True_v_b         = R_mat' * Heli.v;
end