clear all;

% Initialized target
x0 = [-300; 4; -200; 5];
cov0 = diag([10, 0.5, 10, 0.5]);
x0_est = x0+chol(cov0)*randn(4,1);

% Time for simulation
dt = 1;
t_end = 50;
time = 1:dt:t_end;

% Area of simulation
x_lim = [-1000 1000];
y_lim = [-1000 1000];
V = (x_lim(2) - x_lim(1))*(y_lim(2) - y_lim(1));     % Area (m^2)

% Kalman filter stuff
q = 0.25;                       % Process noise strength squared
r = 50;
F = [1, dt, 0, 0; 0, 1, 0, 0; 0, 0, 1, dt; 0, 0, 0, 1];
H = [1, 0, 0, 0; 0, 0, 1, 0];
R = r*eye(2);
G = [dt^2/2; dt];
Q = q*[G*G', zeros(2); zeros(2), G*G'];
K = length(time);
gamma_g = 9.21;                 % Gate threshold
c = 2;                          % Normalization constant

%PDA constants
PG = 0.99;
PD = 0.9;

% Empty state and covariance vectors
x_true = zeros(4,K);
x_est_prior = zeros(4,K);
x_est_posterior = zeros(4,K);
cov_prior = zeros(4,4,K);
cov_posterior = zeros(4,4,K);

% Measurement vectors
z_true = zeros(2, K);
z_gate = zeros(2, 10*K);           % Inne i loop? Variabel st�rrelse
z_count = 1;

% Structs
z_strength = struct('x', {});
z_k_gate = struct('x', {});
nr_gate = 1;

% -----------------------------------------------
% Main loop
for k = 1:K
	% Dynamic model
	if k == 1
		x_true(:,k) = x0;
		x_est_prior(:,k) = x0_est;
		cov_prior(:,:,k) = cov0;
	else
		x_true(:,k) = F*x_true(:,k-1);
		x_est_prior(:,k) = F*x_est_posterior(:,k-1);
		cov_prior(:,:,k) = F*cov_posterior(:,:,k-1)*F'+Q;
    end
    
	% Generate measurement for real target
	noise = chol(R)*randn(2,1);
	z(:,k) = H*x_true(:,k)+noise;
    
    % Add clutter measurements and signal strength
    clut_h = 0.001;
    clut_l = 0.00001;
    lambda = clut_l + (clut_h - clut_l)*rand();     % Clutter density
    num_clut = int16(lambda*V);                % Number of clutter points
    z_all = [z(1,k) randi([x_lim(1) x_lim(2)],1,num_clut);
        z(2,k) randi([y_lim(1) y_lim(2)],1,num_clut);
        rand() rand(1,num_clut)];
    
    % Test signal strength
    z_strength = [];
    j = 1;
    for i = 1:length(z_all(1,:))
        if (z_all(3,i) < PD)
            z_strength(j).x = [z_all(1,i); z_all(2,i)];
            j = j + 1;
        end
    end
    
	% Find measurements within validation region
    S = H*cov_prior(:,:,k)*H'+R;             % Covariance of the innovation
    W = cov_prior(:,:,k)*H'/S;               % Gain
    m_k = 0;
    z_k_gate = [];
    beta_i = [0];
    v_i = [0; 0];
    for i = 1:length(z_strength)
        v_ik = z_strength(i).x - H*x_est_prior(:,k);            % Measurement innovation
        NIS_temp = v_ik'/S*v_ik;
        if NIS_temp < gamma_g          % Within validation region
            %z_gate(:,z_count) = z_strength(i).x;
            z_k_gate(length(z_k_gate)+1).x = z_strength(i).x;
            z_count = z_count + 1;
            
            v_i = [v_i v_ik];
            beta_i = [beta_i exp(-0.5*NIS_temp)];
            m_k = m_k + 1;
        end
    end
    beta_i(1) = 2*(1-PD*PG)*m_k/(gamma_g);
    if beta_i ~= 0
        beta_i = beta_i/sum(beta_i);            % Normalize
    end
    
    if m_k ~= 0
        v_k = [dot(v_i(1,:), beta_i);
            dot(v_i(2,:), beta_i)];
    else
        v_k = [0; 0];
    end
    
    % Covariance of correct measurement and SOI
    P_c = cov_prior(:,:,k)-W*S*W';
    part_result = 0;
    for i = 2:length(beta_i)
        part_result = part_result + beta_i(i)*(v_i(:,i)*v_i(:,i)');
    end
    part_result = part_result - v_k*v_k';
    SOI = W*(part_result)*W';

	x_est_posterior(:,k) = x_est_prior(:,k)+W*v_k;
    %cov_posterior(:,:,k) = beta_i(1)*cov_prior(:,:,k)-(1-beta_i(1))*P_c+SOI;
	cov_posterior(:,:,k) = cov_prior(:,:,k)-W*S*W';
end

%--------------------------------------------------
% Plotting
figure; hold on;
plot(x_true(3,:), x_true(1,:), 'k');
plot(x_est_posterior(3,:), x_est_posterior(1,:))
plot(x_est_prior(3,:), x_est_prior(1,:))
legend('True', 'posterior', 'prior');
plot(z(2,:), z(1,:), 'o');
plot(z_k_gate.x(2), z_k_gate.x(1), 'x');
%plot(z_gate(2,1:z_count-1), z_gate(1,1:z_count-1), 'x');
title('North-east position plot');
figure; subplot(2,1,1); hold on;
plot(time, x_true(2,:), 'k');
plot(time, x_est_posterior(2,:));
plot(time, x_est_prior(2,:));
title('North velocity');
subplot(2,1,2); hold on;
plot(time, x_true(4,:), 'k');
plot(time, x_est_posterior(4,:));
plot(time, x_est_prior(4,:));
title('East velocity');
%figure; plot(time,NIS)