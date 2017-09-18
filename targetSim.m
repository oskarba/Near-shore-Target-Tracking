x0 = [-300; 4; -200; 5];
cov0 = diag([10, 0.5, 10, 0.5]);
x0_est = x0+chol(cov0)*randn(4,1);
dt = 1;
t_end = 50;
q = 0.25;                          % Process noise strength squared
r = 50;
time = 1:dt:t_end;
F = [1, dt, 0, 0; 0, 1, 0, 0; 0, 0, 1, dt; 0, 0, 0, 1];
H = [1, 0, 0, 0; 0, 0, 1, 0];
R = r*eye(2);
G = [dt^2/2; dt];
Q = q*[G*G', zeros(2); zeros(2), G*G'];
K = length(time);
% Empty state and covariance vectors
x_true = zeros(4,K);
x_est_prior = zeros(4,K);
x_est_posterior = zeros(4,K);
cov_prior = zeros(4,4,K);
cov_posterior = zeros(4,4,K);
z = zeros(2, K);
%-----------------------------------------------
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
	% Generate measurement
	noise = chol(R)*randn(2,1);
	z(:,k) = H*x_true(:,k)+noise;
	% Measurement update
	S = H*cov_prior(:,:,k)*H'+R;              % Covariance of the innovation
	W = cov_prior(:,:,k)*H'/S;                % Gain
    v = z(:,k)-H*x_est_prior(:,k);            % Innovation
	x_est_posterior(:,k) = x_est_prior(:,k)+W*v;
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