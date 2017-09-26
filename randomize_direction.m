function x_new = randomize_direction(x)
[theta,rho] = cart2pol(x(2),x(4));
varTheta = 0.2;
theta = theta -varTheta+rand()*2*varTheta;
[x1,y1] = pol2cart(theta,rho);
x_new = [x(1); x1; x(3); y1];
end