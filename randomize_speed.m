function x_new = randomize_speed(x)
varMax = 0.7;
rand_speed = [-varMax+rand()*2*varMax -varMax+rand()*2*varMax];
x_new = [x(1); x(2)+rand_speed(1); x(3); x(4)+rand_speed(2)];
end