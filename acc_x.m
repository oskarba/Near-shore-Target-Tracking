function x_new = acc_x(x, accx, accy)
global dt
x_new = [x(1); x(2) +accx*dt; x(3); x(4)+accy*dt];
end