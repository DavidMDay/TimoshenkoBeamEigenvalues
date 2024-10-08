% MATLAB script generating a figure
load('param.data');
load('iter.data');
k = param(:,1);
plot( k, iter, 'o--');
z = axis;
z(3) = 0;
axis(z);
xlabel('Radius of Gyration k');
ylabel('Newton Iterations');
title('Continuation problem for 4th mode of free-free Timoshenko beam');
