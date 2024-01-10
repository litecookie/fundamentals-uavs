%% a line to call simulink model\
qd11Array = out.quadrotor;
% Create a 3D plot using the values from 'muthuArray'
figure;
plot3(qd11Array(:,2), qd11Array(:,1), qd11Array(:,3), 'o-', 'LineWidth', 2);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Plot of Quadrotor Data');
grid on;