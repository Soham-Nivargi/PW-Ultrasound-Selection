qm0 = load('Queue-based/Exact paper method/Results/Simulation/25_unwindowed.mat');
qm1 = load('Queue-based/My method/Results/Simulation/25_unwindowed.mat');
qm1_win = load('Queue-based/My method/Results/Simulation/25_windowed.mat');
qm2 = load('Queue-based/My method with modifications to Thresholding/Results/Simulation/25_unwindowed.mat');
qm2_win = load('Queue-based/My method with modifications to Thresholding/Results/Simulation/25_windowed.mat');

nqm0 = load('Not-queue-based/Exact paper method/Results/Simulation/25_unwindowed.mat');
nqm1 = load('Not-queue-based/My method/Results/Simulation/25_unwindowed.mat');
nqm1_win = load('Not-queue-based/My method/Results/Simulation/25_windowed.mat');
nqm2 = load('Not-queue-based/My method with modifications to Thresholding/Results/Simulation/25_unwindowed.mat');
nqm2_win = load('Not-queue-based/My method with modifications to Thresholding/Results/Simulation/25_windowed.mat');

figure();
plot(qm0.j, qm0.con_angle, 'r', 'LineWidth', 2, 'DisplayName', 'Queue-Paper'); hold on;
plot(qm1.j, qm1.con_angle, 'g', 'LineWidth', 2, 'DisplayName', 'Queue-mod unwindowed');
plot(qm1_win.j, qm1_win.con_angle, 'b', 'LineWidth', 2, 'DisplayName', 'Queue-mod windowed');
plot(nqm0.j, nqm0.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'Unqueue-Paper');
plot(nqm1.j, nqm1.con_angle, 'm', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod unwindowed');
plot(nqm1_win.j, nqm1_win.con_angle, 'y', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod windowed');
hold off;

legend;
xlabel('Angle Limit');
ylabel('CNR');

figure();
plot(qm0.j, qm0.con_angle, 'r', 'LineWidth', 2, 'DisplayName', 'Queue-Paper'); hold on;
plot(qm1.j, qm1.con_angle, 'g', 'LineWidth', 2, 'DisplayName', 'Queue-mod unwindowed');
% plot(qm1_win.j, qm1_win.con_angle, 'b', 'LineWidth', 2, 'DisplayName', 'Queue-mod windowed');
plot(nqm0.j, nqm0.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'Unqueue-Paper');
plot(nqm1.j, nqm1.con_angle, 'm', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod unwindowed');
% plot(nqm1_win.j, nqm1_win.con_angle, 'y', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod windowed');
hold off;

legend;
xlabel('Angle Limit');
ylabel('CNR');

% figure();
% % plot(qm0.j, qm0.con_angle, 'r', 'LineWidth', 2, 'DisplayName', 'Queue-Paper'); hold on;
% plot(qm1.j, qm1.con_angle, 'r', 'LineWidth', 2, 'DisplayName', 'Queue-mod unwindowed');
% plot(qm1_win.j, qm1_win.con_angle, 'g', 'LineWidth', 2, 'DisplayName', 'Queue-mod windowed');
% % plot(nqm0.j, nqm0.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'Unqueue-Paper');
% plot(nqm1.j, nqm1.con_angle, 'b', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod unwindowed');
% plot(nqm1_win.j, nqm1_win.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod windowed');
% % plot(qm2.j, qm2.con_angle, 'm', 'LineWidth', 2, 'DisplayName', 'Queue-mod cumu unwindowed');
% % plot(qm2_win.j, qm2_win.con_angle, 'y', 'LineWidth', 2, 'DisplayName', 'Queue-mod cumu windowed');
% % plot(nqm0.j, nqm0.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'Unqueue-Paper');
% plot(nqm2.j, nqm2.con_angle, 'k', 'LineWidth', 2, 'DisplayName', 'UnQueue-mod cumu unwindowed');
% % plot(nqm0.j, nqm0.con_angle, 'c', 'LineWidth', 2, 'DisplayName', 'Unqueue-Paper');
% plot(nqm2_win.j, nqm2_win.con_angle, 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'UnQueue-mod cumu windowed');
% hold off;
% 
% legend;
% xlabel('Angle Limit');
% ylabel('CNR');