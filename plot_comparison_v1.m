% Generate comparison plots for selected cases
clear all; close all; clc;

load PanelvExact_v1.mat
% 1st col: Mach, 2nd col: Panel, 3rd col: Exact

%% Plots
figure(1)
plot(disk_perp(:,1),disk_perp(:,3),'r')
hold on
plot(sphere(:,1),sphere(:,3),'b')
plot(sharp_cone_5_deg(:,1),sharp_cone_5_deg(:,3),'m')
plot(circular_cylinder(:,1),circular_cylinder(:,3),'g')
plot(disk_perp(:,1),disk_perp(:,3),'r')
plot(disk_perp(:,1),disk_perp(:,2),'rs')
plot(sphere(:,1),sphere(:,3),'b')
plot(sphere(:,1),sphere(:,2),'bs')
plot(sharp_cone_5_deg(:,1),sharp_cone_5_deg(:,3),'m')
plot(sharp_cone_5_deg(:,1),sharp_cone_5_deg(:,2),'ms')
plot(circular_cylinder(:,1),circular_cylinder(:,3),'g')
plot(circular_cylinder(:,1),circular_cylinder(:,2),'gs')
hold off
xlabel('Mach number, nd')
ylabel('Drag coefficient')
ylim([1.9 2.2])
grid on
legend('Disk (perp.)','Sphere','Sharp-Cone (5 deg.)','Circular Cyl.')

%% Plots
figure(2)
plot(disk_perp(:,1),(disk_perp(:,2)-disk_perp(:,3)).*100./disk_perp(:,3),'r')
hold on
plot(sphere(:,1),(sphere(:,2)-sphere(:,3)).*100./sphere(:,3),'b')
plot(sharp_cone_5_deg(:,1),(sharp_cone_5_deg(:,2)-sharp_cone_5_deg(:,3)).*100./sharp_cone_5_deg(:,3),'m')
plot(circular_cylinder(:,1),(circular_cylinder(:,2)-circular_cylinder(:,3)).*100./circular_cylinder(:,3),'g')
hold off
xlabel('Mach number, nd')
ylabel('% Difference from Amalytical')
grid on
legend('Disk (perp.)','Sphere','Sharp-Cone (5 deg.)','Circular Cyl.')