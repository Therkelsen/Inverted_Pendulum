% Plot af teoretiske gains
pos1 = table2array(readtable("Excel data.xlsx", 'Range','A2:A56'));
ang1 = table2array(readtable("Excel data.xlsx", 'Range','B2:B56'));
x1 = (10:10:550)
figure 
hold on
plot(x1,ang1)
plot(x1,pos1)
legend("Pendulum angle [rad]","Cart position [m]")

x = (10:10:10000);
% Plot af teoretiske pendul gains
pos2 = table2array(readtable("Excel data.xlsx", 'Range','D2:D1001'));
ang2 = table2array(readtable("Excel data.xlsx", 'Range','E2:E1001'));
figure
hold on
plot(x,ang2)
plot(x,pos2)
legend("Pendulum angle [rad]","Cart position [m]")


% Plot af endelige gains
pos3 = table2array(readtable("Excel data.xlsx", 'Range','G2:G1001'));
ang3 = table2array(readtable("Excel data.xlsx", 'Range','H2:H1001'));
figure
hold on
plot(x,ang3)
plot(x,pos3)
legend("Pendulum angle [rad]","Cart position [m]")


% PLot af endelige gains med halv vognlængdes støj
pos4 = table2array(readtable("Excel data.xlsx", 'Range','J2:J1001'));
ang4 = table2array(readtable("Excel data.xlsx", 'Range','K2:K1001'));
figure
hold on
plot(x,ang4)
plot(x,pos4)
legend("Pendulum angle [rad]","Cart position [m]")

% 
% figure
% hold on
% plot(x,pos3)
% plot(x,pos4)
% 
% figure
% hold on
% plot(x,ang3)
% plot(x,ang4)

