clear all

SLEG = dlmread('SLEG.txt');
TOSCA = dlmread('TOSCA');

HeatMin1 = dlmread('HeatMin1');
HeatMin3 = dlmread('HeatMin3');


t_SLEG = SLEG(:,1);
t_TOSCA = TOSCA(:,1);
t_HeatMin1 = HeatMin1(:,1);
t_HeatMin3 = HeatMin3(:,1);


alt_SLEG = SLEG(:,11);
alt_TOSCA = TOSCA(:,5);
alt_HeatMin1 = HeatMin1(:,5);
alt_HeatMin3 = HeatMin3(:,5);

figure()
hold on
% plot(t_SLEG,alt_SLEG/1000);
plot(t_TOSCA,alt_TOSCA);
plot(t_HeatMin1,alt_HeatMin1/1000);
plot(t_HeatMin3,alt_HeatMin3/1000);


% heat_SLEG = SLEG(:,);
heat_TOSCA = TOSCA(:,14);
heat_HeatMin1 = HeatMin1(:,9);
heat_HeatMin3 = HeatMin3(:,9);


figure()
hold on
% plot(t_SLEG,alt_SLEG/1000);
plot(t_TOSCA,heat_TOSCA);
plot(t_HeatMin1,heat_HeatMin1/1000000);
plot(t_HeatMin3,heat_HeatMin3/1000000);
