function output = SpaceLinerEndpoint(input)
auxdata = input.auxdata;

% vF = input.phase(8).finalstate(4);
% 
% output.objective = -vF;




% if isempty(input.phase(1).integral)
%     input.phase(1).integral = 0;
% end
% 
% if isempty(input.phase(2).integral)
%     input.phase(2).integral = 0;
% end
% 
% if isempty(input.phase(3).integral)
%     input.phase(3).integral = 0;
% end
% 
% if isempty(input.phase(4).integral)
%     input.phase(4).integral = 0;
% end
% 
% if isempty(input.phase(5).integral)
%     input.phase(5).integral = 0;
% end
% 
% if isempty(input.phase(6).integral)
%     input.phase(6).integral = 0;
% end
% 
% if isempty(input.phase(7).integral)
%     input.phase(7).integral = 0;
% end
% 
% if isempty(input.phase(8).integral)
%     input.phase(8).integral = 0;
% end


% input
% input.phase

output.objective = input.phase(1).integral+...
    input.phase(2).integral+...
    input.phase(3).integral+...
    input.phase(4).integral+...
    input.phase(5).integral+...
    input.phase(6).integral+...
    input.phase(7).integral+...
    input.phase(8).integral;



% gammaF = input.phase(1).finalstate(5);
% 
% output.objective = gammaF^2*1000;

%%
t01 = input.phase(1).initialtime;
tf1 = input.phase(1).finaltime;
x01 = input.phase(1).initialstate;
xf1 = input.phase(1).finalstate;

t02 = input.phase(2).initialtime;
tf2 = input.phase(2).finaltime;
x02 = input.phase(2).initialstate;
xf2 = input.phase(2).finalstate;

t03 = input.phase(3).initialtime;
tf3 = input.phase(3).finaltime;
x03 = input.phase(3).initialstate;
xf3 = input.phase(3).finalstate;

t04 = input.phase(4).initialtime;
tf4 = input.phase(4).finaltime;
x04 = input.phase(4).initialstate;
xf4 = input.phase(4).finalstate;

t05 = input.phase(5).initialtime;
tf5 = input.phase(5).finaltime;
x05 = input.phase(5).initialstate;
xf5 = input.phase(5).finalstate;

t06 = input.phase(6).initialtime;
tf6 = input.phase(6).finaltime;
x06 = input.phase(6).initialstate;
xf6 = input.phase(6).finalstate;

t07 = input.phase(7).initialtime;
tf7 = input.phase(7).finaltime;
x07 = input.phase(7).initialstate;
xf7 = input.phase(7).finalstate;

t08 = input.phase(8).initialtime;
tf8 = input.phase(8).finaltime;
x08 = input.phase(8).initialstate;
xf8 = input.phase(8).finalstate;

output.eventgroup(1).event = [x02-xf1 t02-tf1];
output.eventgroup(2).event = [x03-xf2 t03-tf2];
output.eventgroup(3).event = [x04-xf3 t04-tf3];
output.eventgroup(4).event = [x05-xf4 t05-tf4];
output.eventgroup(5).event = [x06-xf5 t06-tf5];
output.eventgroup(6).event = [x07-xf6 t07-tf6];
output.eventgroup(7).event = [x08-xf7 t08-tf7];

output.eventgroup(8).event = [tf1-t01,tf2-t02,tf3-t03,tf4-t04,tf5-t05,tf6-t06,tf7-t07,tf8-t08];
end
