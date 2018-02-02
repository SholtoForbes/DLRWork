function output = SpaceLinerEndpoint(input)
auxdata = input.auxdata;

vF = input.phase(1).finalstate(4);

output.objective = -vF;



%%
% t01 = input.phase(1).initialtime;
% tf1 = input.phase(1).finaltime;
% x01 = input.phase(1).initialstate;
% xf1 = input.phase(1).finalstate;
% 
% t02 = input.phase(2).initialtime;
% tf2 = input.phase(2).finaltime;
% x02 = input.phase(2).initialstate;
% xf2 = input.phase(2).finalstate;
% 
% lon2f  = input.phase(2).finalstate(2);
% lat2f  = input.phase(2).finalstate(3);

% output.eventgroup(1).event = [x02(1:9)-xf1(1:9) t02-tf1];
% output.eventgroup(2).event = [lon2f-auxdata.target_lon lat2f-auxdata.target_lat];
end
