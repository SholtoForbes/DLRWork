function output = SpaceLinerEndpoint(input)

% timeF = input.phase(1).finaltime;
% % 
% output.objective = timeF;

output.objective = input.phase.integral;
end
