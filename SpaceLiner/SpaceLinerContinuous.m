function phaseout = SpaceLinerContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics
states.alt     = input.phase(1).state(:,1);
states.lon     = input.phase(1).state(:,2);
states.lat     = input.phase(1).state(:,3);
states.v       = input.phase(1).state(:,4);
states.gamma   = input.phase(1).state(:,5);
states.zeta    = input.phase(1).state(:,6);
states.mFuel   = input.phase(1).state(:,7);

controls.Alpha = input.phase(1).state(:,8); % Note the 'controls' class here is not the same as the 'control' defined for GPOPS (these are the vehicle controls, not the control theory controls)
% controls1.eta   = input.phase(1).state(:,9);

% controls1.throttle  = 1;

Alphadot  = input.phase(1).control(:,1);
% etadot1 = input.phase(1).control(:,2);

time = input.phase(1).time;

auxdata = input.auxdata;


throttle = ones(length(states.alt),1);

throttle(time>time(end)/20*5) = 0.891; %using throttle to turn engines off

throttle(time>time(end)/20*6) = 0.812;

throttle(time>time(end)/20*7) = .7333;

throttle(time>time(end)/20*8) = .6545;

throttle(time>time(end)/20*9) = .5757;

throttle(time>time(end)/20*10) = 0.496;

throttle(time>=time(end)*0.6) = 1; %after separation


[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T,Isp1,Isp2,m,heating_rate] = SpaceLinerVehicleModel(time,states,controls,throttle,auxdata,time(end));

% [altdot1, londot1, latdot1, vdot1, gammadot1, azidot1, -Fueldt1', Alphadot1]

phaseout(1).dynamics  = [altdot1, londot1, latdot1, vdot1, gammadot1, azidot1, -Fueldt1', Alphadot];
phaseout(1).path = q1;
phaseout(1).integrand = heating_rate;
end

%======================================================