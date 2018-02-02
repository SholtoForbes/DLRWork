function phaseout = SpaceLinerContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics
states1.alt     = input.phase(1).state(:,1);
states1.lon     = input.phase(1).state(:,2);
states1.lat     = input.phase(1).state(:,3);
states1.v       = input.phase(1).state(:,4);
states1.gamma   = input.phase(1).state(:,5);
states1.zeta    = input.phase(1).state(:,6);
states1.mFuel   = input.phase(1).state(:,7);

controls1.Alpha = input.phase(1).state(:,8); % Note the 'controls' class here is not the same as the 'control' defined for GPOPS (these are the vehicle controls, not the control theory controls)
% controls1.eta   = input.phase(1).state(:,9);

controls1.throttle  = 1;

Alphadot1  = input.phase(1).control(:,1);
% etadot1 = input.phase(1).control(:,2);

time1 = input.phase(1).time;

auxdata = input.auxdata;

[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T,Isp1] = SpaceLinerVehicleModel(time1,states1,controls1,auxdata);

phaseout(1).dynamics  = [altdot1, londot1, latdot1, vdot1, gammadot1, azidot1, -Fueldt1, Alphadot1];
phaseout(1).path = q1;

end

%======================================================