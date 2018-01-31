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
controls1.eta   = input.phase(1).state(:,9);

controls1.throttle  = 1;

Alphadot1  = input.phase(1).control(:,1);
etadot1 = input.phase(1).control(:,2);

time1 = input.phase(1).time;

auxdata = input.auxdata;

[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T,Isp1] = VehicleModelCombined(states1,controls1,auxdata,Stage);

phaseout(1).dynamics  = [altdot1, londot1, latdot1, vdot1, gammadot1, azidot1, -Fueldt1, Alphadot1, etadot1];
phaseout(1).path = q1;
%%

states2.alt     = input.phase(2).state(:,1);
states2.lon     = input.phase(2).state(:,2);
states2.lat     = input.phase(2).state(:,3);
states2.v       = input.phase(2).state(:,4);
states2.gamma   = input.phase(2).state(:,5);
states2.azi     = input.phase(2).state(:,6);
states2.mFuel   = input.phase(2).state(:,7);

controls2.Alpha = input.phase(2).state(:,8);
controls2.eta   = input.phase(2).state(:,9);
controls2.throttle  = input.phase(2).state(:,10);

Alphadot2  = input.phase(2).control(:,1);
etadot2 = input.phase(2).control(:,2);
throttledot2 = input.phase(2).control(:,3);

time2 = input.phase(1).time;

[altdot2,londot2,latdot2,gammadot2,vdot2,azidot2, q2, M, Fd, rho,L,Fueldt2,T,Isp2] = VehicleModelCombined(states2,controls2,auxdata,Stage);


phaseout(2).dynamics  = [altdot2, londot2, latdot2, vdot2, gammadot2, azidot2, -Fueldt2, Alphadot2, etadot2, throttledot2];


phaseout(2).path = q2;

end

%======================================================