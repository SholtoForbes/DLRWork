function ydot = VehicleModel_forward(f_t, f_y,auxdata,Alpha,throttle,stage)

phase.state(:,1) = f_y(1);
phase.state(:,2) = f_y(6);
phase.state(:,3) = f_y(5);
phase.state(:,4) = f_y(3);
phase.state(:,5) = f_y(2);
phase.state(:,6) = f_y(4);
phase.state(:,7) = f_y(7);
phase.state(:,8) = Alpha;


[altdot,londot,latdot,gammadot,a,zetadot, q, M, Fd, rho,L,Fueldt] = SpaceLinerVehicleModel(f_t,phase,throttle,auxdata,stage);

ydot = [altdot;gammadot;a;zetadot;latdot;londot;-Fueldt];
% =========================================================================
end








