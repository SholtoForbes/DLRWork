function ydot = VehicleModel_forward(f_t, f_y,auxdata,Alpha,throttle,timeF)

states1.alt = f_y(1);
states1.gamma = f_y(2);
states1.v = f_y(3);
states1.zeta = f_y(4);
states1.lat = f_y(5);
states1.lon = f_y(6);
states1.mFuel = f_y(7);

controls1.Alpha = Alpha;
% controls1.eta = eta;
%

[altdot,londot,latdot,gammadot,a,zetadot, q, M, Fd, rho,L,Fueldt] = SpaceLinerVehicleModel(f_t,states1,controls1,throttle,auxdata,timeF);

ydot = [altdot;gammadot;a;zetadot;latdot;londot;-Fueldt];
% =========================================================================
end








