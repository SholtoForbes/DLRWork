function ydot = VehicleModel_forward(f_t, f_y,auxdata,Alpha,eta)

states1.alt = f_y(1);
states1.gamma = f_y(2);
states1.v = f_y(3);
states1.zeta = f_y(4);
states1.lat = f_y(5);
states1.lon = f_y(6);

controls1.Alpha = Alpha;
controls1.eta = eta;
%

[altdot,londot,latdot,gammadot,a,zetadot] = SpaceLinerVehicleModel(states1,controls1,auxdata);

ydot = [altdot;gammadot;a;zetadot;latdot;londot];
% =========================================================================
end








