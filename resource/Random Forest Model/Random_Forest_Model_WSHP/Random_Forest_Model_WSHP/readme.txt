%{ 
How to use predicition model:

reg1 = random forest model

y=speed_ratio output from MPC 
Tz=zone_temperature output from MPC
Tz_set=predicted_zone_temperature_setpoint
%}


Tz_set=predict(reg1,[y,Tz]) % Can input y & Tz as scalars or column vectors