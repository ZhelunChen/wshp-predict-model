function cap = predict_hp_cool_capacity(T_loa_in, T_sou_in, y_speed)
%     q_cool_nom=7.210; % Nominal Cooling Capacity (kW)
%     coeff= [1.4560, -0.0434, 0.0016, -0.0046, 0.0007, -0.0016];
%     cap_wo_speed = biquadratic(T_loa_in, T_sou_in, coeff)*q_cool_nom; % Capacity w/o Speed
%     cap=y_speed.*cap_wo_speed; % Capacity scaled by speed
    
    % retrain by Yicheng (2024/02/29)
    % Model Parameters
    coeff = [1.62445116768452;0.542640588815011;9.23751454293513;-0.0244700650250555;-0.0174478084335477;0.00907381862714126;-0.00159204216924559;-0.00433225327038308;-4.62774065168163;-36.4388897427878];
    % Features
    x1 = T_loa_in;
    x2 = T_sou_in;
    x3 = y_speed;

    % Use matrix multiplication for applying coefficients
    cap = coeff(1)*x1.*x3+coeff(2)*x2.*x3+coeff(3)*x3.*x3+...
        coeff(4)*x1.*x1.*x3+coeff(5)*x1.*x2.*x3+coeff(6)*x1.*x3.*x3+...
        coeff(7)*x2.*x2.*x3+coeff(8)*x2.*x3.*x3+coeff(9)*x3.*x3.*x3+coeff(10)*1.*x3;
end