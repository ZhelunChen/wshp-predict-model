%{
% High-level MPC for TAMU WSHP System
% author: Caleb Calfa
% email: cjcalfa@tamu.edu
% Revisions:
    Rev1 - Added Zone Temperature Upper & Lower Bounds to Public Properties
%}

classdef WSHP_MPC_rev2< matlab.System
    properties
        PH = 1 % scalar: prediction horizon
        CH = 1 % scalar: control horizon
        dt = 900 % scalar: time step

        occ_start = 5 % scalar: hour index of occupancy start time
        occ_end = 21 % scalar: hour index of occupany end time

        w = [1, 0, 0] % vector: weights in the objective function
        u_lb = [0, 0, 0, 0]
        u_ub = [1, 0.1, 20, 5]

        Tz_lb_unocc=12;
        Tz_ub_unocc=30;
        Tz_lb_occ=22;
        Tz_ub_occ=26;

        t_start= 237*24*3600

        Location=0 % 0 = Atlanta, Buffalo, Tucson; 1 = New York

        price_tou = [0.02987, 0.02987, 0.02987, 0.02987,...
            0.02987, 0.02987, 0.04667, 0.04667,...
            0.04667, 0.04667, 0.04667, 0.04667,...
            0.15877, 0.15877, 0.15877, 0.15877,...
            0.15877, 0.15877, 0.15877, 0.04667,...
            0.04667, 0.04667, 0.02987, 0.02987]

        Tb=17.8
        dep=112.8
    end

    % Hidden Properties
    properties (Access = private)        
        Rg=6.693542227015886
        Re=2.3677482157900611
        Ri=1.6861300148536231
        Rw=11.106467926191005
        Cwe=50000
        Cwi=37515.136743737348
        Cai=4881.0882226682534

        nu = 4 % input vector size: [y, eps, q, p]
       
        mwi=0.3155
        di=0.03 
        da=0.033 
        s=0.05
        
        db=0.2 
        k_w=0.611 
        k_p=0.5
        k_g=2.8

        Cp_g=800 
        rho_g=1600 
        rho_w=996 
        Cp_w=4.18e3 

        mew_w=8.67*10^-4 

        Nu=3.66

        timestep
    end

    methods (Static)
        % WSHP Cooling Capacity
        function cap = predict_hp_cool_capacity(obj,T_loa_in, T_sou_in, y_speed)
            q_cool_nom=-7.210; % Nominal Cooling Capacity (kW)
            coeff = [1.4560,-0.0434,0.0018,-0.0045,0.0003,-0.0011];
            cap_wo_speed = obj.biquadratic(T_loa_in, T_sou_in, coeff)*q_cool_nom; % Capacity w/o Speed
            cap=y_speed.*cap_wo_speed; % Capacity scaled by speed
        end

        % WSHP Cooling Power Consumption
        function pow = predict_hp_cool_power(obj,T_loa_in, T_sou_in, y_speed)
            p_cool_nom=1.335; % Nominal Cooling Power (kW)
            p_shift=0.12; % Power Offset (kW) -> Due to continuous fan operation
            coeff = [1.4560,-0.0435,-0.0005,-0.0044,0.0004,0.0009];
            pow_wo_speed = obj.biquadratic(T_loa_in, T_sou_in, coeff)*p_cool_nom; % Power w/o Speed
            pow=y_speed.*pow_wo_speed+p_shift; % Power scaled by speed
        end

        % WSHP Waterside dT computation
        function [dTw] = predict_dTw(obj,q_cool)
            coeff = [0,1/0.97751724];
            dTw = obj.linear(q_cool, coeff);
        end

        function res = biquadratic(x1, x2, coeff)
            res = coeff(1) + (coeff(2)+coeff(3)*x1).*x1 + (coeff(4)+coeff(5)*x2).*x2 + coeff(6)*x1.*x2;
        end

        function res = linear(x1, coeff)
            res = coeff(1) + coeff(2)*x1;
        end

        % RC_Network + GLHE
        function [new_state_PH] = RC_GLHE(obj,states,predictor,q_cool_ph)

            import casadi.*

            % predictor: 7-by-PH: [out_temp,qint, qhvac, qsolout, qradin, To, Tb]
            % state: 7-by-1: [Tz, Twin, Twout, Tg1, Tg2, T1, T2]

            % GLHE Resistance & Capacitance Calculations
            Rconv=1/(obj.Nu*obj.k_w*pi())/obj.dep; % Convection Resistance from water to pipe (K/W)
            Rcond_p=log(obj.da/obj.di)/(2*pi()*obj.k_p)/obj.dep; % Conduction Resistance thru pipe (K/W)

            Rg=acosh((obj.db^2+obj.da^2-obj.s^2)/(2*obj.db*obj.da))/(2*pi()*obj.k_g)*(1.601-0.888*obj.s/obj.db)/obj.dep; % Total grout resistance (K/W)
            x=log(sqrt(obj.db^2+2*obj.da^2)/(2*obj.da))/log(obj.db/(sqrt(2)*obj.da))/obj.dep; % Grout resistance separation factor

            Rcond_g=x*Rg;  % Conduction resistance from pipe to grout (K/W)
            Rgb=(1-x)*Rg;  % Conduction resistance from grout to borehole outer diameter (K/W)

            Rfg=Rconv+Rcond_p+Rcond_g; % Conduction resistance from water to grout (K/W)

            Rar=acosh((2*obj.s^2-obj.da^2)/obj.da^2)/(2*pi()*obj.k_g)/obj.dep; % Resistance between pipe outer walls (K/W)

            Rgg=2*Rgb*(Rar-2*x*Rg)/(2*Rgb-Rar+2*x*Rg); % Resistance between two grout zones (K/W)

            % According to Bauer et al. 2011 Rgg can be negative and still satisfy 2nd Law of Thermodynamics as long as:
            %z=(1/Rgg+1/(2*Rgb))^-1; % is greater then zero

            Cg=obj.rho_g*obj.Cp_g*pi()/4*(obj.db^2/2-obj.da^2)*obj.dep; % Capacitance of one grout zone
            Cw=obj.rho_w*obj.Cp_w*pi()/4*(obj.di^2)*obj.dep; % heat capacitance of u-tube

            % Calculate Delta Tw Array
            dTw=obj.predict_dTw(obj,q_cool_ph); % Delta Waterside Temperature

            % Matrix Creation
            A = zeros(7,7);
            B = zeros(7,7);

            % RC Coefficients
            A(1,1) = -1/obj.Cai*(1/obj.Rg+1/obj.Ri);
            A(1,3) = 1/(obj.Cai*obj.Ri);
            A(2,2) = -1/obj.Cwe*(1/obj.Re+1/obj.Rw);
            A(2,3) = 1/(obj.Cwe*obj.Rw);
            A(3,1) = 1/(obj.Cwi*obj.Ri);
            A(3,2) = 1/(obj.Cwi*obj.Rw);
            A(3,3) = -1/obj.Cwi*(1/obj.Rw+1/obj.Ri);

            B(1,1) = 1/(obj.Cai*obj.Rg);
            B(1,2) = 1/obj.Cai;
            B(1,3) = 1/obj.Cai;
            B(2,1) = 1/(obj.Cwe*obj.Re);
            B(2,4) = 1/obj.Cwe;
            B(3,5) = 1/obj.Cwi;

            % GLHE Coeffients
            A(4,4) = -1/Cg*(1/Rfg+1/Rgg+1/Rgb);
            A(4,5) = 1/(Cg*Rgg);
            A(4,6) = 1/(Cg*Rfg);
            A(5,4) = 1/(Cg*Rgg);
            A(5,5) = -1/Cg*(1/Rfg+1/Rgg+1/Rgb);
            A(5,7) = 1/(Cg*Rfg);
            A(6,4) = 1/(Cw*Rfg);
            A(6,6) = -1/Cw*(obj.mwi*obj.Cp_w+1/Rfg);
            A(7,5) = 1/(Cw*Rfg);
            A(7,6) = 1/Cw*(obj.mwi*obj.Cp_w);
            A(7,7) = -1/Cw*(obj.mwi*obj.Cp_w+1/Rfg);

            B(4,7) = 1/(Cg*Rgb);
            B(5,7) = 1/(Cg*Rgb);
            B(6,6) = 1/Cw*(obj.mwi*obj.Cp_w);

            % Container Intialization
            new_state_PH = [];

            % Symbolic Solving Framework
            S = MX.sym('s', 7,1); % [Tz, Two, Twi, Tg1, Tg2, T1, T2]
            D = MX.sym('d', 7,1); % [out_temp,qint, qhvac, qsolout, qradin, To, Tb]
            Sdot = A*S + B*D; 
            ode = struct('x', S, 'p', D, 'ode', Sdot);% 'quad', C*S)
            opts = struct('tf', obj.dt);
            F = integrator('F', 'cvodes', ode, opts);

            % Initial State
            x0=states;
            for  i = 1:obj.PH % Solve for zone temperature using prediction inputs for entire PH
                % Calculate Outlet Temperature
                % Initialize Disturbance
                disturb_i = MX(predictor(1:5,i)); 
                disturb_i(3)=-q_cool_ph(i); % Update Disturbance Matrix (Negative Sign Cooling Load)
                % Calculate Outlet Temperature
                Tw_o=dTw(i)+x0(7);
                % Append Tw_o and Tb
                disturb_i=[disturb_i;Tw_o;obj.Tb];
                % Evaluate F
                Fk = F('x0', x0, 'p', disturb_i);
                x0 = Fk.xf;
                % save predictions
                new_state_PH = [new_state_PH,x0];
            end
        end

        % RC_Network + GLHE Numeric Vector
        function [new_state_PH] = RC_GLHE_full(obj,states,predictor,q_cool_ph)

            import casadi.*

            % predictor: 7-by-PH: [out_temp,qint, qhvac, qsolout, qradin, To, Tb]
            % state: 7-by-1: [Tz, Twin, Twout, Tg1, Tg2, T1, T2]

            % GLHE Resistance & Capacitance Calculations
            Rconv=1/(obj.Nu*obj.k_w*pi())/obj.dep; % Convection Resistance from water to pipe (K/W)
            Rcond_p=log(obj.da/obj.di)/(2*pi()*obj.k_p)/obj.dep; % Conduction Resistance thru pipe (K/W)

            Rg=acosh((obj.db^2+obj.da^2-obj.s^2)/(2*obj.db*obj.da))/(2*pi()*obj.k_g)*(1.601-0.888*obj.s/obj.db)/obj.dep; % Total grout resistance (K/W)
            x=log(sqrt(obj.db^2+2*obj.da^2)/(2*obj.da))/log(obj.db/(sqrt(2)*obj.da))/obj.dep; % Grout resistance separation factor

            Rcond_g=x*Rg;  % Conduction resistance from pipe to grout (K/W)
            Rgb=(1-x)*Rg;  % Conduction resistance from grout to borehole outer diameter (K/W)

            Rfg=Rconv+Rcond_p+Rcond_g; % Conduction resistance from water to grout (K/W)

            Rar=acosh((2*obj.s^2-obj.da^2)/obj.da^2)/(2*pi()*obj.k_g)/obj.dep; % Resistance between pipe outer walls (K/W)

            Rgg=2*Rgb*(Rar-2*x*Rg)/(2*Rgb-Rar+2*x*Rg); % Resistance between two grout zones (K/W)

            % According to Bauer et al. 2011 Rgg can be negative and still satisfy 2nd Law of Thermodynamics as long as:
            %z=(1/Rgg+1/(2*Rgb))^-1; % is greater then zero

            Cg=obj.rho_g*obj.Cp_g*pi()/4*(obj.db^2/2-obj.da^2)*obj.dep; % Capacitance of one grout zone
            Cw=obj.rho_w*obj.Cp_w*pi()/4*(obj.di^2)*obj.dep; % heat capacitance of u-tube

            % Calculate Delta Tw Array
            dTw=obj.predict_dTw(obj,q_cool_ph); % Delta Waterside Temperature

            % Matrix Creation
            A = zeros(7,7);
            B = zeros(7,7);

            % RC Coefficients
            A(1,1) = -1/obj.Cai*(1/obj.Rg+1/obj.Ri);
            A(1,3) = 1/(obj.Cai*obj.Ri);
            A(2,2) = -1/obj.Cwe*(1/obj.Re+1/obj.Rw);
            A(2,3) = 1/(obj.Cwe*obj.Rw);
            A(3,1) = 1/(obj.Cwi*obj.Ri);
            A(3,2) = 1/(obj.Cwi*obj.Rw);
            A(3,3) = -1/obj.Cwi*(1/obj.Rw+1/obj.Ri);

            B(1,1) = 1/(obj.Cai*obj.Rg);
            B(1,2) = 1/obj.Cai;
            B(1,3) = 1/obj.Cai;
            B(2,1) = 1/(obj.Cwe*obj.Re);
            B(2,4) = 1/obj.Cwe;
            B(3,5) = 1/obj.Cwi;

            % GLHE Coeffients
            A(4,4) = -1/Cg*(1/Rfg+1/Rgg+1/Rgb);
            A(4,5) = 1/(Cg*Rgg);
            A(4,6) = 1/(Cg*Rfg);
            A(5,4) = 1/(Cg*Rgg);
            A(5,5) = -1/Cg*(1/Rfg+1/Rgg+1/Rgb);
            A(5,7) = 1/(Cg*Rfg);
            A(6,4) = 1/(Cw*Rfg);
            A(6,6) = -1/Cw*(obj.mwi*obj.Cp_w+1/Rfg);
            A(7,5) = 1/(Cw*Rfg);
            A(7,6) = 1/Cw*(obj.mwi*obj.Cp_w);
            A(7,7) = -1/Cw*(obj.mwi*obj.Cp_w+1/Rfg);

            B(4,7) = 1/(Cg*Rgb);
            B(5,7) = 1/(Cg*Rgb);
            B(6,6) = 1/Cw*(obj.mwi*obj.Cp_w);

            % Container Intialization
            new_state_PH = [];
           
            % Symbolic Solving Framework
            S = MX.sym('s', 7,1); % [Tz, Two, Twi, Tg1, Tg2, T1, T2]
            D = MX.sym('d', 7,1); % [out_temp,qint, qhvac, qsolout, qradin, To, Tb]
            Sdot = A*S + B*D; 
            ode = struct('x', S, 'p', D, 'ode', Sdot);% 'quad', C*S)
            opts = struct('tf', obj.dt);
            F = integrator('F', 'cvodes', ode, opts);

            % Initial State
            x0=states;
            for  i = 1:obj.PH % Solve for zone temperature using prediction inputs for entire PH
                % Calculate Outlet Temperature
                % Initialize Disturbance
                disturb_i = MX(predictor(1:5,i)); 
                disturb_i(3)=-q_cool_ph(i); % Update Disturbance Matrix (Negative Sign Cooling Load)
                % Calculate Outlet Temperature
                Tw_o=dTw(i)+x0(7);
                % Append Tw_o and Tb
                disturb_i=[disturb_i;Tw_o;obj.Tb];
                % Evaluate F
                Fk = F('x0', x0, 'p', disturb_i);
                x0 = Fk.xf;
                x0=full(evalf(x0)); % Convert to dense array
                % save predictions
                new_state_PH = [new_state_PH,x0];
            end
        end
    end

    methods (Access = protected)
        % Block Interface
        % specify type, size, complex .. for multiple inputs.
        function num = getNumInputsImpl(~)
            % number of block inputs:
            % prediction: 5-by-PH 
            % state: 7-by-1: [Tz, Twe, Twi, Tg1, Tg2, T1, T2]
            % time: 1-by-1
            % u_prev: 1-by-obj.nu*obj.PH
            num = 4;
        end

        function num = getNumOutputsImpl(~)
            num = 7; % [Tz_opt,y_opt,e_opt,Tub_i,Tlb_i,status,x_opt]
        end
        function varargout = getOutputDataTypeImpl(obj)
            % multiple outputs
            varargout = cell(1, getNumOutputs(obj));
            for i = 1:getNumOutputs(obj)
                varargout{i} = 'double';
            end
        end
        function varargout = getInputDataTypeImpl(~)
            varargout = cell(1, getNumInputs(obj));
            for i = 1:getNumInputs(obj)
                varargout{i} = 'double';
            end
        end
        function varargout = getOutputSizeImpl(obj)
            % vector
            % sz1 = [1,obj.nu];
            varargout = cell(1, getNumOutputs(obj));
            for i = 1:getNumOutputs(obj)
                varargout{i} = [1, 1];
            end
            varargout{end}=[1,obj.nu*obj.PH]; % xopt size
        end

        function [sz1, sz2, sz3, sz4] = getInputSizeImpl(~)
            % predictor: 5-by-PH: [out_temp,qint, qhvac, qsolout, qradin]
            % state: 7-by-1: [Tz, Twe, Twi, Tg1, Tg2, T1, T2
            % time 1-by-1
            % u_prev 1-by-obj.nu*obj*PH
            sz1=[5, obj.PH];
            sz2 = [7, 1];
            sz3 = [1, 1];
            sz4= [1,obj.nu*obj.PH];
        end

        function varargout = isInputComplexImpl(~)
            varargout = cell(1, getNumInputs(obj));
            for i = 1:getNumInputs(obj)
                varargout{i} = false;
            end
        end
        function varargout = isOutputComplexImpl(obj)
            varargout = cell(1, getNumOutputs(obj));
            for i = 1:getNumOutputs(obj)
                varargout{i} = false;
            end
        end
        function varargout = isInputFixedSizeImpl(~)
            % need to guard fixed size for each input
            varargout = cell(1, getNumInputs(obj));
            for i = 1:getNumInputs(obj)
                varargout{i} = true;
            end
        end
        function varargout = isOutputFixedSizeImpl(obj)
            % need to guard fixed size for each output
            varargout = cell(1, getNumOutputs(obj));
            for i = 1:getNumOutputs(obj)
                varargout{i} = true;
            end
        end

        % Set Up Optimization
        function setupImpl(obj, ~, ~)
            obj.timestep=0;
        end

        % Overarching Function (Simulink Block Interface in final script)
        function [Tz_opt,y_opt,e_opt,Tub_i,Tlb_i,status,x_opt] = stepImpl(obj,predictors,state,time,u_prev) 

            import casadi.*

            % Define Temperature Upper & Lower Bounds
            Tub = obj.Tz_ub_unocc*ones(1,24);
            Tlb = obj.Tz_lb_unocc*ones(1,24);
            Tub(obj.occ_start:obj.occ_end) = obj.Tz_ub_occ;
            Tlb(obj.occ_start:obj.occ_end) = obj.Tz_lb_occ;

            % Get Price and Temp Setpoint lower & upper bounds over PH
            price_ph = zeros(1, obj.PH);
            Tub_ph = zeros(1, obj.PH);
            Tlb_ph = zeros(1, obj.PH);
            for i = 1:obj.PH
                t = time + (i-1)*obj.dt;
                h = floor(mod(t, 86400)/3600) + 1;
                price_ph(i) = obj.price_tou(h); % TOU price
                Tub_ph(i)=Tub(h); % Temperature Setpoint Upper Bounds
                Tlb_ph(i)=Tlb(h); % Temperature Setpoint Lower Bounds
            end

            Tub_i=Tub_ph(1,1);
            Tlb_i=Tlb_ph(1,1);

            % Variables: - use single quote. double quote failed here.
            nvar=obj.nu; % [y,e,q,p]
            nvars=obj.nu*obj.PH; % Number of Variables

            u = MX.sym('u', 1, nvars); % Creation of symbolic matrix in casadi
            y=u(1:nvar:nvars); % Compressor Speed
            e=u(2:nvar:nvars); % Slack Variable
            q=u(3:nvar:nvars); % WSHP Capacity
            p=u(4:nvar:nvars); % WSHP Power

            % ----------------------------
            % objective function
            % ----------------------------
            % Control Timestep Changes
            du=(u-u_prev).^2; 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Yicheng/ Aaron Add Test Case Information Here
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if obj.Location>0
                % calculate the billable demand cost for NewYork case (YC 2023/1030)
                AllPH_temp=[1:obj.PH];
                DailyPeak=[(12*(3600/obj.dt)+1) : (19*(3600/obj.dt))];
                DailyOffpeak=AllPH_temp(~ismember(AllPH_temp,DailyPeak));

                % determine the peak or offpeak period in the ph
                dT_index=floor(mod(t, 86400)/(obj.dt)) + 1;
                PeakIndex=DailyPeak-dT_index;
                for k=1:size(PeakIndex,2)
                    if PeakIndex(k)<1
                        PeakIndex(k)=PeakIndex(k)+96;
                    end
                end
                PeakIndex=sort(PeakIndex);
                OffpeakIndex=DailyOffpeak-dT_index;
                for k=1:size(OffpeakIndex,2)
                    if OffpeakIndex(k)<1
                        OffpeakIndex(k)=OffpeakIndex(k)+96;
                    end
                end
                OffpeakIndex=sort(OffpeakIndex);

                % calculate the rolling average

                % for peak period
                % check discontinuities
                discont_indices = find(diff(PeakIndex) > 1);
                PeakSepa_vectors = cell(length(discont_indices) + 1, 1);
                % Split the vector at discontinuities
                start_index = 1;
                for i = 1:length(discont_indices)
                    end_index = discont_indices(i);
                    PeakSepa_vectors{i} = PeakIndex(start_index:end_index);
                    start_index = end_index + 1;
                end
                PeakSepa_vectors{end} = PeakIndex(start_index:end);
                % calculate the rolling average power
                for s=1:size(PeakSepa_vectors,1)
                    if size(PeakSepa_vectors{s},2)<4
                        PeakMovingAvePower_temp{s,1}=mean(p(PeakSepa_vectors{s,1}));
                    else
                        for k=1:(size(PeakSepa_vectors{s},2)-3)
                            PeakMovingAvePower_temp{s,k}=(p(PeakSepa_vectors{s,1}(k))+p(PeakSepa_vectors{s,1}(k+1))+p(PeakSepa_vectors{s,1}(k+2))+p(PeakSepa_vectors{s,1}(k+3)))/4;
                        end
                    end
                end
                PeakMovingAvePower=[];
                for r=1:size(PeakMovingAvePower_temp,1)
                    for c=1:size([PeakMovingAvePower_temp{r,:}],2)
                        PeakMovingAvePower=[PeakMovingAvePower PeakMovingAvePower_temp{r,c}];
                    end
                end

                % for offpeak period
                % check discontinuities
                discont_indices = find(diff(OffpeakIndex) > 1);
                OffpeakSepa_vectors = cell(length(discont_indices) + 1, 1);
                % Split the vector at discontinuities
                start_index = 1;
                for i = 1:length(discont_indices)
                    end_index = discont_indices(i);
                    OffpeakSepa_vectors{i} = OffpeakIndex(start_index:end_index);
                    start_index = end_index + 1;
                end
                OffpeakSepa_vectors{end} = OffpeakIndex(start_index:end);
                % calculate the rolling average power
                for s=1:size(OffpeakSepa_vectors,1)
                    if size(OffpeakSepa_vectors{s},2)<4
                        OffpeakMovingAvePower_temp{s,1}=mean(p(OffpeakSepa_vectors{s,1}));
                    else
                        for k=1:(size(OffpeakSepa_vectors{s},2)-3)
                            OffpeakMovingAvePower_temp{s,k}=(p(OffpeakSepa_vectors{s,1}(k))+p(OffpeakSepa_vectors{s,1}(k+1))+p(OffpeakSepa_vectors{s,1}(k+2))+p(OffpeakSepa_vectors{s,1}(k+3)))/4;
                        end
                    end
                end
                OffpeakMovingAvePower=[];
                for r=1:size(OffpeakMovingAvePower_temp,1)
                    for c=1:size([OffpeakMovingAvePower_temp{r,:}],2)
                        OffpeakMovingAvePower=[OffpeakMovingAvePower OffpeakMovingAvePower_temp{r,c}];
                    end
                end

                % calculate the maximum rolling average power
                alp=20;  % the coefficient of approximation
                MaxOffpeak_temp=0;
                MaxPeak_temp=0;
                for i_MovPeak=1:size(PeakMovingAvePower,2)
                    MaxPeak_temp=MaxPeak_temp+exp(alp*PeakMovingAvePower{i_MovPeak});
                end
                MaxPeakMovingAvePower=1/alp*log(MaxPeak_temp);

                for i_MovOffpeak=1:size(OffpeakMovingAvePower,2)
                    MaxOffpeak_temp=MaxOffpeak_temp+exp(alp*OffpeakMovingAvePower{i_MovOffpeak});
                end
                MaxOffpeakMovingAvePower=1/alp*log(MaxOffpeak_temp);

                Cost=MaxPeakMovingAvePower*21.54/30+MaxOffpeakMovingAvePower*8.37/30;
                J_1=sum((Cost.*obj.dt/3600));
            else
                J_1=sum((price_ph.*p*obj.dt/3600));
            end

            % Calculate Objective Function
            J_2=sum(e.^2);
            J_3=sum(du);
            fval_sum=obj.w(1)*J_1+obj.w(2)*J_2+obj.w(3)*J_3;
            fv = Function('fval', {u}, {fval_sum}); % Function definition in casidi u-> Input, fval_sum-> Output
            f = fv(u);  

            % ---------------------------------
            % Constraints
            % ---------------------------------
            C_Tz=zeros(1,7);
            C_Tz(1)=1;

            C_Tw=zeros(1,7);
            C_Tw(7)=1;

            [new_state_PH] = obj.RC_GLHE(obj,state,predictors,q);
            Tz_pred_PH=C_Tz*new_state_PH;
            Tw_pred_PH=C_Tw*new_state_PH;
            q_pred_PH=obj.predict_hp_cool_capacity(obj,Tz_pred_PH,Tw_pred_PH,y); % Ensure signs are correct!!!
            p_pred_PH=obj.predict_hp_cool_power(obj,Tz_pred_PH,Tw_pred_PH,y);
 
            % Define Inequalities, Equalities, & Bounds
            g = {};
            lbg = []; % Inequality Lower Bound
            ubg = []; % Inequality Upper Bound
            lbu = []; % Optimization Vector Lower Bound
            ubu = []; % Optimization Vector Upper Bound

            for k=1:obj.PH
                % Zone Temperature & Slack Variable Inequality
                g = {g{:}, Tz_pred_PH(k)+e(k), Tz_pred_PH(k)-e(k)}; 
                lbg = [lbg; Tlb_ph(k); 0];
                ubg = [ubg; inf; Tub_ph(k)];

                lbu=[lbu;obj.u_lb'];
                ubu=[ubu;obj.u_ub'];
            end
           
            % HVAC Load Equality
            g = {g{:}, ((q+q_pred_PH).^2)'}; 
            lbg = [lbg; zeros(obj.PH,1)];
            ubg = [ubg; zeros(obj.PH,1)];

            % HVAC Power Equality
            g = {g{:}, ((p-p_pred_PH).^2)'}; 
            lbg = [lbg; zeros(obj.PH,1)];
            ubg = [ubg; zeros(obj.PH,1)];

            % Define Optimization in Casadi
            prob = struct('f',f, 'x',u, 'g',vertcat(g{:}));
            options = struct;
            options.print_time = true;
            options.error_on_fail = false;
            opt_solver='ipopt';

            options.ipopt.tol=1e-8;
            options.ipopt.acceptable_tol=1e-7;
            options.ipopt.constr_viol_tol=1e-8;
            options.ipopt.linear_solver='mumps'; % Only Package currently available
            options.ipopt.max_iter=1e7; 
            options.ipopt.max_cpu_time=obj.dt;
            options.ipopt.bound_relax_factor=1e-4; 
            options.ipopt.warm_start_init_point='yes';
            solver = nlpsol('solver', opt_solver, prob, options);
            
            % Optimization Initialization
            u_start=(obj.u_ub+obj.u_lb)./2;
            u_in=repmat(u_start,1,obj.PH);
           
            % Solve Optimization Problem
            sol = solver('x0', u_in, 'lbx', lbu, 'ubx', ubu, 'lbg', lbg, 'ubg', ubg);
            
            % Record Optimization Feasibility
            stats = solver.stats();
            message=stats.return_status % Optimization Feasibility

            if message == "Solve_Succeeded"
                int=0;
            elseif message == "Solved_To_Acceptable_Level"
                int=0;
            elseif message == "Infeasible_Problem_Detected"
                int=2;
            elseif message == "Maximum_Iterations_Exceeded"
                int=3;
            elseif message == "Restoration_Failed"
                int=4;
            else
                int=5;
            end
            
            status=int;

            % Record Full Optimization Solution
            x = full(sol.x);
            f = full(sol.f);
            
            % Next Step Action in case of Infeasibility
            if status>0
                x = [u_prev(5:obj.nu*obj.PH),u_start];
            end
            
            % Optimal Outputs
            y_opt=x(1,1); % Optimal Compressor Speed Ratio
            e_opt=x(1,2); % Optimal Slack Variable
            x_opt=x; % Optimal Solution Vector
            f_opt=f; % Objective Function Value

            q_opt_ph=x(3:nvar:nvars);
            [new_state_PH] = obj.RC_GLHE_full(obj,state,predictors,q_opt_ph);
            Tz_opt=new_state_PH(1,1); 
        end

        function resetImpl(obj)
            % Initialize discrete-state properties.
        end

    end
end