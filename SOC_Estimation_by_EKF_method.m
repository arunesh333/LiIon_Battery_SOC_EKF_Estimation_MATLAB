clc; clear;

% Accessing the Data Points varying in time
load('Datasets\Battery_datapoint.mat') 

% Storing the Values into new Structured called 'Li_ion_Battery' for better
% understanding
Li_ion_Battery.RecordingTime = Measured_values.Time;
Li_ion_Battery.Measured_Voltage = Measured_values.Voltage;
Li_ion_Battery.Measured_Current = Measured_values.Current;
Li_ion_Battery.Measured_Temperature = Measured_values.Battery_Temp_degC;

% Calculating SOC (Measured SOC) by Coulomb Counting Method for refernce
Rated_Capacity = 4.81;
Remaining_Capacity = Rated_Capacity + Measured_values.Ah; 
Li_ion_Battery.Measured_SOC = (Remaining_Capacity / Rated_Capacity) * 100;

% Downsampling the Data points (1 : 5)
Li_ion_Battery.RecordingTime = Li_ion_Battery.RecordingTime(1:5:end);
Li_ion_Battery.Measured_Voltage = Li_ion_Battery.Measured_Voltage(1:5:end);
Li_ion_Battery.Measured_Current = Li_ion_Battery.Measured_Current(1:5:end);
Li_ion_Battery.Measured_Temperature = Li_ion_Battery.Measured_Temperature(1:5:end);
Li_ion_Battery.Measured_SOC = Li_ion_Battery.Measured_SOC(1:5:end);

% Flip the sign of current: In this dataset, discharge current is positive. 
% For modeling, we often assume discharge as negative — hence the sign change.
Li_ion_Battery.Measured_Current_R = -Li_ion_Battery.Measured_Current;

% Convert time from seconds to hours for better readability in plots
Li_ion_Battery.RecordingTime_Hours = Li_ion_Battery.RecordingTime / 3600;

% Calling the EKF Estimation funcition
[Estimated_SOC, Estimated_Vt, Vt_Error] = EKF_function(Li_ion_Battery.Measured_Current_R, Li_ion_Battery.Measured_Voltage, Li_ion_Battery.Measured_Temperature);

% --- Plot Results ---
% 1. Measured Terminal Voltage
figure;
plot(Li_ion_Battery.RecordingTime_Hours, Li_ion_Battery.Measured_Voltage, 'b');
xlabel('Time [hr]');
ylabel('Voltage [V]');
title('Measured Terminal Voltage');
grid on;

% 2. EKF Estimated Terminal Voltage
figure;
plot(Li_ion_Battery.RecordingTime_Hours, Estimated_Vt, 'r');
xlabel('Time [hr]');
ylabel('Voltage [V]');
title('EKF Estimated Terminal Voltage');
grid on;

% 3.Voltage Estimation Error (Measured - Estimated)
figure;
plot(Li_ion_Battery.RecordingTime_Hours, (Li_ion_Battery.Measured_Voltage - Estimated_Vt)*100, 'k');
xlabel('Time [hr]');
ylabel('Voltage Error [V] [in %]');
title('Terminal Voltage Estimation Error');
grid on;
% 4. SOC Comparison
figure;
plot(Li_ion_Battery.RecordingTime_Hours, Li_ion_Battery.Measured_SOC, 'b');
hold on;
plot(Li_ion_Battery.RecordingTime_Hours, Estimated_SOC * 100, 'r');
legend('Coulomb Counting', 'EKF Estimated');
ylabel('SOC [%]');
xlabel('Time [hr]');
title('SOC Comparison');
grid on;

% 5. SOC Error
figure;
plot(Li_ion_Battery.RecordingTime_Hours, Li_ion_Battery.Measured_SOC - Estimated_SOC * 100, 'k');
ylabel('SOC Error [%]');
xlabel('Time [hr]');
title('SOC Estimation Error');
grid on;

function [Estimated_SOC, Estimated_Vt, Vt_Error] = EKF_function(Current, Vt_Actual, Temperature)
    
    % Accessing the Battery Parameters such as R0, R1, C1, Temp etc
    load 'Datasets\BatteryModel.mat';  

    % Accessing the SOC-OCV curve points
    load 'Datasets\SOC-OCV.mat';        
    
    % Initialisation of State of charge, Terminal Voltage and its error
    Estimated_SOC = [];
    Estimated_Vt = [];
    Vt_Error = [];
    
    % Initial SOC = 100% 
    SOC_Init = 1;   

    % Xk = State Vector Matrix [State_of_Charge; V1 (RC branch voltage drop)]
    Xk = [SOC_Init; 0]; 

    % Sampling time (Delta T) in seconds
    Sample_time = 1;        

    % Qn_rated = Battery rated capacity in Ampere-seconds (Ah * 3600)
    Qn_rated = 4.81 * 3600;    

    % f_R0, f_R1, f_C1 = Interpolation functions that return R0, R1, C1 based on current Temperature and SOC. These reflect Thevenin model parameters.
    f_R0 = scatteredInterpolant(Parameters_of_Battery.T, Parameters_of_Battery.SOC, Parameters_of_Battery.R0);
    f_R1 = scatteredInterpolant(Parameters_of_Battery.T, Parameters_of_Battery.SOC, Parameters_of_Battery.R1);
    f_C1 = scatteredInterpolant(Parameters_of_Battery.T, Parameters_of_Battery.SOC, Parameters_of_Battery.C1);
    
    % SOCOCV = Polynomial fit for the OCV vs SOC curve (11th-order fit)
    SOCOCV = polyfit(SOC_OCV.SOC, SOC_OCV.OCV, 11);

    % dSOC_OCV = Derivative of OCV-SOC polynomial: used in Jacobian    
    dSOC_OCV = polyder(SOCOCV);

    % Nx = Length of the State Vector
    n_x = length(Xk);

    % Rx = Measurement noise Covariance (assuming to be 0.2)
    R_x = 2e-1;   

    % Px = State error Covariance Matrix (Assuming to be 0.025 and 0.01 )
    P_x = [0.025, 0;
        0, 0.01];    

    % Qx = Process Noise covariance matric (Assuming to be 0.0004)
    Q_x = [1.0e-4, 0;
        0, 1.0e-4];

    % Finding length of the current values to run the for loop
    ik = length(Current); 

    for k = 1:ik
        T = Temperature(k);
        I_k = Current(k);
        SOC = Xk(1);
        V1 = Xk(2);
  
        % Value of R0, R1, C1 in partiular itertative data substituted in
        % interpolated function
        R0 = f_R0(T, SOC);
        R1 = f_R1(T, SOC);
        C1 = f_C1(T, SOC);

        OCV = polyval(SOCOCV, SOC);
        dOCV = polyval(dSOC_OCV, SOC);

        % Time_constant is Tau = RC
        Time_constant = R1*C1;
        a1 = exp(-Sample_time / Time_constant);
        b1 = R1 * (1 - a1);

        % Terminal Voltage - Vt
        TerminalVoltage = OCV - R0 * I_k - V1;
    
        % Assume that charging and Discharging take place 100% Efficiently
        % (No losses)
        eff = 1; 

        % Jacobian of measurement function OCV
        C_x = [dOCV, -1];
        Error_x = Vt_Actual(k) - TerminalVoltage;

        Estimated_SOC = [Estimated_SOC; Xk(1)];
        Estimated_Vt = [Estimated_Vt; TerminalVoltage];
        Vt_Error = [Vt_Error; Error_x];

        % A = State transition matrix
        A = [1, 0;
             0, a1];

        % B = Control input matrix based on battery model (current input)
        B = [-(eff * Sample_time / Qn_rated); b1];

        % Xk = A * Xk + B * I_k updates the state prediction based on system dynamics
        Xk = A * Xk + B * I_k;

        % P_x = Covariance prediction update step
        P_x = A * P_x * A' + Q_x;

        % Calulaction of Kalman Gain
        K = P_x * C_x' / (C_x * P_x * C_x' + R_x);

        % Based on Kalman Gain the previous state is adjusted
        Xk = Xk + K * Error_x;

        % P_x update = Covariance matrix correction after measurement assimilation
        P_x = (eye(n_x) - K * C_x) * P_x;
    end
end
