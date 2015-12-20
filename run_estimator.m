

% format long
orthogonal_check = menu('Estimator Menu','1) Standard State Estimator','2) Orthogonal State Estimator'); 
detect_bad_data = menu('Estimator Menu','1) Detect Bad Data','2) Ignore Bad Data'); 


% Read power flow input data and power flow solution data
load('PowerFlowOutput')
% Read state estimator data input file
load ('StateEstimatorData')

baseMVA    = PowerFlowInputData.baseMVA;
tolerance  = PowerFlowInputData.powerflow_tolerance;
numbus     = PowerFlowInputData.numbus;
numline    = PowerFlowInputData.numline;
numgen     = PowerFlowInputData.numgen;
numarea    = PowerFlowInputData.numarea;
refbus      = PowerFlowInputData.refbus;
frombus    = PowerFlowInputData.frombus;
tobus      = PowerFlowInputData.tobus;
R          = PowerFlowInputData.R;
X          = PowerFlowInputData.X;
Bcap       = PowerFlowInputData.Bcap;
Bustype    = PowerFlowInputData.Bustype;
Psched     = PowerFlowInputData.Psched;
Qsched     = PowerFlowInputData.Qsched;
Vsched     = PowerFlowInputData.Vsched;
Y          = PowerFlowInputData.Y;
G          = PowerFlowInputData.G;
B          = PowerFlowInputData.B;

%Power Flow Solution:
Pinj =    PowerFlowSolution.Pinj;
Qinj =    PowerFlowSolution.Qinj;
vr   =    PowerFlowSolution.vr;
Vmag =    PowerFlowSolution.Vmag;
Theta =   PowerFlowSolution.Theta;

numVmeas =   SCADA_voltage_meas_data.number_voltage_meas;
Vmeasbus =   SCADA_voltage_meas_data.voltage_meas_bus;
VmeasNAME =  SCADA_voltage_meas_data.voltage_meas_NAMES;
Vmeassigma = SCADA_voltage_meas_data.voltage_meas_sigma;
Vmeasstatus = SCADA_voltage_meas_data.voltage_meas_status;
Vmeas_base_value = SCADA_voltage_meas_data.voltage_meas_base_value;
Vmeas_measurement_value = SCADA_voltage_meas_data.voltage_meas_measurement_value;

% numVmeas
% Vmeasbus
% VmeasNAMErun_estimator
% Vmeassigma
% Vmeasstatus
% Vmeas_base_value
% Vmeas_measurement_value


numAmeas =   SCADA_angle_meas_data.number_angle_meas;
Ameasbus =   SCADA_angle_meas_data.angle_meas_bus;
AmeasNAME =  SCADA_angle_meas_data.angle_meas_NAMES;
Ameassigma = SCADA_angle_meas_data.angle_meas_sigma;
Ameasstatus = SCADA_angle_meas_data.angle_meas_status;
Ameas_base_value = SCADA_angle_meas_data.angle_meas_base_value;
Ameas_measurement_value = SCADA_angle_meas_data.angle_meas_measurement_value;

numImeas    = SCADA_injection_meas_data.number_injection_meas;
Imeasbus    = SCADA_injection_meas_data.injection_meas_bus;
ImeasNAME =  SCADA_injection_meas_data.injection_meas_NAMES;
ImeasPsigma = SCADA_injection_meas_data.injection_meas_Psigma;
ImeasQsigma = SCADA_injection_meas_data.injection_meas_Qsigma;
Imeasstatus = SCADA_injection_meas_data.injection_meas_status;
Imeas_base_Pvalue = SCADA_injection_meas_data.injection_meas_base_Pvalue;
Imeas_base_Qvalue = SCADA_injection_meas_data.injection_meas_base_Qvalue;
Imeas_measurement_Pvalue = SCADA_injection_meas_data.injection_meas_measurement_Pvalue;
Imeas_measurement_Qvalue = SCADA_injection_meas_data.injection_meas_measurement_Qvalue;

numFmeas     = SCADA_flow_meas_data.number_flow_meas;
Fmeasbranch  = SCADA_flow_meas_data.flow_meas_branch;
Fmeasfrombus = SCADA_flow_meas_data.flow_meas_frombus;
Fmeastobus   = SCADA_flow_meas_data.flow_meas_tobus;
FmeasNAME    =  SCADA_flow_meas_data.flow_meas_NAMES;
FmeasPsigma  = SCADA_flow_meas_data.flow_meas_Psigma;
FmeasQsigma  = SCADA_flow_meas_data.flow_meas_Qsigma;
Fmeasstatus  = SCADA_flow_meas_data.flow_meas_status;
Fmeas_base_Pvalue = SCADA_flow_meas_data.flow_meas_base_Pvalue;
Fmeas_base_Qvalue = SCADA_flow_meas_data.flow_meas_base_Qvalue;
Fmeas_measurement_Pvalue = SCADA_flow_meas_data.flow_meas_measurement_Pvalue;
Fmeas_measurement_Qvalue = SCADA_flow_meas_data.flow_meas_measurement_Qvalue;

% num_measurements is the total number of measurements, but they might
% not all be active
num_measurements = numVmeas + numAmeas + 2*numImeas + 2*numFmeas;

% State estimator convergence tolerance
Est_tolerance = 0.0001;
%significant level of the hypothesis test
alpha=0.1;

Vmag = ones(numbus,1);
Theta = zeros(numbus,1);

Norm_resid = zeros(num_measurements,1);

fprintf(' Estimator iteration summary \n');
fprintf(' Iteration      Residual        Number Active    Degrees of    Bad Data Threshold    Largest       Bad\n');
fprintf('                   J            Measurements     Freedom               Tj            Normalized    Measurement\n');
fprintf('                                                                                     Residual      at\n');



% run first estimation then determine if bad data is detected
if orthogonal_check==1
    Estimation_Loop
    DetectandIdentifyBadData
else
    Estimation_Orthogonal_Loop
    %DetectandIdentifyBadData_Orthogonal
    DetectandIdentifyBadData
end
    


%Start Bad Data Detection Loop, if bad data detect falg is set
 % Detect and Identify bad measurements.

if detect_bad_data == 1
    while J > TJ
      if orthogonal_check==1  
        Estimation_Loop
        DetectandIdentifyBadData
      else
        Estimation_Orthogonal_Loop
        %DetectandIdentifyBadData_Orthogonal
        DetectandIdentifyBadData
      end
    end
 end

display('Final State Estimator Result')

%------------------------------------------------------------------------------------------------

fprintf(' %s \n',' ');
fprintf(' %s \n','Measurement     Base Case Value    Measured Value        Estimated Value ');
fprintf(' %s \n','Name Status   kV     MW    MVAR   kV     MW    MVAR   kV      MW    MVAR ');
for i = 1:numbus
   fprintf(' %s %3d     \n','              Bus',i);   
   Vbase = 230.;
   MVABase = 100.; 
   rad2deg = 180/pi;
   
   if numVmeas > 0
       for j = 1:numVmeas
           if Vmeasbus(j) == i 
            fprintf(' %s    %d     %6.1f              %6.1f              %6.1f  \n',VmeasNAME(j,:),Vmeasstatus(j),...
                Vmeas_base_value(j)*Vbase,...
                Vmeas_measurement_value(j)*Vbase, ...
                Vmeas_estimated_value(j)*Vbase);
              end
       end
   end
   
   if numAmeas > 0
       for j = 1:numAmeas
           if Ameasbus(j) == i 
            fprintf(' %s    %d     %6.1f              %6.1f              %6.1f \n',AmeasNAME(j,:),Ameasstatus(j),...
                Ameas_base_value(j)*rad2deg,...
                Ameas_measurement_value(j)*rad2deg, ...
                Ameas_estimated_value(j)*rad2deg);
           end
       end
   end
   
   if numImeas > 0
     for j = 1:numImeas
         if Imeasbus(j) == i 
          fprintf(' %s    %d          %6.1f %6.1f       %6.1f %6.1f        %6.1f %6.1f \n',ImeasNAME(j,:),Imeasstatus(j),...
                Imeas_base_Pvalue(j)*MVABase, Imeas_base_Qvalue(j)*MVABase, ...
                Imeas_measurement_Pvalue(j)*MVABase, Imeas_measurement_Qvalue(j)*MVABase, ...
                Imeas_estimated_Pvalue(j)*MVABase, Imeas_estimated_Qvalue(j)*MVABase);
           end
       end
   end
   
 if numFmeas >0  
   for j = 1:numFmeas
       if Fmeasfrombus(j) == i ;
        fprintf(' %s    %d          %6.1f %6.1f       %6.1f %6.1f        %6.1f %6.1f \n',FmeasNAME(j,:),Fmeasstatus(j),...
                Fmeas_base_Pvalue(j)*MVABase, Fmeas_base_Qvalue(j)*MVABase, ...
                Fmeas_measurement_Pvalue(j)*MVABase, Fmeas_measurement_Qvalue(j)*MVABase, ...
                Fmeas_estimated_Pvalue(j)*MVABase, Fmeas_estimated_Qvalue(j)*MVABase);
      end
   end
 end                                   
  
 end

% end


