% SCADA Program

% load in the power flow data and solution
load('PowerFlowOutput.mat')

%Power Flow Input data
%-------------------------------------------------------------------
Maxiter    = PowerFlowInputData.Maxiter;
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

%Power Flow Solution Data:
%-------------------------------------------------------------------
Pinj =    PowerFlowSolution.Pinj;
Qinj =    PowerFlowSolution.Qinj;
vr   =    PowerFlowSolution.vr;
Vmag =    PowerFlowSolution.Vmag;
Theta =   PowerFlowSolution.Theta;
Bustype = PowerFlowSolution.Bustype;

% Read State estimator measurement data
%-------------------------------------------------------------------
measurement_datainput_Excel

% Create measurement value tables
%-------------------------------------------------------------------
Vmeasvalue =  zeros(numVmeas,1);
Ameasvalue =  zeros(numAmeas,1);
ImeasPvalue = zeros(numImeas,1);
ImeasQvalue = zeros(numImeas,1);
FmeasPvalue = zeros(numFmeas,1);
FmeasQvalue = zeros(numFmeas,1);

%-------------------------------------------------------------------
% get the measurement values from the power flow

RefBus_checkoff = 0;
% store voltage measurements
for i_vmeas = 1:numVmeas
    Vmeas_base_value(i_vmeas)=Vmag(Vmeasbus(i_vmeas));
    if Vmeasstatus(i_vmeas) == 1
        Vmeasvalue(i_vmeas)=Vmag(Vmeasbus(i_vmeas));
        if Vmeasbus(i_vmeas) == refbus
            RefBus_checkoff = 1;
        end
    end
end

if RefBus_checkoff ~= 1
    fprintf(' %s \n','*************************************************');
    fprintf(' %s \n','Ref Bus must be one of the Voltage measurements');
    fprintf(' %s \n','It must have a status of 1');
    fprintf(' %s \n','*************************************************');
end

%get the angle measurements from the power flow
for i_ameas=1:numAmeas
    Ameas_base_value(i_ameas)=Theta(Ameasbus(i_ameas));
    if Ameasstatus(i_ameas) == 1
         Ameasvalue(i_ameas)=Theta(Ameasbus(i_ameas));
    end
end

%get the injection measurements from the power flow
for i_injmeas=1:numImeas
    Imeas_base_Pvalue(i_injmeas)=Pinj(Imeasbus(i_injmeas));
    Imeas_base_Qvalue(i_injmeas)=Qinj(Imeasbus(i_injmeas));
    if Imeasstatus(i_injmeas) == 1
        ImeasPvalue(i_injmeas)=Pinj(Imeasbus(i_injmeas));
        ImeasQvalue(i_injmeas)=Qinj(Imeasbus(i_injmeas));
    end
end

%get the line P and Q flow values from the power flow
for i_flowmeas = 1:numFmeas
        i=Fmeasfrombus(i_flowmeas);
        j=Fmeastobus(i_flowmeas);
        ibranch = Fmeasbranch(i_flowmeas);
        Gline =  R(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bline = -X(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bcapline = Bcap(ibranch)/2;
        Fmeas_base_Pvalue(i_flowmeas) = Gline*Vmag(i)^2 - Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j))-Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));
        Fmeas_base_Qvalue(i_flowmeas)= -1*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j))*Gline - Vmag(i)^2*Bline + Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j))*Bline - (Vmag(i)^2)*Bcapline;
    if Fmeasstatus(i_flowmeas) == 1
        FmeasPvalue(i_flowmeas) = Fmeas_base_Pvalue(i_flowmeas);
        FmeasQvalue(i_flowmeas) = Fmeas_base_Qvalue(i_flowmeas);
    end
end

noise=0;

Scadamenu; % alter measurements for demo

Vmeas_measurement_value  = Vmeasvalue;
Ameas_measurement_value  = Ameasvalue;
Imeas_measurement_Pvalue = ImeasPvalue;
Imeas_measurement_Qvalue = ImeasQvalue;
Fmeas_measurement_Pvalue = FmeasPvalue;
Fmeas_measurement_Qvalue = FmeasQvalue;

% Printout of measurement input and measurements with noise

fprintf(' %s \n','                                                      ');
fprintf(' %s \n','Measurement    Base Case Value      Measured Value          Measurement      Measurement ');
fprintf(' %s \n','Name Status   kV     MW    MVAR   kV      MW      MVAR        Sigma (P)        Sigma (Q) ');
for i = 1:numbus
   fprintf(' %s %3d     \n','              Bus',i);    
   Vbase = 230.;
   MVABase = 100.; 
   rad2deg = 180/pi;
   
   if numVmeas > 0
       for j = 1:numVmeas
           if Vmeasbus(j) == i 
            fprintf(' %s    %d     %6.1f              %6.1f                     %10.1E \n',voltage_meas_NAMES(j,:),Vmeasstatus(j), Vmeas_base_value(j)*Vbase,...
                                                                                          Vmeas_measurement_value(j)*Vbase, Vmeassigma(j) );
           end
       end
   end
   
   if numAmeas > 0
       for j = 1:numAmeas
           if Ameasbus(j) == i 
            fprintf(' %s    %d     %6.1f              %6.1f                     %10.1E \n',angle_meas_NAMES(j,:),Ameasstatus(j), Ameas_base_value(j)*rad2deg,...
                                                                                          Ameas_measurement_value(j)*rad2deg, Ameassigma(j));
           end
       end
   end
   
   if numImeas > 0
       for j = 1:numImeas
           if Imeasbus(j) == i 
            fprintf(' %s    %d          %6.1f %6.1f        %6.1f   %6.1f      %10.1E      %10.1E \n',...
                ImeasNAME(j,:),Imeasstatus(j), Imeas_base_Pvalue(j)*MVABase, Imeas_base_Qvalue(j)*MVABase, ...
                Imeas_measurement_Pvalue(j)*MVABase, Imeas_measurement_Qvalue(j)*MVABase, ImeasPsigma(j), ImeasQsigma(j) );
           end
       end
   end
   
 if numFmeas >0  
   for j = 1:numFmeas
      if Fmeasfrombus(j) == i ;
          fprintf(' %s    %d          %6.1f %6.1f        %6.1f   %6.1f      %10.1E      %10.1E  \n',...
                FmeasNAME(j,:),Fmeasstatus(j), Fmeas_base_Pvalue(j)*MVABase, Fmeas_base_Qvalue(j)*MVABase, ...
                Fmeas_measurement_Pvalue(j)*MVABase, Fmeas_measurement_Qvalue(j)*MVABase, FmeasPsigma(j), FmeasQsigma(j) );
      end
   end
 end        
                                          
end
 

%-----------------------------------------------------------------------------
% Add measurement value tables to measurement data structures to 
% contain the measurement values and store in mat file

SCADA_voltage_meas_data.voltage_meas_status =             Vmeasstatus;
SCADA_voltage_meas_data.voltage_meas_base_value =         Vmeas_base_value;
SCADA_voltage_meas_data.voltage_meas_measurement_value =  Vmeas_measurement_value;

SCADA_angle_meas_data.angle_meas_status =             Ameasstatus;
SCADA_angle_meas_data.angle_meas_base_value =         Ameas_base_value;
SCADA_angle_meas_data.angle_meas_measurement_value =  Ameas_measurement_value;

SCADA_injection_meas_data.injection_meas_status =             Imeasstatus;
SCADA_injection_meas_data.injection_meas_base_Pvalue =        Imeas_base_Pvalue;
SCADA_injection_meas_data.injection_meas_base_Qvalue =        Imeas_base_Qvalue;
SCADA_injection_meas_data.injection_meas_measurement_Pvalue = Imeas_measurement_Pvalue;
SCADA_injection_meas_data.injection_meas_measurement_Qvalue = Imeas_measurement_Qvalue;

SCADA_flow_meas_data.flow_meas_status =             Fmeasstatus;
SCADA_flow_meas_data.flow_meas_base_Pvalue =        Fmeas_base_Pvalue;
SCADA_flow_meas_data.flow_meas_base_Qvalue =        Fmeas_base_Qvalue;
SCADA_flow_meas_data.flow_meas_measurement_Pvalue = Fmeas_measurement_Pvalue;
SCADA_flow_meas_data.flow_meas_measurement_Qvalue = Fmeas_measurement_Qvalue;


%  Save StateEstimatorData.mat file to send to State Estimator

save('StateEstimatorData','SCADA_voltage_meas_data',...
                          'SCADA_angle_meas_data',...
                          'SCADA_injection_meas_data',...
                          'SCADA_flow_meas_data');%who -file StateEstimatorData

