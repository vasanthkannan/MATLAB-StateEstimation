%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Reading Measurement Data File:
% Input Data File 
[file,pathname] = uigetfile('*measurementdata.xls*','Select Spreadsheet File');
if (pathname == 0),
    error('You Must Select A Valid Data File')
end
S=file;          % Name of the File that we need to read

fprintf(' Case ID: %s \n',file);
%--------------------------------------------------------------------------
% Read measurement data from the spreadsheet
[voltage_measurements, voltage_meas_NAMES_CHAR ]   = xlsread(S, 'voltage_measurements');

[numrows,N] = size(voltage_meas_NAMES_CHAR);
for i=2:numrows
name                     = cell2mat(voltage_meas_NAMES_CHAR(i,2));
[M, numcols] = size(name);
    for j = 1:numcols
        voltage_meas_NAMES(i-1, j) = name(1,j);
    end
end
numVmeas = numrows - 1;

%--------------------------------------------------------------------------
[angle_measurements, angle_meas_NAMES_CHAR]     = xlsread(S, 'angle_measurements');

[numrows,N] = size(angle_meas_NAMES_CHAR);
for i=2:numrows
name                     = cell2mat(angle_meas_NAMES_CHAR(i,2));
[M, numcols] = size(name);
    for j = 1:numcols
        angle_meas_NAMES(i-1, j) = name(1,j);
    end
end
numAmeas = numrows - 1;

%--------------------------------------------------------------------------
[injection_measurements, injection_meas_NAMES_CHAR] = xlsread(S, 'injection_measurements');

[numrows,N] = size(injection_meas_NAMES_CHAR);
for i=2:numrows
name                     = cell2mat(injection_meas_NAMES_CHAR(i,2));
[M, numcols] = size(name);
    for j = 1:numcols
        injection_meas_NAMES(i-1, j) = name(1,j);
    end
end
numImeas = numrows - 1;

%--------------------------------------------------------------------------
[flow_measurements, flow_meas_NAMES_CHAR]      = xlsread(S, 'flow_measurements');

[numrows,N] = size(flow_meas_NAMES_CHAR);
for i=2:numrows
name                     = cell2mat(flow_meas_NAMES_CHAR(i,2));
[M, numcols] = size(name);
    for j = 1:numcols
        flow_meas_NAMES(i-1, j) = name(1,j);
    end
end
numFmeas = numrows - 1;

%--------------------------------------------------------------------------



% Get the table dimensions from the tables read off spreadsheet
[numVmeas,numVvariables]=size(voltage_measurements);
[numAmeas,numAvariables]=size(angle_measurements);
[numImeas,numIvariables]=size(injection_measurements);
[numFmeas,numFvariables]=size(flow_measurements);

%NOTE: This state estimator depends upon the user placing a Vmeas and
% an Ameas at the slack bus

% Create the data structure for all measurements and fill with spreadsheet data
% NOTE THAT ALL DATA GOES IN AT THIS POINT BUT NOT THE MEAS VALUES
SCADA_voltage_meas_data = struct('number_voltage_meas', numVmeas,...
                                 'voltage_meas_bus',              [voltage_measurements(1:numVmeas,1)],...
                                 'voltage_meas_NAMES',            [voltage_meas_NAMES(1:numVmeas,:)],...
                                 'voltage_meas_sigma',            [voltage_measurements(1:numVmeas,3)],...
                                 'voltage_meas_status',           [voltage_measurements(1:numVmeas,4)],...
                                 'voltage_meas_base_value',       [voltage_measurements(1:numVmeas,5)],...
                                 'voltage_meas_measurement_value',[voltage_measurements(1:numVmeas,6)]    );
 


SCADA_angle_meas_data =   struct('number_angle_meas', numAmeas,...
                                 'angle_meas_bus',              [angle_measurements(1:numAmeas,1)],...
                                 'angle_meas_NAMES',            [angle_meas_NAMES(1:numAmeas,:)],...
                                 'angle_meas_sigma',            [angle_measurements(1:numAmeas,3)],...
                                 'angle_meas_status',           [angle_measurements(1:numVmeas,4)],...
                                 'angle_meas_base_value',       [angle_measurements(1:numVmeas,5)],...
                                 'angle_meas_measurement_value',[angle_measurements(1:numVmeas,6)]    );
if numImeas==0
   SCADA_injection_meas_data = struct('number_injection_meas', numImeas,...
                                       'injection_meas_bus',   [],...
                                       'injection_meas_NAMES', [],...                                       
                                       'injection_meas_Psigma',[],...
                                       'injection_meas_Qsigma',[],...
                                       'injection_meas_status',[],...
                                       'injection_meas_base_Pvalue',[],...
                                       'injection_meas_base_Qvalue',[],...
                                       'injection_meas_measurement_Pvalue',[],...
                                       'injection_meas_measurement_Qvalue',[]    );
                                       
else
    SCADA_injection_meas_data = struct('number_injection_meas', numImeas,...
                                       'injection_meas_bus',               [injection_measurements(1:numImeas,1)],...
                                       'injection_meas_NAMES',             [injection_meas_NAMES(1:numImeas,:)],...
                                       'injection_meas_Psigma',            [injection_measurements(1:numImeas,3)],...
                                       'injection_meas_Qsigma',            [injection_measurements(1:numImeas,4)],...
                                       'injection_meas_status',            [injection_measurements(1:numImeas,5)],...
                                       'injection_meas_base_Pvalue',       [injection_measurements(1:numImeas,6)],...
                                       'injection_meas_base_Qvalue',       [injection_measurements(1:numImeas,7)],...
                                       'injection_meas_measurement_Pvalue',[injection_measurements(1:numImeas,8)],...
                                       'injection_meas_measurement_Qvalue',[injection_measurements(1:numImeas,9)]    );
end

if numFmeas==0
   SCADA_flow_meas_data = struct('number_flow_meas', numFmeas,...
                              'flow_meas_branch', [],...
                              'flow_meas_NAMES', [],... 
                              'flow_meas_frombus',[],...
                              'flow_meas_tobus',  [],...
                              'flow_meas_Psigma', [],...
                              'flow_meas_Qsigma', [],...
                              'flow_meas_status', [],...
                              'flow_meas_base_Pvalue',[],...
                              'flow_meas_base_Qvalue',[],...
                              'flow_meas_measurement_Pvalue',[],...
                              'flow_meas_measurement_Qvalue',[]    );
else
    
    SCADA_flow_meas_data = struct('number_flow_meas', numFmeas,...
                              'flow_meas_branch',            [flow_measurements(1:numFmeas,1)],...
                              'flow_meas_NAMES',             [flow_meas_NAMES(1:numFmeas,:)],...
                              'flow_meas_frombus',           [flow_measurements(1:numFmeas,3)],...
                              'flow_meas_tobus',             [flow_measurements(1:numFmeas,4)],...
                              'flow_meas_Psigma',            [flow_measurements(1:numFmeas,5)],...
                              'flow_meas_Qsigma',            [flow_measurements(1:numFmeas,6)],...
                              'flow_meas_status',            [flow_measurements(1:numFmeas,7)],...
                              'flow_meas_base_Pvalue',       [flow_measurements(1:numFmeas,8)],...
                              'flow_meas_base_Qvalue',       [flow_measurements(1:numFmeas,9)],...
                              'flow_meas_measurement_Pvalue',[flow_measurements(1:numFmeas,10)],...
                              'flow_meas_measurement_Qvalue',[flow_measurements(1:numFmeas,11)]    );
end

% Place data in local arrays for use in calculating values
Vmeasbus =    SCADA_voltage_meas_data.voltage_meas_bus;
VmeasNAME =   SCADA_voltage_meas_data.voltage_meas_NAMES;
Vmeassigma =  SCADA_voltage_meas_data.voltage_meas_sigma;
Vmeasstatus = SCADA_voltage_meas_data.voltage_meas_status;

Ameasbus =    SCADA_angle_meas_data.angle_meas_bus;
AmeasNAME =   SCADA_angle_meas_data.angle_meas_NAMES;
Ameassigma =  SCADA_angle_meas_data.angle_meas_sigma;
Ameasstatus = SCADA_angle_meas_data.angle_meas_status;


Imeasbus    = SCADA_injection_meas_data.injection_meas_bus;
ImeasNAME  =  SCADA_injection_meas_data.injection_meas_NAMES;
ImeasPsigma = SCADA_injection_meas_data.injection_meas_Psigma;
ImeasQsigma = SCADA_injection_meas_data.injection_meas_Qsigma;
Imeasstatus = SCADA_injection_meas_data.injection_meas_status;


Fmeasbranch  = SCADA_flow_meas_data.flow_meas_branch;
FmeasNAME   =  SCADA_flow_meas_data.flow_meas_NAMES;
Fmeasfrombus = SCADA_flow_meas_data.flow_meas_frombus;
Fmeastobus   = SCADA_flow_meas_data.flow_meas_tobus;
FmeasPsigma  = SCADA_flow_meas_data.flow_meas_Psigma;
FmeasQsigma  = SCADA_flow_meas_data.flow_meas_Qsigma;
Fmeasstatus  = SCADA_flow_meas_data.flow_meas_status; 

