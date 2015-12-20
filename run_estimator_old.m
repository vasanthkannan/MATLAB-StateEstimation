format long

% Read state estimator data input file
load ('StateEstimatorData')


tolerance  = PowerFlowInputData.powerflow_tolerance;
numbus     = PowerFlowInputData.numbus;
numline    = PowerFlowInputData.numline;
numgen     = PowerFlowInputData.numgen;
numarea    = PowerFlowInputData.numarea;
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

% % Initial values for Y supporting parallel lines
% Y=zeros(numbus,numbus);
% for j = 1:numline
%   Y(frombus(j),tobus(j))=Y(frombus(j),tobus(j))-1/(R(j)+sqrt(-1)*X(j));
%   Y(tobus(j),frombus(j))=Y(tobus(j),frombus(j))-1/(R(j)+sqrt(-1)*X(j));
% end
% for j = 1:numbus
%   Y(j,j) = -sum(Y(j,:));
% end
% 
% % Getting the G, B before counting in Bcap. Since in our flow measurement,
% % H_flow, the B and G in the formula do not include Bcap.
% G0=-real(Y);
% B0=-imag(Y);
% 
% for j = 1:numline
%   Y(frombus(j),frombus(j))=Y(frombus(j),frombus(j))+sqrt(-1)*Bcap(j)/2;
%   Y(tobus(j),tobus(j))=Y(tobus(j),tobus(j))+sqrt(-1)*Bcap(j)/2;
% end

% G=real(Y);
% %B=imag(Y); we used defintion Y=G-jB in this program.
% B=-imag(Y);

G0 = G;
B0 = B;

numVmeas =   SCADA_voltage_meas_data.number_voltage_meas;
Vmeasbus =   SCADA_voltage_meas_data.voltage_meas_bus;
Vmeassigma = SCADA_voltage_meas_data.voltage_meas_sigma;
Vmeasvalue = SCADA_voltage_meas_data.voltage_meas_value;

numAmeas =   SCADA_angle_meas_data.number_angle_meas;
Ameasbus =   SCADA_angle_meas_data.angle_meas_bus;
Ameassigma = SCADA_angle_meas_data.angle_meas_sigma;
Ameasvalue = SCADA_angle_meas_data.angle_meas_value;

numImeas    = SCADA_injection_meas_data.number_injection_meas;
Imeasbus    = SCADA_injection_meas_data.injection_meas_bus;
ImeasPsigma = SCADA_injection_meas_data.injection_meas_Psigma;
ImeasQsigma = SCADA_injection_meas_data.injection_meas_Qsigma;
ImeasPvalue = SCADA_injection_meas_data.injection_meas_Pvalue;
ImeasQvalue = SCADA_injection_meas_data.injection_meas_Qvalue;


numFmeas     = SCADA_flow_meas_data.number_flow_meas;
Fmeasbranch  = SCADA_flow_meas_data.flow_meas_branch;
Fmeasfrombus = SCADA_flow_meas_data.flow_meas_frombus;
Fmeastobus   = SCADA_flow_meas_data.flow_meas_tobus;
FmeasPsigma  = SCADA_flow_meas_data.flow_meas_Psigma;
FmeasQsigma  = SCADA_flow_meas_data.flow_meas_Qsigma;
FmeasPvalue = SCADA_flow_meas_data.flow_meas_Pvalue;
FmeasQvalue = SCADA_flow_meas_data.flow_meas_Qvalue;


%significant level of the hypothesis test
alpha=0.001;

numFPmeas=numFmeas;
numFQmeas=numFmeas;
numIPmeas=numImeas;
numIQmeas=numImeas;
FPmeasbranch =Fmeasbranch;
FQmeasbranch =Fmeasbranch;
FPmeasfrombus=Fmeasfrombus;
FQmeasfrombus=Fmeasfrombus;
FPmeastobus=Fmeastobus;
FQmeastobus=Fmeastobus;
IPmeasbus=Imeasbus;  
IQmeasbus=Imeasbus;  

%Give ramdom value to TJ and J, and make sure J>TJ, so that we can enter
%the threshold testing iteration loop. Later on TJ and J are determined by
%the state estimator and the bad measurement detector.

buildHmatrix

%stop

TJ=0;
J=1;

while J>TJ
   %Begin the detection and Identification Iteration Loop
   Nm=numFPmeas+numFQmeas+numIPmeas+numIQmeas;    % number of measuremets
   Ns=2*(numbus-1);                               % number of states
   K=Nm-Ns;                                       %degree of freedom
   %Calculate TJ, the threshold for J(x).
   TJ=2*K-chi2inv(alpha,K);
 

   maxiter=20;
   % degrees_per_radian= 57.29577951;
   iteration=1;
   badmeasurementbus=0;
   % Build data structure for voltage mag and angle and install flat Start
   for i=1:numbus
       z.vmag(i)=1.0;
       z.vang(i)=0.0;
   end

   for i = 1:numVmeas 
       z.vmag(Vmeasbus(i)) = Vmeasvalue(i); 
   end 

   % Begin the iteration loop
   while (iteration<maxiter)
       % Build residual J
       J=zeros(4,1); %Since Vmeasbus=1, Ameasbus=1, which means they are 
                     %reference bus, there's no need to have residual for 
                     %voltage and angle measurement.
    
       % Residual for flow measurements
       for i=1:numbus
           vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
       end
       for i = 1:numFPmeas
           a=FPmeasfrombus(i);
           b=FPmeastobus(i);
           Iline = (vr(a) - vr(b))/(R(FPmeasbranch(i))+sqrt(-1)*X(FPmeasbranch(i)));
           Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FPmeasbranch(i))*0.5);
           Sline = vr(a) * conj(Iline);
           FcalcPvalue(i) = real(Sline);
           J1(i)=(FmeasPvalue(i)- FcalcPvalue(i))^2/FmeasPsigma(i)^2;
           J(1)=J1(i)+J(1);
       end
       for i = 1:numFQmeas
           a=FQmeasfrombus(i);
           b=FQmeastobus(i);
           Iline = (vr(a) - vr(b))/(R(FQmeasbranch(i))+sqrt(-1)*X(FQmeasbranch(i)));
           Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FQmeasbranch(i))*0.5);
           Sline = vr(a) * conj(Iline);
           FcalcQvalue(i)= imag(Sline);
           J2(i)=(FmeasQvalue(i)-FcalcQvalue(i))^2/FmeasQsigma(i)^2;
           J(2)=J(2)+J2(i);
       end    
       
       % Residual for ImeasPvalue
       IcalcPvalue=zeros(numIPmeas,1);
       IcalcQvalue=zeros(numIQmeas,1);
       for i=1:numIPmeas
           a = IPmeasbus(i);
           for b=1:numbus
                Pinj =z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
                IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
           end
           J3(i)=(ImeasPvalue(i)-IcalcPvalue(i))^2/ImeasPsigma(i)^2;
           J(3)=J(3)+J3(i);
       end
       for i=1:numIQmeas
           a = IQmeasbus(i);
           for b=1:numbus
               Qinj =z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
               IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
           end
           J4(i)=(ImeasQvalue(i)-IcalcQvalue(i))^2/ImeasQsigma(i)^2;
           J(4)=J(4)+J4(i);
       end
       
       J=sum(J);
 
       numbus2=2*numbus;
       nummeas=numIPmeas+numIQmeas+numFPmeas+numFQmeas;
       % Build Rweighting matrix
       Rweighting=zeros(nummeas);
       for i=1:numFPmeas
            Rweighting(i,i)=1/FmeasPsigma(i);
       end
       for i=(numFPmeas+1):(numFPmeas+numFQmeas)
            Rweighting(i,i)=1/FmeasQsigma(i-numFPmeas);
       end
       for i=(numFPmeas+numFQmeas+1):(numFPmeas+numFQmeas+numIPmeas)
            Rweighting(i,i)=1/ImeasPsigma(i-numFPmeas-numFQmeas);
       end
       for i=(numFPmeas+numFQmeas+numIPmeas+1):(numFPmeas+numFQmeas+numIPmeas+numIQmeas)
            Rweighting(i,i)=1/ImeasQsigma(i-numFPmeas-numFQmeas-numIPmeas);
       end
    
       % Build h
       h=sparse(nummeas,(numbus-1)*2);
       % Build P flow part of h matrix
       H_pflow=sparse(numFPmeas,(numbus-1)*2);
       for i=1:numFPmeas
            a=FPmeasfrombus(i);
            b=FPmeastobus(i);
            if a==1  % reference bus
               H_pflow(i,2*(b-1)-1)=-cos(-z.vang(b))*G0(a,b)-sin(-z.vang(b))*B0(a,b);
               H_pflow(i,2*(b-1))=-z.vmag(b)*(sin(-z.vang(b))*G0(a,b)-cos(-z.vang(b))*B0(a,b));
            elseif b==1
               H_pflow(i,2*(a-1)-1)=2*z.vmag(a)*G0(a,b)-(cos(z.vang(a))*G0(a,b)+sin(z.vang(a))*B0(a,b));
               H_pflow(i,2*(a-1))=z.vmag(a)*(sin(z.vang(a))*G0(a,b)-cos(z.vang(a))*B0(a,b));
            else
               H_pflow(i,2*(a-1)-1)=2*z.vmag(a)*G0(a,b)-z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G0(a,b)+sin(z.vang(a)-z.vang(b))*B0(a,b));
               H_pflow(i,2*(a-1))=-z.vmag(a)*z.vmag(b)*(-sin(z.vang(a)-z.vang(b))*G0(a,b)+cos(z.vang(a)-z.vang(b))*B0(a,b));
               H_pflow(i,2*(b-1)-1)=-z.vmag(a)*(cos(z.vang(a)-z.vang(b))*G0(a,b)+sin(z.vang(a)-z.vang(b))*B0(a,b));
               H_pflow(i,2*(b-1))=-z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G0(a,b)-cos(z.vang(a)-z.vang(b))*B0(a,b));
            end
       end
    
       % Build Q flow part of h matrix
       H_Qflow=sparse(numFQmeas,(numbus-1)*2);
       for i=1:numFQmeas
             a=FQmeasfrombus(i);
             b=FQmeastobus(i);
             if a==1  % reference bus
                H_Qflow(i,2*(b-1)-1)=-sin(-z.vang(b))*G0(a,b)+cos(-z.vang(b))*B0(a,b);
                H_Qflow(i,2*(b-1))=-z.vmag(b)*(-cos(-z.vang(b))*G0(a,b)-sin(-z.vang(b))*B0(a,b));
             elseif b==1
                H_Qflow(i,2*(a-1)-1)=-2*z.vmag(a)*(B0(a,b)+Bcap(FQmeasbranch(i)))-(sin(z.vang(a))*G0(a,b)-cos(z.vang(a))*B0(a,b));
                H_Qflow(i,2*(a-1))=-z.vmag(a)*(cos(z.vang(a))*G0(a,b)+sin(z.vang(a))*B0(a,b));
             else
                H_Qflow(i,2*(a-1)-1)=-2*z.vmag(a)*(B0(a,b)+Bcap(FQmeasbranch(i)))-z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G0(a,b)-cos(z.vang(a)-z.vang(b))*B0(a,b));
                H_Qflow(i,2*(a-1))=-z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G0(a,b)+sin(z.vang(a)-z.vang(b))*B0(a,b));
                H_Qflow(i,2*(b-1)-1)=-z.vmag(a)*sin(z.vang(a)-z.vang(b))*G0(a,b)+cos(z.vang(a)-z.vang(b))*B0(a,b);
                H_Qflow(i,2*(b-1))=-z.vmag(a)*z.vmag(b)*(-cos(z.vang(a)-z.vang(b))*G0(a,b)-sin(z.vang(a)-z.vang(b))*B0(a,b));
             end
       end
    % Build injection part of h matrix
    H_Pinj=zeros(numIPmeas,numbus2-2);
    magnitudeP=zeros(numbus,1);
    angleP=zeros(numbus,1);
    % Build Pinj part of h.
    for i=1:numIPmeas
        a = IPmeasbus(i);
        for b=2:numbus
            if b==a
               for m=1:numbus
                   if m==b
                      magnitudeP(m)=2*z.vmag(m)*G(b,b);
                      angleP(m)=0; 
                   else
                      magnitudeP(m)=z.vmag(m)*(cos(z.vang(a)-z.vang(m))*G(a,m)-sin(z.vang(a)-z.vang(m))*B(a,m));
                      angleP(m)=z.vmag(a)*z.vmag(m)*(-sin(z.vang(b)-z.vang(m))*G(b,m)-cos(z.vang(b)-z.vang(m))*B(b,m));
                   end
               end
               H_Pinj(i,2*(b-1)-1)=sum(magnitudeP);
               H_Pinj(i,2*(b-1))=sum(angleP);
            else
               H_Pinj(i,2*(b-1)-1)=z.vmag(a)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
               H_Pinj(i,2*(b-1))=z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
            end
        end
    end  
        
    H_Qinj=zeros(numIQmeas,numbus2-2);
    magnitudeQ=zeros(numbus,1);
    angleQ=zeros(numbus,1);
    % Build Qinj part of h.
    for i=1:numIQmeas
        a = IQmeasbus(i);
            for b=2:numbus
                if b==a 
                    for m=1:numbus
                        if m==b
                           magnitudeQ(m)=2*z.vmag(m)*B(b,b);
                           angleQ(m)=0; 
                        else
                           magnitudeQ(m)=z.vmag(m)*(sin(z.vang(a)-z.vang(m))*G(a,m)+cos(z.vang(a)-z.vang(m))*B(a,m));
                           angleQ(m)=z.vmag(a)*z.vmag(m)*(cos(z.vang(a)-z.vang(m))*G(b,m)-sin(z.vang(a)-z.vang(m))*B(a,m));
                        end
                    end
                    H_Qinj(i,2*(b-1)-1)=sum(magnitudeQ);
                    H_Qinj(i,2*(b-1))=sum(angleQ);
                else
                    H_Qinj(i,2*(b-1)-1)=z.vmag(a)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
                    H_Qinj(i,2*(b-1))=z.vmag(a)*z.vmag(b)*(-cos(z.vang(a)-z.vang(b))*G(a,b)+sin(z.vang(a)-z.vang(b))*B(a,b));
                end
            end
        end
        
        h=[ H_pflow
            H_Qflow
            H_Pinj
            H_Qinj];
        
        % build_h_transpose_r_inverse_h
        htrih=transpose(h)*inv(Rweighting)*h;
     
        % build delta z
        dZ_Pflow=zeros(numFPmeas,1);
        dZ_Qflow=zeros(numFQmeas,1);
        for i=1:numFPmeas
            a=Fmeasfrombus(i);
            b=Fmeastobus(i);
            dZ_Pflow(i)=FmeasPvalue(i)-FcalcPvalue(i);
        end
        for i=1:numFQmeas
            a=Fmeasfrombus(i);
            b=Fmeastobus(i); 
            dZ_Qflow(i)=FmeasQvalue(i)-FcalcQvalue(i);
        end
        % injection part
        dZ_Pinj=zeros(numIPmeas,1);
        dZ_Qinj=zeros(numIQmeas,1);
        for i=1:numIPmeas
            dZ_Pinj(i)=ImeasPvalue(i)-IcalcPvalue(i);
        end
        for i=1:numIQmeas
            dZ_Qinj(i)=ImeasQvalue(i)-IcalcQvalue(i);
        end
        dz=[dZ_Pflow
            dZ_Qflow
            dZ_Pinj
            dZ_Qinj];
        
        
        %Solve for delta x
        dX=inv(htrih)*transpose(h)*inv(Rweighting)*dz;
        
        %measurement error covariance matrix
        E=diag(Rweighting-h*inv(htrih)*transpose(h));
        for i=1:numFPmeas
            sigmaPflow(i)=E(i);
        end
        for i=1:numFQmeas
            sigmaQflow(i)=E(i+numFPmeas);
        end
        for i=1:numIPmeas
            sigmaPinjec(i)=E(i+numFPmeas+numFQmeas);
        end
        for i=1:numIQmeas
            sigmaQinjec(i)=E(i+numFPmeas+numFQmeas+numIPmeas);
        end
        
        
        %update_voltage_vector
        
        if max(abs(dX))>0.0000000000000001          %set the critiriar which helps us decide when to stop the iterations.
            
            
            for i=2:numbus
                z.vmag(i)=z.vmag(i)+dX((i-1)*2-1);
                z.vang(i)=z.vang(i)+dX((i-1)*2);
            end
            
        else break
        end
        z.vang;
        z.vmag;
        
        
    % Caculate Residual again at the end of the iteration loop
    % Build residual J
    J=zeros(4,1); %Since Vmeasbus=1, Ameasbus=1, which means they are reference bus, there's no need to have residual for voltage and angle measurement.
    
    % Residual for flow measurements
    for i=1:numbus
         vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
    end
    for i = 1:numFPmeas
        a=FPmeasfrombus(i);
        b=FPmeastobus(i);
        Iline = (vr(a) - vr(b))/(R(FPmeasbranch(i))+sqrt(-1)*X(FPmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FPmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcPvalue(i) = real(Sline);
        J1(i)=(FmeasPvalue(i)- FcalcPvalue(i))^2/FmeasPsigma(i)^2;
        J(1)=J1(i)+J(1);
    end
    for i = 1:numFQmeas
        a=FQmeasfrombus(i);
        b=FQmeastobus(i);
        Iline = (vr(a) - vr(b))/(R(FQmeasbranch(i))+sqrt(-1)*X(FQmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FQmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcQvalue(i)= imag(Sline);
        J2(i)=(FmeasQvalue(i)-FcalcQvalue(i))^2/FmeasQsigma(i)^2;
        J(2)=J(2)+J2(i);
    end
    
    % Residual for ImeasPvalue
    IcalcPvalue=zeros(numIPmeas,1);
    IcalcQvalue=zeros(numIQmeas,1);
    for i=1:numIPmeas
        a = IPmeasbus(i);
        for b=1:numbus
            Pinj =z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
            IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
        end
        J3(i)=(ImeasPvalue(i)-IcalcPvalue(i))^2/ImeasPsigma(i)^2;
        J(3)=J(3)+J3(i);
    end
    for i=1:numIQmeas
        a = IQmeasbus(i);
        for b=1:numbus
            Qinj =z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
            IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
        end
        J4(i)=(ImeasQvalue(i)-IcalcQvalue(i))^2/ImeasQsigma(i)^2;
        J(4)=J(4)+J4(i);
    end
    J=sum(J);
   
    
    iteration=iteration+1;
        
    drawnow;
        
    end %end of state estimator iteration loop
    
    if J<=TJ
        break
    else
      % Detect and Identify bad measurements.
      % Calculate all f(i),y_norm(i) values for Nm measurements.
      for i=1:numbus
         vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
      end
      for i = 1:numFPmeas
        a=FPmeasfrombus(i);
        b=FPmeastobus(i);
        Iline = (vr(a) - vr(b))/(R(FPmeasbranch(i))+sqrt(-1)*X(FPmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FPmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcPvalue(i) = real(Sline);
        if sigmaPflow(i)>0.000001
            y1_norm(i)=(FmeasPvalue(i)-FcalcPvalue(i))/sigmaPflow(i);
        else 
            y1_norm(i)=0;
        end
      end
      y1_norm=y1_norm(1:numFPmeas);
    for i = 1:numFQmeas
        a=FQmeasfrombus(i);
        b=FQmeastobus(i);
        Iline = (vr(a) - vr(b))/(R(FQmeasbranch(i))+sqrt(-1)*X(FQmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FQmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcQvalue(i)= imag(Sline);
        if sigmaQflow(i)>0.000001
            y2_norm(i)=(FmeasQvalue(i)-FcalcQvalue(i))/sigmaQflow(i);
        else
            y2_norm(i)=0;
        end
    end
    y2_norm=y2_norm(1:numFQmeas);
    
    if numImeas==0
        y3_norm=[];
        y4_norm=[];
    else
        IcalcPvalue=zeros(numIPmeas,1);
        IcalcQvalue=zeros(numIQmeas,1);
        for i=1:numIPmeas
            a = IPmeasbus(i);
            for b=1:numbus
                Pinj =z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
                IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
            end
            if sigmaPinjec(i)>0.000001
               y3_norm(i)=(ImeasPvalue(i)-IcalcPvalue(i))/sigmaPinjec(i);
            else
               y3_norm(i)=0;
            end
        end
        y3_norm=y3_norm(1:numIPmeas);
        for i=1:numIQmeas
            a = IQmeasbus(i);
            for b=1:numbus
                Qinj = z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
                IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
            end
            if sigmaQinjec(i)>0.000001
               y4_norm(i)=(ImeasQvalue(i)-IcalcQvalue(i))/sigmaQinjec(i);
            else
                y4_norm(i)=0;
            end
        end
         y4_norm=y4_norm(1:numIQmeas);
   end
    % Assume that numFPmeas, numFQmeas never go to zero.
    if numIPmeas==0
        y_norm=[y1_norm,y2_norm,y4_norm];
    elseif numIQmeas==0
        y_norm=[y1_norm,y2_norm,y3_norm];
    elseif numIPmeas==0&&numIQmeas==0
        y_norm=[y1_norm,y2_norm];
    else
        y_norm=[y1_norm,y2_norm,y3_norm,y4_norm];
    end
    
    
     % Display
        fprintf(' %s \n',' ');
        fprintf(' %s \n','Bus Voltage and Phase Angle Measurements: ');
        fprintf(' %s \n','-------------------------------------------------------------------------------------');
        fprintf(' %s \n',' From      Measured        Calculated      Difference     Covariance of      ');
        fprintf(' %s \n',' Bus       Value           Value                          Meansurement       ');
        fprintf(' %s \n','-------------------------------------------------------------------------------------');
        
         fprintf(' %s \n','Bus Voltage: ');
        for i=1:numVmeas
            fprintf(' %5d   %8.3f        %8.3f        %8.3f      %8.3f  \n',Vmeasbus(i), Vmeasvalue(i), z.vmag(Vmeasbus(i)),(Vmeasvalue(i) - z.vmag(Vmeasbus(i))), Vmeassigma(i));
        end
        
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');
        
        fprintf(' %s \n','Bus Phase Angle: ');
        
        for i=1:numAmeas  
            fprintf(' %5d   %8.3f        %8.3f        %8.3f      %8.3f  \n',Ameasbus(i), Ameasvalue(i), z.vang(Ameasbus(i)),(Ameasvalue(i) - z.vang(Ameasbus(i))), Ameassigma(i));
        end
 
 
        fprintf(' %s \n',' ');
        fprintf(' %s \n','-------------------------------------------------------------------------------------');
        fprintf(' %s \n',' From      To      Measured        Calculated      Difference     Covariance of      ');
        fprintf(' %s \n',' Bus       Bus     Value           Value                          Meansurement       ');
        fprintf(' %s \n','-------------------------------------------------------------------------------------');
        
        fprintf(' %s \n','Pflow: ');
        for i=1:numFPmeas
            fprintf(' %5d    %5d    %8.3f        %8.3f       %8.3f       %8.3f \n',FPmeasfrombus(i),FPmeastobus(i),FmeasPvalue(i)*baseMVA, FcalcPvalue(i)*baseMVA,(FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,FmeasPsigma(i));
        end
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');
        
        fprintf(' %s \n','Qflow: ');
        for i=1:numFQmeas
            fprintf(' %5d    %5d    %8.3f        %8.3f       %8.3f       %8.3f \n',FQmeasfrombus(i),FQmeastobus(i),FmeasQvalue(i)*baseMVA, FcalcQvalue(i)*baseMVA,(FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,FmeasQsigma(i));
        end
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');

       
        fprintf(' %s \n','Pinjection: ');
        for i=1:numIPmeas
            fprintf(' %5d           %8.3f        %8.3f        %8.3f      %8.3f  \n',IPmeasbus(i), ImeasPvalue(i)*baseMVA, IcalcPvalue(i)*baseMVA,(ImeasPvalue(i)-IcalcPvalue(i))*baseMVA, ImeasPsigma(i));
        end
        
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');
        fprintf(' %s \n','Qinjection: ');
        
        for i=1:numIQmeas  
            fprintf(' %5d           %8.3f        %8.3f        %8.3f      %8.3f  \n',IQmeasbus(i), ImeasQvalue(i)*baseMVA, IcalcQvalue(i)*baseMVA,(ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,ImeasQsigma(i));
        end
 

       fprintf(' %s \n','');
       fprintf(' %s \n','');
       fprintf(' %s \n','');
       fprintf(' %s \n','-------------------');
       fprintf(' %s \n','Bad Data Analysis');
       fprintf(' %s \n','-----------------------------------------------------------------------');
       fprintf(' %s \n','Measurement  Measurement    Caculated  Difference   Covariance   y_norm ');
       fprintf(' %s \n','Number       Value          Value                   of error            ');
       fprintf(' %s \n','-----------------------------------------------------------------------');
       fprintf(' %s\n','Pflow measurements:');
       for i=1:numFPmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasPvalue(i)*baseMVA,  FcalcPvalue(i)*baseMVA, (FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,sigmaPflow(i),y1_norm(i)*baseMVA);
       end
       fprintf(' %s\n','Qflow measurements:');
       for i=1:numFQmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasQvalue(i)*baseMVA,  FcalcQvalue(i)*baseMVA, (FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,sigmaQflow(i),y2_norm(i)*baseMVA);
       end 
       fprintf(' %s\n','Pinjection measurements:');
       for i=1:numIPmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', IPmeasbus(i), ImeasPvalue(i)*baseMVA,  IcalcPvalue(i)*baseMVA, (ImeasPvalue(i)-IcalcPvalue(i))*baseMVA,sigmaPinjec(i),y3_norm(i)*baseMVA);
       end
       fprintf(' %s\n','Qinjection measurements:');
       for i=1:numIQmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', IQmeasbus(i), ImeasQvalue(i)*baseMVA,  IcalcQvalue(i)*baseMVA, (ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,sigmaQinjec(i),y4_norm(i)*baseMVA);
       end 
    %find the maximum y_norm(i), which is the prime suspect of error,both
    %the value and the index, putting them into Primesuspectvalue and
    %Primesuspectindex seperately.
    
    
    [Primesuspectvalue, Primesuspectindex]=max(abs(y_norm))
    if Primesuspectvalue<0.03
        fprintf('There is no bad measurement anymore!\n');
        break;
    end
        
    % Remove the bad measurement and display the bad measurement index and
    % value.
    if Primesuspectindex<=numFPmeas    % the removed measurements is a P flow measurement
        fprintf('The bad measurement is a P flow measurement. %e\n');
        fprintf('The bad measurement occurs at(the index of P flow measurements): %2d    %s', Primesuspectindex,'');
        fprintf('Its value is:%e',Primesuspectvalue);
        fprintf(' %s \n',' ');
        numFPmeas=numFPmeas-1;    
        if numFPmeas~=0
          for i=Primesuspectindex:(numFPmeas)
            FmeasPvalue(i)=FmeasPvalue(i+1);
            FmeasPsigma(i)=FmeasPsigma(i+1);
            FPmeasbranch(i)=FPmeasbranch(i+1);
            FPmeasfrombus(i)=FPmeasfrombus(i+1);
            FPmeastobus(i)=FPmeastobus(i+1);
          end
        end
    elseif Primesuspectindex>numFPmeas&&Primesuspectindex<=(numFPmeas+numFQmeas)
        fprintf('The bad measurement is a Q flow measurement. %e\n');
        fprintf('The bad measurement occurs at (the index of Q flow measurements):  %2d    %s', (Primesuspectindex-numFPmeas),'');
        fprintf('Its value is:%e',Primesuspectvalue);
        fprintf(' %s \n',' ');
        numFQmeas=numFQmeas-1;
        if numFQmeas~=0
          for i=(Primesuspectindex-numFPmeas):(numFQmeas)
            FmeasQvalue(i)=FmeasQvalue(i+1);
            FmeasQsigma(i)=FmeasQsigma(i+1);
            FQmeasbranch(i)=FQmeasbranch(i+1);
            FQmeasfrombus(i)=FQmeasfrombus(i+1);
            FQmeastobus(i)=FQmeastobus(i+1);
          end
        end
    elseif Primesuspectindex>(numFPmeas+numFQmeas) && Primesuspectindex<=(numFPmeas+numFQmeas+numIPmeas)
        badmeasurementbus=Primesuspectindex-numFPmeas-numFQmeas;
        fprintf('The bad measurement is a P injection measurement. %e\n');
        fprintf('The bad measurement occurs at bus: %2d    %s',  IPmeasbus(badmeasurementbus),'');
        fprintf('Its value is:%2d',Primesuspectvalue);
        fprintf(' %s \n',' ');
        numIPmeas=numIPmeas-1;
        if numIPmeas~=0
          for i=(Primesuspectindex-numFPmeas-numFQmeas):(numIPmeas)
            ImeasPvalue(i)=ImeasPvalue(i+1);
            ImeasPsigma(i)=ImeasPsigma(i+1);
            IPmeasbus(i)=IPmeasbus(i+1);
          end
          IPmeasbus=IPmeasbus(1:numIPmeas);
        end
    else
        badmeasurementbus=Primesuspectindex-numFPmeas-numFQmeas-numIPmeas;
        fprintf('The bad measurement is a Q injection measurement. %e\n');
        fprintf('The bad measurement occurs at bus: %2d    %s',  IQmeasbus(badmeasurementbus), '');
        fprintf('Its value is:%2d',Primesuspectvalue);
        fprintf(' %s \n',' ');
        numIQmeas=numIQmeas-1;
        if numIQmeas~=0
          for i=(Primesuspectindex-numFPmeas-numFQmeas-numIPmeas):(numIQmeas)
            ImeasQvalue(i)=ImeasQvalue(i+1);
            ImeasQsigma(i)=ImeasQsigma(i+1);
            IQmeasbus(i)=IQmeasbus(i+1);
          end
          IQmeasbus=IQmeasbus(1:numIQmeas);
        end
    end
          % Calculate all f(i),y_norm(i) values for Nm measurements.
      for i=1:numbus
         vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
      end
      for i = 1:numFPmeas
        a=FPmeasfrombus(i);
        b=FPmeastobus(i);
%         for j = 1:numline
%           if frombus(j) == a & tobus(j)==b
%               break
%           elseif frombus(j) == b & tobus(j)==a
%               break
%           end
%         end
        Iline = (vr(a) - vr(b))/(R(FPmeasbranch(i))+sqrt(-1)*X(FPmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FPmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcPvalue(i) = real(Sline);
        y1_norm(i)=(FmeasPvalue(i)-FcalcPvalue(i))/sigmaPflow(i);
       
      end
       y1_norm=y1_norm(1:numFPmeas);
       
    for i = 1:numFQmeas
        a=FQmeasfrombus(i);
        b=FQmeastobus(i);
%         for j = 1:numline
%           if frombus(j) == a & tobus(j)==b
%               break
%           elseif frombus(j) == b & tobus(j)==a
%               break
%           end
%         end
        Iline = (vr(a) - vr(b))/(R(FQmeasbranch(i))+sqrt(-1)*X(FQmeasbranch(i)));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(FQmeasbranch(i))*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcQvalue(i)= imag(Sline);
        %if FmeasQvalue(i)-FcalcQvalue(i)<0.00001
         % y2_norm(i)=0;
       %else
          y2_norm(i)=(FmeasQvalue(i)-FcalcQvalue(i))/sigmaQflow(i);
          
       % end
    end
    y2_norm=y2_norm(1:numFQmeas);
    
    IcalcPvalue=zeros(numIPmeas,1);
    IcalcQvalue=zeros(numIQmeas,1);
    for i=1:numIPmeas
        a = IPmeasbus(i);
        for b=1:numbus
            Pinj =z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
            IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
        end
        y3_norm(i)=(ImeasPvalue(i)-IcalcPvalue(i))/sigmaPinjec(i);
       
    end
     y3_norm=y3_norm(1:numIPmeas);
    for i=1:numIQmeas
        a = IQmeasbus(i);
        for b=1:numbus
            Qinj = z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
            IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
        end
        y4_norm(i)=(ImeasQvalue(i)-IcalcQvalue(i))/sigmaQinjec(i);
        
    end
    y4_norm=y4_norm(1:numIQmeas);
    % Assume that numFPmeas, numFQmeas never go to zero.
    if numIPmeas==0
        y_norm=[y1_norm,y2_norm,y4_norm];
    elseif numIQmeas==0
        y_norm=[y1_norm,y2_norm,y3_norm];
    elseif numIPmeas==0&&numIQmeas==0
        y_norm=[y1_norm,y2_norm];
    else
        y_norm=[y1_norm,y2_norm,y3_norm,y4_norm];
    end
  end       % end "Detect and Identify bad measurements" process
 
end % end threshold iteration loop

%Recalculate y_norm
for i = 1:numFPmeas
    y1_norm(i)=(FmeasPvalue(i)-FcalcPvalue(i))/sigmaPflow(i);
end
for i=1:numFQmeas
    y2_norm(i)=(FmeasQvalue(i)-FcalcQvalue(i))/sigmaQflow(i);
end
for i=1:numIPmeas
    y3_norm(i)=(ImeasPvalue(i)-IcalcPvalue(i))/sigmaPinjec(i);
end
for i=1:numIQmeas
    y4_norm(i)=(ImeasQvalue(i)-IcalcQvalue(i))/sigmaQinjec(i);
end

    fprintf(' %s \n',' ');
    fprintf ('The final value of J is %d\n' , J)
    fprintf ('The final value of TJ is %d\n' , TJ) 
     % Display
        fprintf(' %s \n',' ');
        fprintf(' %s \n','-----------------------------------------------------------------------------------');
        fprintf(' %s \n','From    To       Measured       Calculated        Difference   Covariance of       ');
        fprintf(' %s \n','Bus     Bus      Value          Value                          Meansurement        ');
        fprintf(' %s \n','-----------------------------------------------------------------------------------');
        
        fprintf(' %s \n','Pflow: ');
        for i=1:numFPmeas
            fprintf(' %5d %5d    %8.3f        %8.3f        %8.3f       %8.3f  \n',FPmeasfrombus(i),FPmeastobus(i),FmeasPvalue(i)*baseMVA, FcalcPvalue(i)*baseMVA,(FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,FmeasPsigma(i));
        end
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');
        
        fprintf(' %s \n','Qflow: ');
        for i=1:numFQmeas
            fprintf(' %5d %5d    %8.3f        %8.3f        %8.3f       %8.3f  \n',FQmeasfrombus(i),FQmeastobus(i),FmeasQvalue(i)*baseMVA, FcalcQvalue(i)*baseMVA,(FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,FmeasQsigma(i));
        end
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');

       
        fprintf(' %s \n','Pinjection: ');
        for i=1:numIPmeas
            fprintf(' %5d           %8.3f        %8.3f       %8.3f       %8.3f     \n',IPmeasbus(i),ImeasPvalue(i)*baseMVA, IcalcPvalue(i)*baseMVA,(ImeasPvalue(i)-IcalcPvalue(i))*baseMVA, ImeasPsigma(i));
        end
        
        fprintf(' %s \n',' ');
        fprintf(' %s \n',' ');
        fprintf(' %s \n','Qinjection: ');
        
        for i=1:numIQmeas  
            fprintf(' %5d           %8.3f        %8.3f       %8.3f       %8.3f     \n',IQmeasbus(i),ImeasQvalue(i)*baseMVA, IcalcQvalue(i)*baseMVA,(ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,ImeasQsigma(i));
        end
 

       fprintf(' %s \n','');
       fprintf(' %s \n','');
       fprintf(' %s \n','');
       fprintf(' %s \n','-------------------');
       fprintf(' %s \n','Bad Data Analysis');
       fprintf(' %s \n','-----------------------------------------------------------------------');
       fprintf(' %s \n','Measurement  Measurement    Caculated  Difference   Covariance   y_norm ');
       fprintf(' %s \n','Number       Value          Value                   of error            ');
       fprintf(' %s \n','-----------------------------------------------------------------------');
       fprintf(' %s\n','Pflow measurements:');
       for i=1:numFPmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasPvalue(i)*baseMVA,  FcalcPvalue(i)*baseMVA, (FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,sigmaPflow(i),y1_norm(i)*baseMVA);
       end
       fprintf(' %s\n','Qflow measurements:');
       for i=1:numFQmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasQvalue(i)*baseMVA,  FcalcQvalue(i)*baseMVA, (FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,sigmaQflow(i),y2_norm(i)*baseMVA);
       end 
       fprintf(' %s\n','Pinjection measurements:');
       for i=1:numIPmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', IPmeasbus(i), ImeasPvalue(i)*baseMVA,  IcalcPvalue(i)*baseMVA, (ImeasPvalue(i)-IcalcPvalue(i))*baseMVA,sigmaPinjec(i),y3_norm(i)*baseMVA);
       end
       fprintf(' %s\n','Qinjection measurements:');
       for i=1:numIQmeas
           fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', IQmeasbus(i), ImeasQvalue(i)*baseMVA,  IcalcQvalue(i)*baseMVA, (ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,sigmaQinjec(i),y4_norm(i)*baseMVA);
       end 
       
    %Caculate all the flow and injection values from the voltages magnitudes and angles obtained from the state estimator. 
    for i=1:numbus
        vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
    end
    
  
    for i = 1:numline
        a=frombus(i);
        b=tobus(i);
        Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcPvalue(i) = real(Sline);
    end
    for i=1:numline
        a=tobus(i);
        b=frombus(i);
        Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcPvalue(i+numline) = real(Sline);
    end
    for i = 1:numline
        a=frombus(i);
        b=tobus(i);
        Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcQvalue(i)= imag(Sline);
    end
    
     for i = 1:numline
        a=tobus(i);
        b=frombus(i);
        Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
        Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
        Sline = vr(a) * conj(Iline);
        FcalcQvalue(i+numline)= imag(Sline);
    end
    IcalcPvalue=zeros(numbus,1);
    IcalcQvalue=zeros(numbus,1);
    for i=1:numbus
        for b=1:numbus
            Pinj =z.vmag(i)*z.vmag(b)*(cos(z.vang(i)-z.vang(b))*G(i,b)-sin(z.vang(i)-z.vang(b))*B(i,b));
            IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
        end
    end
    
    for i=1:numbus
        for b=1:numbus
            Qinj = z.vmag(i)*z.vmag(b)*(sin(z.vang(i)-z.vang(b))*G(i,b)+cos(z.vang(i)-z.vang(b))*B(i,b));
            IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
        end
    end
       
       
        % Display AC state estimator result
         fprintf(' %s \n',' ');
         fprintf(' %s \n','State Estimator Solution:');
         fprintf(' %s \n',' ');
         fprintf(' %s \n','Bus       Volt.    Volt.    Mw        Mvar      To      Mw        Mvar       ');
         fprintf(' %s \n','Number    Mag.     Angle    Injec.    Injec.    Bus     Flow      Flow       ');
         fprintf(' %s \n','------    ------   ------  --------  --------   ----   -------   -------     ');          

        for i=1:numbus
            fprintf('%5d  %8.3f  %8.3f  %8.3f  %8.3f \n',i, z.vmag(i),z.vang(i),IcalcPvalue(i)*baseMVA,IcalcQvalue(i)*baseMVA); 
            for j=1:numline
                a=frombus(j);
                b=tobus(j);
               if b==i
                  fprintf('                                               %5d    %8.3f  %8.3f   \n', a, FcalcPvalue(j+numline)*baseMVA, FcalcQvalue(j+numline)*baseMVA);
               elseif a==i
                  fprintf('                                               %5d    %8.3f  %8.3f   \n', b, FcalcPvalue(j)*baseMVA, FcalcQvalue(j)*baseMVA);
               end
            end
        end
   