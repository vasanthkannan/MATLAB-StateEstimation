%Begin the detection and Identification Iteration Loop
%num_active_measurements    % number of active measuremets
num_states = (2*numbus-1);   % number of states
K = num_active_measurements - num_states;                 % degrees of freedom
%Calculate TJ, the threshold for J(x).
TJ=2*K-chi2inv(alpha,K);


fprintf('                                %6d         %6d            %10.3f       \n', num_active_measurements, K, TJ );

if detect_bad_data == 1

        % Calculate the Residual Standardizd Deviation
        Resid_std_dev = sqrt( diag( Rcov - H * inv(H'*inv(Rcov)*H)*H' ) );

        %Set slack bus diagonals to 1
        Resid_std_dev(refbus_Vrow_H) = 1;
        Resid_std_dev(refbus_Arow_H) = 1;

        % Calculate the Normalized Residual
        Norm_resid = inv( diag(Resid_std_dev) )*meas_error;

        measurement_index = 1;
        Largest_Norm_Resid = 0;
%         measurement_index
        % Voltage Measurements
        for i_vmeas = 1:numVmeas
            if Vmeasstatus(i_vmeas) == 1
                if i ~= refbus
                    if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                        Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                        Largest_Norm_Resid_Index = i_vmeas;
                        Largest_Norm_Resid_Type = 'Vmeas';
                    end
                elseif i == refbus
                     if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                        Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                        Largest_Norm_Resid_Index = i_vmeas;
                        Largest_Norm_Resid_Type = 'Vmeas';
                     end
                 end
            measurement_index = measurement_index + 1;
            end
        end
% measurement_index

        % Angle Measurements
        for i_ameas = 1:numAmeas
            if Ameasstatus(i_ameas) == 1
                if i ~= refbus
                    if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                        Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                        Largest_Norm_Resid_Index = i_ameas;
                        Largest_Norm_Resid_Type = 'Ameas';
                    end
                elseif i == refbus
                    if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                        Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                        Largest_Norm_Resid_Index = i_ameas;
                        Largest_Norm_Resid_Type = 'Ameas';
                    end                    
                end
             measurement_index = measurement_index + 1;
             end
        end
% measurement_index

        % Injection Measurements

        for i_injmeas=1:numImeas
          if Imeasbus(i_injmeas) ~= refbus
            if Imeasstatus(i_injmeas) == 1
                if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                    Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                    Largest_Norm_Resid_Index = i_injmeas;
                    Largest_Norm_Resid_Type = 'Imeas';
                end
                measurement_index = measurement_index + 1;
                if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                    Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                    Largest_Norm_Resid_Index = i_injmeas;
                    Largest_Norm_Resid_Type = 'Imeas';
                end
                measurement_index = measurement_index + 1;
            end;
          end;
        end
% measurement_index

        % Flow Measurements

        for i_flowmeas = 1:numFmeas
            if Fmeasstatus(i_flowmeas) == 1
                if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                    Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                    Largest_Norm_Resid_Index = i_flowmeas;
                    Largest_Norm_Resid_Type = 'Fmeas';
                end
                measurement_index = measurement_index + 1;
                if abs(Norm_resid(measurement_index))> Largest_Norm_Resid
                    Largest_Norm_Resid = abs(Norm_resid(measurement_index));
                    Largest_Norm_Resid_Index = i_flowmeas;
                    Largest_Norm_Resid_Type = 'Fmeas';
                end
                measurement_index = measurement_index + 1;
            end;
        end
% measurement_index


        %     Largest_Norm_Resid
        %     Largest_Norm_Resid_Index
        %     Largest_Norm_Resid_Type

        % test to see if bad data has been detected
        if J > TJ

            switch Largest_Norm_Resid_Type
                case 'Vmeas'
                    j = Largest_Norm_Resid_Index;
                    fprintf('Bad Data Detected                                                                     %6.3f      %s\n',...
                        Largest_Norm_Resid, VmeasNAME(j,:) );

                    % disp('Largest error was V Measurement')
                    Vmeasstatus(Largest_Norm_Resid_Index) = 0;
                case 'Ameas'
                    j = Largest_Norm_Resid_Index;
                    fprintf('Bad Data Detected                                                                    %6.3f      %s\n',...
                        Largest_Norm_Resid, AmeasNAME(J,:) );
                    %disp('Largest error was A Measurement')
                    Ameasstatus(Largest_Norm_Resid_Index) = 0;
                case 'Imeas'
                    j = Largest_Norm_Resid_Index;
                    fprintf('Bad Data Detected                                                                    %6.3f      %s\n',...
                        Largest_Norm_Resid, ImeasNAME(j,:) );
                    %disp('Largest error was Inj Measurement')
                    Imeasstatus(Largest_Norm_Resid_Index) = 0;
                case 'Fmeas'
                    j = Largest_Norm_Resid_Index;
                    fprintf('Bad Data Detected                                                                    %6.3f      %s\n',...
                        Largest_Norm_Resid, FmeasNAME(j,:) );
                    %disp('Largest error was Flow Measurement')
                    Fmeasstatus(Largest_Norm_Resid_Index) = 0;
            end

        end

end