% Store ones on the diagonal of the H Matrix

% The H Matrix is built with one Vrow for each Voltage meas and 
% with one Arow for each Angle meas 

% The H Matrix is built with Prows(prow) and Qrows (qrow) for the
% injection measurements and flow measurements 

% All rows have Theta columns (thcol) and V columns (vcol). 
num_rows_H = num_measurements;
num_cols_H = 2*numbus;
H = zeros( num_rows_H, num_cols_H );
Rcov = zeros( num_rows_H, num_rows_H );

H_row_num = 0;

% Build the H Matrix rows for Voltage Measurements
for i_vmeas = 1:numVmeas
        i = Vmeasbus(i_vmeas);
        H_row_num = H_row_num + 1;
        vrow = H_row_num;
        vcol = 2*i;
        if Vmeasstatus(i_vmeas) ~= 0
%            H(vrow,vcol) = 1.0;
             H(vrow,vcol) = Vmag(i);
        end
        Rcov(vrow, vrow) = Vmeassigma(i_vmeas)^2;
        if i == refbus
            refbus_Vrow_H = H_row_num;
        end
end
    

% Build the H Matrix rows for Angle Measurements
for i_ameas = 1:numAmeas
        i = Ameasbus(i_ameas);
        H_row_num = H_row_num + 1;
        arow = H_row_num;
        thcol = 2*i-1;
        if Ameasstatus(i_ameas) ~= 0        
           H(arow,thcol) = 1.0;
        end
        Rcov(arow, arow) = Ameassigma(i_ameas)^2;  
        if i == refbus
            refbus_Arow_H = H_row_num;
        end
end

% Build the H Matrix rows for Injection Measurements

for i_injmeas = 1: numImeas

        i = Imeasbus(i_injmeas);
        H_row_num = H_row_num + 1;
        prow = H_row_num;
        H_row_num = H_row_num + 1;
        qrow = H_row_num;
        
    if Imeasstatus(i_injmeas) ~= 0
                
    %   Build HPTh terms
        for j = 1: numbus
            %if j ~= refbus
                thcol = 2*j -1;
                vcol =  2*j;
                if i == j    
                  H(prow,thcol) = -Qinj(i) - B(i,i)*Vmag(i)^2; 
                else
                  H(prow,thcol) = Vmag(i)*Vmag(j)*(G(i,j)*sin(Theta(i)-Theta(j)) - B(i,j)*cos(Theta(i)-Theta(j)));
                end
           % end
        end % end build HPth terms loop

    %   Build HPV terms
        for j = 1: numbus
            %if j ~= refbus
                thcol = 2*j -1;
                vcol =  2*j;
                if i == j    
                  H(prow,vcol) = Pinj(i) + G(i,i)*Vmag(i)^2; 
                else
                  H(prow,vcol) = Vmag(i)*Vmag(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j)));
                end
            %end
        end % end build HPV term loop

    %   Build HQTh terms
        for j = 1: numbus
            %if j ~= refbus
                thcol = 2*j -1;
                vcol =  2*j;
                if i == j    
                  H(qrow,thcol) = Pinj(i) - G(i,i)*Vmag(i)^2;
                else
                  H(qrow,thcol)= -Vmag(i)*Vmag(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j)));
                end
            %end
         end %  end Build HQTh terms loop

    %   Build HQV  terms
        for j = 1: numbus
            %if j ~= refbus
                thcol = 2*j -1;
                vcol =  2*j;
                if i == j    
                  H(qrow,vcol) = Qinj(i) - B(i,i)*Vmag(i)^2;
                else
                  H(qrow,vcol) = Vmag(i)*Vmag(j)*(G(i,j)*sin(Theta(i)-Theta(j))-B(i,j)*cos(Theta(i)-Theta(j)));
                end
            %end
         end% end build HQV terms loop
         
    end %end Imeasstatus loop
    
          Rcov(prow, prow) = ImeasPsigma(i_injmeas)^2;
          Rcov(qrow, qrow) = ImeasQsigma(i_injmeas)^2;

      
end % end master bus loop building H matrix for Injection Measurements

% Build the H Matrix rows for Flow Measurements



for i_flowmeas = 1: numFmeas


        ibranch = Fmeasbranch(i_flowmeas);
        Gline =  R(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bline = -X(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bcapline = Bcap(ibranch)/2;

        i = Fmeasfrombus(i_flowmeas);
        j = Fmeastobus(i_flowmeas);

        H_row_num = H_row_num + 1;
        prow = H_row_num;
        H_row_num = H_row_num + 1;
        qrow = H_row_num;
        
    if Fmeasstatus(i_flowmeas) ~= 0
                
    %   Build H terms for from end (bus i)

        %if i ~= refbus
            thcol = 2*i -1;
            vcol =  2*i;    

            H(prow,thcol) = Gline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j)) - Bline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j));      % dP/dthi

            H(prow,vcol)  = (2*Gline*Vmag(i) - Gline*Vmag(j)*cos(Theta(i)-Theta(j)) - Bline*Vmag(j)*sin(Theta(i)-Theta(j)))*Vmag(i); % dP/dVi

            H(qrow,thcol) = - Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) - Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));    % dQ/dthi

            H(qrow,vcol)  =  (Bline*Vmag(j)*cos(Theta(i)-Theta(j)) - 2*Bcapline*Vmag(i) - 2*Bline*Vmag(i) - Gline*Vmag(j)*sin(Theta(i)-Theta(j)))*Vmag(i); % dQ/dVi
      %  end

    %   Build H terms for to end (bus j)

      %  if j ~= refbus
            thcol = 2*j -1;
            vcol =  2*j;    

            H(prow,thcol) =  Bline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) - Gline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j)); % dP/dthj

            H(prow,vcol)  =  (- Gline*Vmag(i)*cos(Theta(i)-Theta(j)) - Bline*Vmag(i)*sin(Theta(i)-Theta(j)))*Vmag(j);         % dP/dVj

            H(qrow,thcol) =  Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) + Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j)); % dQ/dthj

            H(qrow,vcol) =  (Bline*Vmag(i)*cos(Theta(i)-Theta(j)) - Gline*Vmag(i)*sin(Theta(i)-Theta(j)))*Vmag(j);            % dQ/dVj
       % end

    end %end Fmeasstatus loop
    
    Rcov(prow, prow) = ImeasPsigma(i_injmeas)^2;
    Rcov(qrow, qrow) = ImeasQsigma(i_injmeas)^2;
    
end % end master bus loop building H matrix for Flow Measurements


