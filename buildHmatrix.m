% Store ones on the diagonal of the H Matrix

% The H Matrix is built with one Vrow for each Voltage meas and 
% with one Arow for each Angle meas 

% The H Matrix is built with Prows(prow) and Qrows (qrow) for the
% injection measurements and flow measurements 

% All rows have Theta columns (thcol) and V columns (vcol). 
num_rows_H = num_active_measurements;
num_cols_H = 2*numbus;
clear H Rcov
H = zeros( num_rows_H, num_cols_H );
Rcov = zeros( num_rows_H, num_rows_H );

H_row_num = 0;

% H
% Rcov
% size(Rcov)
% Build the H Matrix rows for Voltage Measurements
for i_vmeas = 1:numVmeas
    i = Vmeasbus(i_vmeas);
    if i ~= refbus
        if Vmeasstatus(i_vmeas) ~= 0     
            H_row_num = H_row_num + 1;
            vrow = H_row_num;
            vcol = 2*i;

            H(vrow,vcol) = 1.0;
%           H(vrow,vcol) = Vmag(i);

            Rcov(vrow, vrow) = Vmeassigma(i_vmeas)^2;
        end
            
     elseif i == refbus
            
            H_row_num = H_row_num + 1;
            vrow = H_row_num;
            vcol = 2*i;
            H(vrow,vcol) = 1.0;
            refbus_Vrow_H = H_row_num;
            Rcov(vrow, vrow) = 1.0;
            
     end
end
% H
% Rcov
% size(Rcov)
% Build the H Matrix rows for Angle Measurements
for i_ameas = 1:numAmeas
    i = Ameasbus(i_ameas);
    if i ~= refbus
        if Ameasstatus(i_ameas) ~= 0     
            H_row_num = H_row_num + 1;
            arow = H_row_num;
            thcol = 2*i-1;

            H(arow,thcol) = 1.0;


            Rcov(arow, arow) = Ameassigma(i_ameas)^2;
        end
            
    elseif i == refbus
            
            H_row_num = H_row_num + 1;
            arow = H_row_num;
            thcol = 2*i-1;
            H(arow,thcol) = 1.0;
            refbus_Arow_H = H_row_num;
            Rcov(arow, arow) = 1.0;
            
    end
end
% H    
% Rcov
% size(Rcov)

% Build the H Matrix rows for Injection Measurements
% skip the inj at the refbus

for i_injmeas = 1: numImeas
  i = Imeasbus(i_injmeas);
  if i ~= refbus
      
    if Imeasstatus(i_injmeas) ~= 0
              
        H_row_num = H_row_num + 1;
        prow = H_row_num;
        H_row_num = H_row_num + 1;
        qrow = H_row_num;
        
                
    %   Build HPTh terms
        for j = 1: numbus
            thcol = 2*j -1;
            vcol =  2*j;
            if j ~= refbus
                if i == j    
                  H(prow,thcol) = -Qinj(i) - B(i,i)*Vmag(i)^2; 
                else
                  H(prow,thcol) = Vmag(i)*Vmag(j)*(G(i,j)*sin(Theta(i)-Theta(j)) - B(i,j)*cos(Theta(i)-Theta(j)));
                end
            elseif j == refbus
                  H(prow,thcol) = 0.0;
                  H(prow,thcol) = 0.0;
            end
        end % end build HPth terms loop

    %   Build HPV terms
        for j = 1: numbus
            thcol = 2*j -1;
            vcol =  2*j;
            if j ~= refbus
                if i == j    
                  H(prow,vcol) = Pinj(i) + G(i,i)*Vmag(i)^2; 
                else
                  H(prow,vcol) = Vmag(i)*Vmag(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j)));
                end
            elseif j == refbus
                  H(prow,vcol) = 0.0;
                  H(prow,vcol) = 0.0;                
            end
        end % end build HPV term loop

    %   Build HQTh terms
        for j = 1: numbus
            thcol = 2*j -1;
            vcol =  2*j;
            if j ~= refbus
                if i == j    
                  H(qrow,thcol) = Pinj(i) - G(i,i)*Vmag(i)^2;
                else
                  H(qrow,thcol)= -Vmag(i)*Vmag(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j)));
                end
            elseif j == refbus
                  H(qrow,thcol) = 0.0;
                  H(qrow,thcol) = 0.0;                
            end
         end %  end Build HQTh terms loop

    %   Build HQV  terms
        for j = 1: numbus
            thcol = 2*j -1;
            vcol =  2*j;
            if j ~= refbus
                if i == j    
                  H(qrow,vcol) = Qinj(i) - B(i,i)*Vmag(i)^2;
                else
                  H(qrow,vcol) = Vmag(i)*Vmag(j)*(G(i,j)*sin(Theta(i)-Theta(j))-B(i,j)*cos(Theta(i)-Theta(j)));
                end
            elseif j == refbus
                  H(qrow,vcol) = 0.0;
                  H(qrow,vcol) = 0.0;                
            end
         end% end build HQV terms loop

          Rcov(prow, prow) = ImeasPsigma(i_injmeas)^2;
          Rcov(qrow, qrow) = ImeasQsigma(i_injmeas)^2;
          
    end   %end Imeasstatus test       
  end     % end refbus test

      
end % end master bus loop building H matrix for Injection Measurements
% H
% Rcov
% size(Rcov)

% Build the H Matrix rows for Flow Measurements

for i_flowmeas = 1: numFmeas
   if Fmeasstatus(i_flowmeas) ~= 0

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
 

                
    %   Build H terms for from end (bus i)
        thcol = 2*i -1;
        vcol =  2*i;    
        if i ~= refbus


            H(prow,thcol) = Gline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j)) - Bline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j));             % dP/dthi

            H(prow,vcol)  = (2*Gline*Vmag(i) - Gline*Vmag(j)*cos(Theta(i)-Theta(j)) - Bline*Vmag(j)*sin(Theta(i)-Theta(j)))*Vmag(i); % dP/dVi

            H(qrow,thcol) = - Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) - Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));           % dQ/dthi

            H(qrow,vcol)  =  (Bline*Vmag(j)*cos(Theta(i)-Theta(j)) - 2*Bcapline*Vmag(i) - 2*Bline*Vmag(i) - Gline*Vmag(j)*sin(Theta(i)-Theta(j)))*Vmag(i); % dQ/dVi
            
        elseif i == refbus
          
             H(prow,thcol) = 0.0;
             H(prow,vcol)  = 0.0;
             H(qrow,thcol) = 0.0;
             H(qrow,vcol)  = 0.0; 
        end

    %   Build H terms for to end (bus j)

        thcol = 2*j -1;
        vcol =  2*j;    
        if j ~= refbus


            H(prow,thcol) =  Bline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) - Gline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));     % dP/dthj

            H(prow,vcol)  =  (- Gline*Vmag(i)*cos(Theta(i)-Theta(j)) - Bline*Vmag(i)*sin(Theta(i)-Theta(j)))*Vmag(j);         % dP/dVj

            H(qrow,thcol) =  Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j)) + Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));     % dQ/dthj

            H(qrow,vcol) =  (Bline*Vmag(i)*cos(Theta(i)-Theta(j)) - Gline*Vmag(i)*sin(Theta(i)-Theta(j)))*Vmag(j);            % dQ/dVj
       
        elseif j == refbus
            
            H(prow,thcol) = 0.0;
            H(prow,vcol)  = 0.0;
            H(qrow,thcol) = 0.0;
            H(qrow,vcol)  = 0.0;
        end

    Rcov(prow, prow) = FmeasPsigma(i_flowmeas)^2;
    Rcov(qrow, qrow) = FmeasQsigma(i_flowmeas)^2;
    
    end %end Fmeasstatus loop
    
end % end master bus loop building H matrix for Flow Measurements

% H
% Rcov
% size(Rcov)
