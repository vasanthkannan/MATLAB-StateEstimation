% start state estimation orthogonal iterative loop

%run subrouting to calculate number of active measurements
count_active_measurements

iteration = 1;
not_converged = 1;

while not_converged == 1

    % Build the H Matrix
    buildHmatrix

    % Calculate the mesairement Residual J and the measurement error vector
    calculate_residual_and_meas_error
    
    %*********************************************************
    %
    % old code:
    % DelTHV = inv( H'*inv(Rcov)*H)*H'*inv(Rcov)* meas_error;
    % 
    % max_delta = 0.0;
    %
    % Update voltages and phase angles
    % for i = 1: numbus
    %
    %    throw = 2*i-1;
    %
    %    Theta(i) = Theta(i) + DelTHV(throw);
    %    if abs( DelTHV(throw) ) > max_delta
    %        max_delta = abs( DelTHV(throw) );
    %    end
    %    vrow = 2*i;
    %    Vmag(i) = Vmag(i)*(1 + DelTHV(vrow));
    %    if abs( Vmag(i)* DelTHV(vrow) ) > max_delta
    %        max_delta = abs( Vmag(i)* DelTHV(vrow) );
    %    end
    %
    % end;
    %*********************************************************
    %
    % New code for EE8725 project
    % 11/8/2015
    %
    %calculate the new H matrix
    %[H'] = sqrt[R^(-1)] * [H]
    %TODO: create sparse

    H1 = sqrt(inv(Rcov))*H;
    
    %calculate the new R-covariance matrix
    %[Z'] = sqrt[R^(-1)] * Z
    %TODO: create sparse
    Z1 = sqrt(inv(Rcov))*meas_error;
    
    %This function calculates the Q matrix and U matrix of H'
    %TODO: create sparse 
    [Q,U] = qr(H1);
    
%     Q
%     U
    
    DelTHV = U\(Q' *Z1);
    
    max_delta = 0.0;
    
    % Update voltages and phase angles
    for i = 1: numbus
    
       throw = 2*i-1;
   
       Theta(i) = Theta(i) + DelTHV(throw);
       if abs( DelTHV(throw) ) > max_delta
           max_delta = abs( DelTHV(throw) );
       end
       vrow = 2*i;
       Vmag(i) = Vmag(i)*(1 + DelTHV(vrow));
       if abs( Vmag(i)* DelTHV(vrow) ) > max_delta
           max_delta = abs( Vmag(i)* DelTHV(vrow) );
       end
   
    end;
    %*********************************************************



    %fprintf('  % 5d     %12.3f          %15.8f\n', iteration, J, max_delta);
    fprintf('  % 5d     %12.3f          \n', iteration, J );
    iteration=iteration+1;
    %max_delta
    if iteration > 10 | max_delta <= Est_tolerance
        not_converged = 0;
    end
      
end % end of main while loop

