% start state est iterative loop


%run subrouting to calculate number of active measurements
count_active_measurements

iteration = 1;
not_converged = 1;

while not_converged == 1

    % Build the H Matrix
    buildHmatrix

    % Calculate the mesairement Residual J and the measurement error vector
    calculate_residual_and_meas_error
    

    DelTHV = inv( H'*inv(Rcov)*H)*H'*inv(Rcov)* meas_error;

%     Htest =  inv( H'*inv(Rcov)*H)*H'*inv(Rcov);

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


    %fprintf('  % 5d     %12.3f          %15.8f\n', iteration, J, max_delta);
    fprintf('  % 5d     %12.3f          \n', iteration, J );
    iteration=iteration+1;
    %max_delta
    if iteration > 10 | max_delta <= Est_tolerance
        not_converged = 0;
    end
      
end % end of main while loop

