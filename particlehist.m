function [ data ] = particlehist( k, mfp, absorbprob, x )
%A function that simulates the path taken by a total of k neutrons in a
%material of thickness x with a mean free path, 'mfp', and probability to 
% absorb neutrons, 'absorbprob'

absorb_count = 0;
reflect_count = 0;
transmit_count = 0;

% create a structure to append {x, y, x} data to
data = struct('x', 0, 'y', 0, 'z', 0);

for n = 1:k

    % create an initial position vector
    r = zeros(1, 3);
    % the particle enters from outside the shielding (x < 0)
    r(1, 1) = -1;

    % first step needs to be normally incident to shielding (x-direction)
    step = [1, 0, 0];
    % add step to the position vector history
    r = [r; r + step];
    
    % define flag variables for absoption, transmission and reflection
    % flag variables will be set to 1 if the process occurs during the
    % simulation
    is_absorbed = 0;
    is_transmit = 0;
    is_reflect = 0;
    
    % keep taking 'steps' whilst the particle hasn't been absorbed,
    % reflected or transmitted
    while is_absorbed == 0 || is_reflect == 0 || is_transmit == 0

        % generate a step in a random direction, of length sampled from the
        % exponential decay distribution
        step = unitvector(1).*exponentialdecay(1, mfp);
        % add step to particle history
        r = [r; r(end, :) + step];

        % if the particle is trasmitted (x > T)...
        if r(end, 1) > x

            % add one to trasmitted counter
            transmit_count = transmit_count + 1;
            % stop simulation (by setting flag variable to 1)
            is_transmit = 1;
            break

        % if the particle is reflected (x < T)...
        elseif r(end, 1) < 0

            % add one to reflected counter
            reflect_count = reflect_count + 1;
            % stop simulation (by setting flag variable to 1)
            is_reflect = 1;
            break

        else
            
            % particle is either absorbed scattered...
            % generate a random number
            u = rand();

            % if the random number is lower than the absorption probability
            % then the particle is absorbed...
            if u < absorbprob

                % add one to absorb counter
                absorb_count = absorb_count + 1;
                % stop simulation
                is_absorbed = 1;
                break

            end
            
            % if the particle is not absorbed at this point, the step
            % procudure (while loop) repeats
            
        end % matches second for loop
    end % matches while loop
    
        % add particle position history to a structure to plot particle path
    % add 'r' to structure
    [a, d] = size(r(:, 1));
    for m = 1:a
        data.x(m, n) = r(m, 1);
    end
    [b, d] = size(r(:, 2));
    for m = 1:b
        data.y(m, n) = r(m, 2);
    end
    [c, d] = size(r(:, 3));
    for m = 1:c
        data.z(m, n) = r(m, 3);
    end
    
end % matches main for loop
    
% simulation complete
    
end

