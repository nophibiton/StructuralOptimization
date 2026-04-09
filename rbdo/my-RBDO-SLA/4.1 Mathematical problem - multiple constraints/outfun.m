function stop = outfun(dk, optimValues, state)

    global xk_prev xk_buffer last_dk_evaluated

    stop = false;

    if strcmp(state,'iter')

        fmt = ['Iteration: %1d, ' 'x = [%.4f' repmat(' %.4f', 1, numel(dk)-1) '].\n'];
        fprintf(fmt,optimValues.iteration,dk);
        % disp(['Iteration ', num2str(optimValues.iteration), ...
        %       ', x = ', num2str(dk)]);
        
        xk_prev = xk_buffer;

        % reset tracker for next iteration
        last_dk_evaluated = [];
    end
end