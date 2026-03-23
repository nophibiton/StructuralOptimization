function stop = outfun(dk, optimValues, state)

    global xk_prev xk_buffer last_dk_evaluated

    stop = false;

    if strcmp(state,'iter')

        disp(['Iteration ', num2str(optimValues.iteration), ...
              ', x = ', num2str(dk)]);
        
        xk_prev = xk_buffer;

        % reset tracker for next iteration
        last_dk_evaluated = [];
    end
end