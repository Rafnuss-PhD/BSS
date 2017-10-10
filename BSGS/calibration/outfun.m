function stop = outfun(x, optimValues, state)
stop = false;
switch state
    case 'init'
        header = ' Iteration   Func-count     min f(x)      x        Procedure';
        disp(header)
    case 'iter'
        fprintf(' %5.0f        %5.0f     %12.6g    %12.6g    %s\n', optimValues.iteration, optimValues.funccount, optimValues.fval, x, optimValues.procedure);
    case 'done'
        fprintf(' %5.0f        %5.0f     %12.6g    %12.6g    %s\n', optimValues.iteration, optimValues.funccount, optimValues.fval, x, optimValues.procedure);
end


