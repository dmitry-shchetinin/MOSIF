function df = mosif_gradient(x, auxdata)
    if auxdata.isNLP
        df = feval(auxdata.funcs.grad, x, auxdata);
    else
        df = 2 * auxdata.Q * x + auxdata.q;
    end
end