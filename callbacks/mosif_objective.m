function f = mosif_objective(x, auxdata)
if auxdata.isNLP
    f = feval(auxdata.funcs.obj, x, auxdata);
else
    f = x' * auxdata.Q * x + auxdata.q' * x + auxdata.p;
end
end

