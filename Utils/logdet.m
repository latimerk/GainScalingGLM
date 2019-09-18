function logdetA = logdet(G)

try 
    G_det = chol(G);
    logdetA = 2*sum(log(diag(G_det)));
catch
    logdetA = log(max(1e-20,det(G)));
end