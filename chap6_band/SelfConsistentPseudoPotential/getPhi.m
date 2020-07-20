function Phi_r = getPhi(V, K, r)
    KDotR = K * r;
    structureFactor = exp(1j*KDotR);
    Phi_r = V' * structureFactor;
end