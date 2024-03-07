function res = biquadratic(x1, x2, coeff)
    res = coeff(1) + (coeff(2)+coeff(3)*x1).*x1 + (coeff(4)+coeff(5)*x2).*x2 + coeff(6)*x1.*x2;
end