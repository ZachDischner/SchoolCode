function Pd = delayPRN(P,n)
    Pd = [P(end-n+1:end),P(1:end-n)];
end