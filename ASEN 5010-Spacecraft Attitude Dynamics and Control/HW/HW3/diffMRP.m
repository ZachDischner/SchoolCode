function sigmaprime = diffMRP(t,sigma)

    % Check for shadow-Schaub's methods are awesome
    sigma=MRPswitch(sigma,1);
    
    if isrow(sigma)
        sigma = sigma';
    end
    
    % Calculate [B(sigma)] Matrix
    w   = [1 0.5 -0.7]';
    B = BmatMRP(sigma);
%     sigmaprime = (1/4) * B * w;
    sigmaprime = dMRP(sigma,w); % Schaub's cool method!
    
end