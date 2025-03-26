function x = polar_transform(u)
    if(length(u)==1)
        x = u;
    else
        u1u2 = mod(u(1:2:end)+u(2:2:end),2);
        u2 = u(2:2:end);

        x = [polar_transform(u1u2) polar_transform(u2)];
    end