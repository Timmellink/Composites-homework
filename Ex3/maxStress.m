function fail = maxStress(sig1,sig2,sig6,s1c,s1t,s2c,s2t,s6)
if -s1c<sig1 && sig1<s1t
    disp("No failure in 1 direction");
    if -s2c<sig2 && sig2<s2t
        disp("No failure in 2 direction")
        if sig6<s6
            disp("no failure in 6 direction");
            fail = false;
        else
            disp("Failure in 6 direction");
            fail = true;
        end
    else
        disp("failure in 2 direction");
        fail = true;
    end
else
    disp("failure in 1 direction");
    fail = true;
end
end

