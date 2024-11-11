function fail = maxStrain(eps1,eps2,eps6,e1c,e1t,e2c,e2t,e6)
if -e1c<eps1 && eps1<e1t
    disp("No failure in 1 direction");
    if -e2c<eps2 && eps2<e2t
        disp("No failure in 2 direction");
        if eps6<e6
            disp("no failure in 6 direction");
            fail = false;
        else
            disp("Failure in 6 direction");
            fail = true;
        end
    else 
        disp("Failure in 2 direction");
        fail = true;
    end
else
    disp("failure in 1 direction");
    fail = true;
end

end
