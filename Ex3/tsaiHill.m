function fail = tsaiHill(sig1,sig2,sig6,s1c,s1t,s2c,s2t,s6)
if sig1<0 && sig2<0 % compressive cases
    val = sig1^2/s1c^2-sig1*sig2/s1c^2+sig2^2/s2c^2+sig6^2/s6^2;
    if val<1
        disp("No failure");
        fail = false;
    else
        disp("Failure");
        fail = true;
    end
elseif sig1>0 && sig2>0 % tensile cases
    val = sig1^2/s1t^2-sig1*sig2/s1t^2+sig2^2/s2t^2+sig6^2/s6^2;
    if val<1
        disp("No failure");
        fail = false;
    else
        disp("Failure");
        fail = true;
    end
elseif sig1>0 && sig2<0 % 1 tensile, 2 compressive
    val = sig1^2/s1t^2-sig1*sig2/s1t^2+sig2^2/s2c^2+sig6^2/s6^2;
    if val<1
        disp("No failure");
        fail = false;
    else
        disp("Failure");
        fail = true;
    end
elseif sig1<0 && sig2>0 % 1 compressive, 2 tensile
    val = sig1^2/s1c^2-sig1*sig2/s1c^2+sig2^2/s2t^2+sig6^2/s6^2;
    if val<1
        disp("No failure");
        fail = false;
    else
        disp("Failure");
        fail = true;
    end
end
end