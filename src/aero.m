function [aero1] = aero(t,z,motionpath)

Z = 1;
tolerance = 1e-4;
size1 = height(motionpath);
angle = motionpath.thet;
A = motionpath.A;

while Z == 1
    for i = 1:size1
        if angle(i) > z(1) - tolerance && angle(i) < z(1) + tolerance
            aero1 = A(i);
            return
        end
    end
    tolerance = tolerance + tolerance/2;
end

end





