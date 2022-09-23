function [radius1] = radius(t,z,motionpath)

Z = 1;
tolerance = 1e-4;
size1 = height(motionpath);
angle = motionpath.thet;
rad = motionpath.r;

    while Z == 1
        for i = 1:size1
            if angle(i) > z(1) - tolerance && angle(i) < z(1) + tolerance
                radius1 = rad(i)*0.001;
                return
            end
        end
        tolerance = tolerance + tolerance/2;
    end

end
