function [damp] = damping(t,z,motionpath,m_d)

angle = motionpath.thet;

p = pi/2;

direction = m_d;

% ## Linear Extension ##

% c = 4429.133;
% c = 7086.614;
% c = 14173.228;
c = 17716.552;



while direction == 1 % deploying
    
    if z(1) > angle(332) 
        damp = c;
        return   
    elseif z(1) <= angle(332) 
        damp = 0;        
        return        
    end
    
end

while direction == -1 % retracting
    
    if z(1) > angle(332) 
        damp = 0;
        return
    elseif z(1) <= angle(332) 
        damp = c;
        return
    end

end

