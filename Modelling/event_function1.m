function [value,isterminal,direction] = event_function1(t,z)

if z(2) >= 0
    
    value = z(1) - (pi/2 + 2.6213734);
    
else
    
    value = z(1) - (pi/2 + 0.4186630);
    
end

isterminal = 1;
direction = 0;

end