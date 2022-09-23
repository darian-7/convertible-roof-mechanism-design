function [value,isterminal,direction] = event_function2(t,z)

thet_in = pi/2 + 0.4186632; % rad - finish angle

value = z(1) - thet_in;
isterminal = 1;
direction = 0;

end
