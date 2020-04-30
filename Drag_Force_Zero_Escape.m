% Acoustophoretic Volumetric Display - Ultimate Capability (AVD-UC).
% 
% Codes for: "What is the Ultimate Capability of Acoustophoretic Volumetric
% Displays?".
% 
% Authors: Tatsuki Fushimi, Bruce W. Drinkwater, and Thomas L. Hill.
% 
% Journal: Applied Physics Letters.
% 
% Created: 30-Apr-2020.
% 
% Please contact Tatsuki FUSHIMI (t.fushimi@bristol.ac.uk) for any inquiry.
% Released under MIT License

function [outputArg1] = Drag_Force_Zero_Escape(eq_type,r,rho_0,vel_z,visc)

switch eq_type
    case 1 %Single Direction vel_z >= 0
        outputArg1 = F_drag(r,rho_0,vel_z,visc);
        outputArg1(vel_z==0) = 0;
    case 2 %Bidirectional any real vel_z
        outputArg1 = F_drag_abs(r,rho_0,vel_z,visc);
        outputArg1(vel_z==0) = 0;
    otherwise
        warning('Invalid Equation type')
end

end

