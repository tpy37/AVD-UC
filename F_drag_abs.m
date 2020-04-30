function F_drag_abs = F_drag_abs(r,rho_0,vel_z,visc)
%F_DRAG_ABS
%    F_DRAG_ABS = F_DRAG_ABS(R,RHO_0,VEL_Z,VISC)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    19-Apr-2020 19:34:28

t2 = r.*rho_0.*vel_z;
t3 = abs(t2);
t4 = abs(visc);
t5 = t4.^2;
t6 = t5.^2;
t7 = 1.0./t4;
F_drag_abs = r.^2.*rho_0.*vel_z.^2.*pi.*sign(vel_z).*((1.0./(t3.*t7.*7.604562737642586e-6).^(3.97e2./5.0e1).*(4.11e2./1.0e3))./(1.0./t3.^8.*t6.^2.*8.941410269742584e40+1.0)+(t4.*1.2e1)./t3+(t3.*t7.*5.0e-7)./(t3.*t7.*2.0e-6+1.0)+(t3.*t7.*(2.6e1./2.5e1))./((t3.*t7.*(2.0./5.0)).^(3.8e1./2.5e1)+1.0)).*(-1.0./2.0);