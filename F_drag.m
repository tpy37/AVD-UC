function F_drag = F_drag(r,rho_0,vel_z,visc)
%F_DRAG
%    F_DRAG = F_DRAG(R,RHO_0,VEL_Z,VISC)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    19-Apr-2020 19:34:28

t2 = visc.^2;
t3 = t2.^2;
t4 = 1.0./visc;
F_drag = r.^2.*rho_0.*vel_z.^2.*pi.*((1.0./(r.*rho_0.*t4.*vel_z.*7.604562737642586e-6).^(3.97e2./5.0e1).*(4.11e2./1.0e3))./(1.0./r.^8.*1.0./rho_0.^8.*t3.^2.*1.0./vel_z.^8.*8.941410269742584e40+1.0)+(visc.*1.2e1)./(r.*rho_0.*vel_z)+(r.*rho_0.*t4.*vel_z.*5.0e-7)./(r.*rho_0.*t4.*vel_z.*2.0e-6+1.0)+(r.*rho_0.*t4.*vel_z.*(2.6e1./2.5e1))./((r.*rho_0.*t4.*vel_z.*(2.0./5.0)).^(3.8e1./2.5e1)+1.0)).*(-1.0./2.0);
