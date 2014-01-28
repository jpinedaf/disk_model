import math

def disk_model(x, y, P):
#
#  x and y are in [au]. The disk is inclinated on the plane of the sky around the x-axis
#  T_10 is temperature at 10au, 
#  q is the temperature power-law exponent
#  Mstar is the central star's mass in [Msun]
#  i_d is the inclination angle of the disk, with 90deg for face-on disk
#  Sigma_c is the column density at r_c, 
#  while r_c is the characteristic radius
#
#  sigma_NT is the non-thermal component of the linewidth, in m/s
#
  mu=2.37 # is the mean molecular weight 
  mu_co=28. # CO molecular weight 
  i_d=45. * math.pi/180.  # disk inclination angle (in the meantime)
  nu_0=115.e9
  sigma_NT=100. # in [m/s]
  Mh=1.660538921e-27 # Hydrogen mass in kg
  c = 2.99792458e8 #speed of light in [m/s]
  m_co=Mh*mu_co
  X_CO=1e-4  # abundance of CO/H_2 is fixed, because it is degenerate with Sigma_c
  M_CO=X_CO*Sigma_c*(2*math.pi * r_c^2)/(1-gamma)

  r=math.sqrt( x^2 + (y/math.cos(i_d))^2 )
  phi = math.atan2( y, x)

  Sigma=Sigma_c *(r/r_c)^(-gamma) * math.exp(-(r/r_c)^(2-gamma)) # Column density
  T = T_10 *(r/10.)^(-q)
  V = math.sqrt(G*Mstar/r) 
  Hp= math.sqrt(k*T/(mu*Mh) *(r^3/(G*Mstar)))
  rho=Sigma/( math.sqrt(2*pi) * Hp) * math.exp(-0.5*(z/Hp)^2)
  sigma_v=math.sqrt(sigma_NT^2 + 2*k*T/m_co)

  v_los=V*math.cos(phi)*sin(i_d) 
  return 0
