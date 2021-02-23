from astropy.cosmology import Planck15 as cosmo

nhalf = 10.

# Theta priors
theta_max = 6.
theta_u = dict(method='uniform', max=theta_max)
theta_c = dict(method='core', max=theta_max)
theta_e = dict(method='exp', max=theta_max)

# Combined priors
conservative = dict(theta=theta_u, O='identical', U=0, name='Conservative', nhalf=nhalf)
adopted = dict(theta=theta_e, O='inverse', U=0., name='Adopted', nhalf=nhalf)
#extreme = dict(theta=theta_c, O='inverse', U=-1, name='Extreme', nhalf=nhalf)

# FRB list
frb_list = ['FRB121102', 'FRB180916', 'FRB180924', 'FRB181112',
            'FRB190102', 'FRB190523', 'FRB190608', 'FRB190611', 'FRB190614',
            'FRB190711', 'FRB190714', 'FRB191001', 'FRB200430']

# Secure P(O|x)
POx_secure = 0.95
