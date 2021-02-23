from frb.associate import frbassociate
from frb.associate import frbs

from IPython import embed

"""
RUN THIS
"""
#frb_name = 'FRB121102'
#frb_name = 'FRB180916'
frb_name = 'FRB180924'
#frb_name = 'FRB190611'
#frb_name = 'FRB190614'
#frb_name = 'FRB181112'
#frb_name = 'FRB190102'
#frb_name = 'FRB190711'
#frb_name = 'FRB190714'
#frb_name = 'FRB191001'
#frb_name = 'FRB200430'

config = getattr(frbs, frb_name.lower())
config['skip_bayesian'] = True

frbA = frbassociate.run_individual(config, show=True, verbose=True)

embed()
