

****************
Offset Functions
****************

This doc describes the offset functions that are coded
in astropath.  A future refactor will allow for user-supplied
functions.

In our formalism, the offset function is written as
:math:`p(\omega|O)`
and it describes the *true* distribution of the transient within
its host galaxy.

Definitions
===========

Here are the key items related to :math:`p(\omega|O)`:

* :math:`\omega`: 3-D true position of the transient
* :math:`O`: galaxy properties (e.g. position, angular size)
* :math:`\theta`: offset between the
  transient and the center of a galaxy.
* :math:`\phi`: galaxy angular size
* :math:`\theta_{max}`: the maximum separation between
  the transient and the galaxy allowed in units of
  :math:`\phi`.  :math:`p(\omega|O) = 0` for
  :math:`\theta > \theta_{max}`

theta_prior
===========

The code base uses a *dict* named *theta_prior* to handle the
offset function.  It currently holds three items:

* method -- A string defining the method used.  core, uniform, exp
* max -- :math:`{\theta_{max}}`
  in units of :math:`\phi`.  :math:`p(\omega|O) = 0` for
  :math:`\theta > \theta_{max}`
* ang_size -- An array of angular sizes :math:`\phi`
  for all of the candidate galaxies.

Below are the 3 offset functions currently coded.
Note that each of these are normalized to unit total
probabilty when integrated over the full sphere.

uniform
+++++++

Here, :math:`p(\omega|O)` is given equal weighting for all locations with
:math:`\theta < \theta_{max}`.

core
++++

Here, :math:`p(\omega|O)` is proportional to :math:`\phi / (\theta + \phi)`

exp
+++

Here, :math:`p(\omega|O)` is proportional to
:math:`(\theta/\phi) \, \exp [-(\theta/\phi)]`

