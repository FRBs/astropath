****************
Offset Functions
****************

This doc describes the offset functions that are coded
in astropath.  A future refactor will allow for user-supplied
functions.

In our formalism, the offset function is written as
:math:`p(w|O)`
and it describes the *true* distribution of the transient within
its host galaxy.

Throughout, :math:`{\theta}` refers to the offset between the
transient and the center of a galaxy.

theta_prior
===========

The code base uses a *dict* named *theta_prior* to handle the
offset function.  It currently holds three items:

* method -- A string defining the method used.  core, uniform, exp
* max -- The maximum separation :math:`{\theta_{max}}`
  in units of angular size.  :math:`p(w|O) = 0` for
  :math:`\theta > \theta_{max}`
* ang_size -- An array of angular sizes :math:`\phi`
  for all of the galaxies.

Here are the 3 offset functions currently coded:

uniform
+++++++

Here, :math:`p(w|O)` is given equal weighting for all locations with
:math:`\theta < \theta_{max}`.

core
++++

Here, :math:`p(w|O)` is proportional to :math:`\phi / (\theta + \phi)`

exp
+++

Here, :math:`p(w|O)` is proportional to
:math:`(\theta/\phi) \, \exp [-(\theta/\phi)]`

