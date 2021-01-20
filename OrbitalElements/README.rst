Orbital elements of Proxima b
=================

===================   ============
**Date**              1/19/21
**Author**            Joseph Livesey
**Approx. runtime**   8 CPU hours (simulations)
                      < 90 minutes (plotting)
===================   ============

The orbital plane of Proxima b is presently unconstrained. In the paper, we claim that the initial longitude of ascending node (Ω) and longitude of pericenter (ϖ) of Proxima b have no bearing on our major result: the timescale on which Proxima b is captured into the synchronous rotation state. Here, we present 100 evolutions of Proxima b's rotational period, assembled with [VSPACE](https://github.com/VitualPlanetaryLaboratory/vplanet/tree/master/vspace) using the CPL tidal model and an initial mutual inclination of 2º between Proxima b and Proxima c.

To execute:

.. code-block:: bash
    vspace vspace.in

Expected output:

.. figure:: initialspace.png
   :width: 600px
   :align: center

The evolutions are virtually identical under our model. Therefore, we can comfortably set both parameters to zero initially.
