Tidal Heating
=================

Overview
--------

Tidal heating of Proxima b due to the CPL and CTL tidal evolutions presented in the paper.

===================   ============
**Date**              12/30/20
**Author**            Joseph Livesey
**Modules**           DistOrb
                      EqTide
                      STELLAR
**Approx. runtime**   5 minutes
===================   ============

To run this example
-------------------

.. code-block:: bash

    vplanet vpl.in
    python makeplot.py <pdf | png>

NB: These are the same simulations as presented in the "CPL_Model" and "CTL_Model" directories, with the addition of the "-SurfEnFluxEqtide" output argument (minus sign puts output in units of W/m^2).

Expected output
---------------

.. figure:: CPLTidalHeating-2deg.pdf
   :width: 150px
   :align: center

.. figure:: CPLTidalHeating-20deg.pdf
   :width: 150px
   :align: center

.. figure:: CTLTidalHeating-2deg.pdf
   :width: 150px
   :align: center

.. figure:: CTLTidalHeating-20deg.pdf
   :width: 150px
   :align: center

The figures show the semi-major axis, eccentricity, and tidal surface heat flux evolutions for Proxima b in each simulation.
