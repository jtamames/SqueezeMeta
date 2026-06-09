**************
combineSQMlite
**************

============== ===============
combineSQMlite R Documentation
============== ===============

Combine several SQM or SQMlite objects
--------------------------------------

Description
~~~~~~~~~~~

Combine an arbitrary number of SQM or SQMlite objects into a single
SQMlite object. This function accepts objects originating from different
projects (i.e. different SqueezeMeta runs).

Usage
~~~~~

.. code:: R

   combineSQMlite(...)

Arguments
~~~~~~~~~

+---------+------------------------------------------------------------+
| ``...`` | an arbitrary number of SQM or SQMlite objects.             |
|         | Alternatively, a single list containing an arbitrary       |
|         | number of SQMlite objects.                                 |
+---------+------------------------------------------------------------+

Value
~~~~~

A SQMlite object

See Also
~~~~~~~~

``combineSQM``

Examples
~~~~~~~~

.. code:: R

   ## Not run: 
   data(Hadza)
   # Load data coming from a different run
   other = loadSQMlite("/path/to/other/project/tables") # e.g. if the project was run using sqm_reads
   # (We could also use loadSQM to load the data as long as the data comes from a SqueezeMeta run)
   combined = combineSQMlite(Hadza, other)
   # Now we can plot together the samples from Hadza and the second project
   plotTaxonomy(combined, 'family')

   ## End(Not run)
