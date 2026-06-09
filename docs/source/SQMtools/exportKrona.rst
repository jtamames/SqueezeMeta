***********
exportKrona
***********

=========== ===============
exportKrona R Documentation
=========== ===============

Export the taxonomy of a SQM object into a Krona Chart
------------------------------------------------------

Description
~~~~~~~~~~~

Generate a krona chart containing the full taxonomy from a SQM object.

Usage
~~~~~

.. code:: R

   exportKrona(SQM, output_name = NA)

Arguments
~~~~~~~~~

+-----------------+----------------------------------------------------+
| ``SQM``         | A SQM, SQMbunch or SQMlite object.                 |
+-----------------+----------------------------------------------------+
| ``output_name`` | character. Name of the output file containing the  |
|                 | Krona charts in html format (default               |
|                 | ``"<project_name>.krona.html")``.                  |
+-----------------+----------------------------------------------------+

Details
~~~~~~~

Original code was kindly provided by Giuseppe D'Auria
(dauria_giu@gva.es).

Value
~~~~~

No return value, but a krona chart is produced in the current working
directory.

See Also
~~~~~~~~

``plotTaxonomy`` for plotting the most abundant taxa of a SQM object.

Examples
~~~~~~~~

.. code:: R


   data(Hadza)
   # Check that kronatools is present.
   ecode = system('ktImportText', ignore.stdout = TRUE, ignore.stderr = TRUE)
   # If so, run.
   if(ecode==0) { exportKrona(Hadza, output_name = file.path(tempdir(), "krona.html")) }
