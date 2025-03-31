***********
exportKrona
***********

.. container::

   =========== ===============
   exportKrona R Documentation
   =========== ===============

   .. rubric:: Export the taxonomy of a SQM object into a Krona Chart
      :name: exportKrona

   .. rubric:: Description
      :name: description

   Generate a krona chart containing the full taxonomy from a SQM
   object.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      exportKrona(SQM, output_name = NA)

   .. rubric:: Arguments
      :name: arguments

   +-----------------+---------------------------------------------------+
   | ``SQM``         | A SQM, SQMbunch or SQMlite object.                |
   +-----------------+---------------------------------------------------+
   | ``output_name`` | character. Name of the output file containing the |
   |                 | Krona charts in html format (default              |
   |                 | ``"<project_name>.krona.html")``.                 |
   +-----------------+---------------------------------------------------+

   .. rubric:: Details
      :name: details

   Original code was kindly provided by Giuseppe D'Auria
   (dauria_giu@gva.es).

   .. rubric:: Value
      :name: value

   No return value, but a krona chart is produced in the current working
   directory.

   .. rubric:: See Also
      :name: see-also

   ``plotTaxonomy`` for plotting the most abundant taxa of a SQM object.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Check that kronatools is present.
      ecode = system('ktImportText', ignore.stdout = TRUE, ignore.stderr = TRUE)
      # If so, run.
      if(ecode==0) { exportKrona(Hadza, output_name = file.path(tempdir(), "krona.html")) }
