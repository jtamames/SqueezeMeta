***********
exportTable
***********

.. container::

   =========== ===============
   exportTable R Documentation
   =========== ===============

   .. rubric:: Export results in tabular format
      :name: exportTable

   .. rubric:: Description
      :name: description

   This function is a wrapper for R's write.table function.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      exportTable(table, output_name)

   .. rubric:: Arguments
      :name: arguments

   +-----------------+---------------------------------------------------+
   | ``table``       | vector, matrix or data.frame. The table to be     |
   |                 | written.                                          |
   +-----------------+---------------------------------------------------+
   | ``output_name`` | Either a character string naming a file or a      |
   |                 | connection open for writing. ‘""’ indicates       |
   |                 | output to the console.                            |
   +-----------------+---------------------------------------------------+

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      Hadza.iron = subsetFun(Hadza, "iron")
      # Write the taxonomic distribution at the genus level of all the genes related to iron.
      exportTable(Hadza.iron$taxa$genus$percent, file.path(tempdir(), "Hadza.ironGenes.genus.tsv"))
      # Now write the distribution of the different iron-related COGs
      #  (Clusters of Orthologous Groups) across samples.
      exportTable(Hadza.iron$functions$COG$tpm, file.path(tempdir(), "Hadza.ironGenes.COG.tsv"))
      # Now write all the information contained in the ORF table.
      exportTable(Hadza.iron$orfs$table, file.path(tempdir(), "Hadza.ironGenes.orftable.tsv"))
