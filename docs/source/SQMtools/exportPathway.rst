*************
exportPathway
*************

.. container::

   ============= ===============
   exportPathway R Documentation
   ============= ===============

   .. rubric:: Export the functions of a SQM object into KEGG pathway
      maps
      :name: exportPathway

   .. rubric:: Description
      :name: description

   This function is a wrapper for the pathview package (Luo *et al.*,
   2017. *Nucleic acids research*, 45:W501-W508). It will generate
   annotated KEGG pathway maps showing which reactions are present in
   the different samples. It will also generate legends with the color
   scales for each sample in separate png files.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      exportPathway(
        SQM,
        pathway_id,
        count = "copy_number",
        samples = NULL,
        split_samples = FALSE,
        sample_colors = NULL,
        log_scale = FALSE,
        fold_change_groups = NULL,
        fold_change_colors = NULL,
        max_scale_value = NULL,
        color_bins = 10,
        output_dir = ".",
        output_suffix = "pathview"
      )

   .. rubric:: Arguments
      :name: arguments

   +------------------------+--------------------------------------------+
   | ``SQM``                | A SQM, SQMbunch or SQMlite object.         |
   +------------------------+--------------------------------------------+
   | ``pathway_id``         | character. The five-number KEGG pathway    |
   |                        | identifier. A list of all pathway          |
   |                        | identifiers can be found in                |
   |                        | https://www.genome.jp/kegg/pathway.html.   |
   +------------------------+--------------------------------------------+
   | ``count``              | character. Either ``"abund"`` for raw      |
   |                        | abundances, ``"percent"`` for percentages, |
   |                        | ``"bases"`` for raw base counts, ``"tpm"`` |
   |                        | for TPM normalized values or               |
   |                        | ``"copy_number"`` for copy numbers         |
   |                        | (default ``"copy_number"``). Note that a   |
   |                        | given count type might not available in    |
   |                        | this object (e.g. TPM or copy number in    |
   |                        | SQMlite objects originating from a SQM     |
   |                        | reads project).                            |
   +------------------------+--------------------------------------------+
   | ``samples``            | character. An optional vector with the     |
   |                        | names of the samples to export. If absent, |
   |                        | all samples will be exported (default      |
   |                        | ``NULL``).                                 |
   +------------------------+--------------------------------------------+
   | ``split_samples``      | logical. Generate a different output file  |
   |                        | for each sample (default ``FALSE``).       |
   +------------------------+--------------------------------------------+
   | ``sample_colors``      | character. An optional vector with the     |
   |                        | plotting colors for each sample (default   |
   |                        | ``NULL``).                                 |
   +------------------------+--------------------------------------------+
   | ``log_scale``          | logical. Use a base 10 logarithmic         |
   |                        | transformation for the color scale. Will   |
   |                        | have no effect if ``fold_change_groups``   |
   |                        | is provided (default ``FALSE``).           |
   +------------------------+--------------------------------------------+
   | ``fold_change_groups`` | list. An optional list containing two      |
   |                        | vectors of samples. If provided, the       |
   |                        | function will generate a single plot       |
   |                        | displaying the log2 fold-change between    |
   |                        | the median abundances of both groups of    |
   |                        | samples ( log(second group / first group)  |
   |                        | ) (default ``NULL``).                      |
   +------------------------+--------------------------------------------+
   | ``fold_change_colors`` | character. An optional vector with the     |
   |                        | plotting colors of both groups in the      |
   |                        | fold-change plot. Will be ignored if       |
   |                        | ``fold_change_group`` is not provided.     |
   +------------------------+--------------------------------------------+
   | ``max_scale_value``    | numeric. Maximum value to include in the   |
   |                        | color scale. By default it is the maximum  |
   |                        | value in the selected samples (if plotting |
   |                        | abundances in samples) or the maximum      |
   |                        | absolute log2 fold-change (if plotting     |
   |                        | fold changes) (default ``NULL``).          |
   +------------------------+--------------------------------------------+
   | ``color_bins``         | numeric. Number of bins used to generate   |
   |                        | the gradient in the color scale (default   |
   |                        | ``10``).                                   |
   +------------------------+--------------------------------------------+
   | ``output_dir``         | character. Directory in which to write the |
   |                        | output files (default ``"."``).            |
   +------------------------+--------------------------------------------+
   | ``output_suffix``      | character. Suffix to be added to the       |
   |                        | output files (default ``"pathview"``).     |
   +------------------------+--------------------------------------------+

   .. rubric:: Value
      :name: value

   A ``ggplot`` if ``split_samples = FALSE`` and the
   `ggpattern <https://CRAN.R-project.org/package=ggpattern>`__ package
   is installed, otherwise nothing. Additionally, Pathview figures will
   be written in the directory specified by ``output_dir``.

   .. rubric:: See Also
      :name: see-also

   ``plotFunctions`` for plotting the most functions taxa of a SQM
   object.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)

      exportPathway(Hadza, "00910", count = 'copy_number',
                    output_dir = tempdir(),
                    output_suffix = "nitrogen_metabolism",
                    sample_colors = c("red", "blue"))
      exportPathway(Hadza, "00250", count = 'tpm',
                    output_dir = tempdir(),
                    output_suffix = "ala_asp_glu_metabolism_FoldChange", 
                    fold_change_groups = list(c("H1"), c("H12")), max_scale_value=2)
