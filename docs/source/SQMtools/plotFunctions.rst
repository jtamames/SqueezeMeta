*************
plotFunctions
*************

.. container::

   ============= ===============
   plotFunctions R Documentation
   ============= ===============

   .. rubric:: Heatmap of the most abundant functions in a SQM object
      :name: plotFunctions

   .. rubric:: Description
      :name: description

   This function selects the most abundant functions across all samples
   in a SQM object and represents their abundances in a heatmap.
   Alternatively, a custom set of functions can be represented.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      plotFunctions(
        SQM,
        fun_level = "KEGG",
        count = "copy_number",
        N = 25,
        fun = NULL,
        samples = NULL,
        display_function_names = TRUE,
        ignore_unmapped = TRUE,
        ignore_unclassified = TRUE,
        gradient_col = c("ghostwhite", "dodgerblue4"),
        rescale_percent = FALSE,
        base_size = 11,
        metadata_groups = NULL
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------+----------------------------------------+
   | ``SQM``                    | A SQM, SQMbunch or SQMlite object.     |
   +----------------------------+----------------------------------------+
   | ``fun_level``              | character. Either ``"KEGG"``,          |
   |                            | ``"COG"``, ``"PFAM"`` or any other     |
   |                            | custom database used for annotation    |
   |                            | (default ``"KEGG"``).                  |
   +----------------------------+----------------------------------------+
   | ``count``                  | character. Either ``"abund"`` for raw  |
   |                            | abundances, ``"percent"`` for          |
   |                            | percentages, ``"bases"`` for raw base  |
   |                            | counts, ``"cpm"`` for coverages per    |
   |                            | million reads, ``"tpm"`` for TPM       |
   |                            | normalized values or ``"copy_number"`` |
   |                            | for copy numbers (default              |
   |                            | ``"copy_number"``). Note that a given  |
   |                            | count type might not available in this |
   |                            | object (e.g. TPM or copy number in     |
   |                            | SQMlite objects originating from a SQM |
   |                            | reads project).                        |
   +----------------------------+----------------------------------------+
   | ``N``                      | integer Plot the ``N`` most abundant   |
   |                            | functions (default ``25``).            |
   +----------------------------+----------------------------------------+
   | ``fun``                    | character. Custom functions to plot.   |
   |                            | If provided, it will override ``N``    |
   |                            | (default ``NULL``).                    |
   +----------------------------+----------------------------------------+
   | ``samples``                | character. Character vector with the   |
   |                            | names of the samples to include in the |
   |                            | plot. Can also be used to plot the     |
   |                            | samples in a custom order. If not      |
   |                            | provided, all samples will be plotted  |
   |                            | (default ``NULL``).                    |
   +----------------------------+----------------------------------------+
   | ``display_function_names`` | logical. Plot function names alongside |
   |                            | function IDs, if available (default    |
   |                            | ``TRUE``).                             |
   +----------------------------+----------------------------------------+
   | ``ignore_unmapped``        | logical. Don't include unmapped reads  |
   |                            | in the plot (default ``TRUE``).        |
   +----------------------------+----------------------------------------+
   | ``ignore_unclassified``    | logical. Don't include unclassified    |
   |                            | ORFs in the plot (default ``TRUE``).   |
   +----------------------------+----------------------------------------+
   | ``gradient_col``           | A vector of two colors representing    |
   |                            | the low and high ends of the color     |
   |                            | gradient (default                      |
   |                            | ``c("ghostwhite", "dodgerblue4")``).   |
   +----------------------------+----------------------------------------+
   | ``rescale_percent``        | logical. Calculate percent counts over |
   |                            | the number of reads in the input       |
   |                            | object, instead of over the total      |
   |                            | number of reads in the original        |
   |                            | project (default ``FALSE``).           |
   +----------------------------+----------------------------------------+
   | ``base_size``              | numeric. Base font size (default       |
   |                            | ``11``).                               |
   +----------------------------+----------------------------------------+
   | ``metadata_groups``        | list. Split the plot into groups       |
   |                            | defined by the user: list('G1' =       |
   |                            | c('sample1', sample2'), 'G2' =         |
   |                            | c('sample3', 'sample4')) default       |
   |                            | ``NULL``).                             |
   +----------------------------+----------------------------------------+

   .. rubric:: Value
      :name: value

   a ggplot2 plot object.

   .. rubric:: See Also
      :name: see-also

   ``plotTaxonomy`` for plotting the most abundant taxa of a SQM object;
   ``plotBars`` and ``plotHeatmap`` for plotting barplots or heatmaps
   with arbitrary data.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      plotFunctions(Hadza)
