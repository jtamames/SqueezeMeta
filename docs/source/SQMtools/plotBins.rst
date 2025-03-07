********
plotBins
********

.. container::

   ======== ===============
   plotBins R Documentation
   ======== ===============

   .. rubric:: Barplot of the most abundant bins in a SQM object
      :name: plotBins

   .. rubric:: Description
      :name: description

   This function selects the most abundant bins across all samples in a
   SQM object and represents their abundances in a barplot.
   Alternatively, a custom set of bins can be represented.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      plotBins(
        SQM,
        count = "percent",
        N = 15,
        bins = NULL,
        others = TRUE,
        samples = NULL,
        ignore_unmapped = FALSE,
        ignore_nobin = FALSE,
        rescale = FALSE,
        color = NULL,
        base_size = 11,
        max_scale_value = NULL,
        metadata_groups = NULL
      )

   .. rubric:: Arguments
      :name: arguments

   +---------------------+-----------------------------------------------+
   | ``SQM``             | A SQM object.                                 |
   +---------------------+-----------------------------------------------+
   | ``count``           | character. Either ``"abund"`` for raw         |
   |                     | abundances, ``"percent"`` for percentages,    |
   |                     | ``"cov"`` for coverages, or ``"cpm"`` for     |
   |                     | coverages per million reads (default          |
   |                     | ``"percent"``).                               |
   +---------------------+-----------------------------------------------+
   | ``N``               | integer Plot the ``N`` most abundant bins     |
   |                     | (default ``15``).                             |
   +---------------------+-----------------------------------------------+
   | ``bins``            | character. Custom bins to plot. If provided,  |
   |                     | it will override ``N`` (default ``NULL``).    |
   +---------------------+-----------------------------------------------+
   | ``others``          | logical. Collapse the abundances of least     |
   |                     | abundant bins, and include the result in the  |
   |                     | plot (default ``TRUE``).                      |
   +---------------------+-----------------------------------------------+
   | ``samples``         | character. Character vector with the names of |
   |                     | the samples to include in the plot. Can also  |
   |                     | be used to plot the samples in a custom       |
   |                     | order. If not provided, all samples will be   |
   |                     | plotted (default ``NULL``).                   |
   +---------------------+-----------------------------------------------+
   | ``ignore_unmapped`` | logical. Don't include unmapped reads in the  |
   |                     | plot (default ``FALSE``).                     |
   +---------------------+-----------------------------------------------+
   | ``ignore_nobin``    | logical. Don't include reads which are not in |
   |                     | a bin in the plot (default ``FALSE``).        |
   +---------------------+-----------------------------------------------+
   | ``rescale``         | logical. Re-scale results to percentages      |
   |                     | (default ``FALSE``).                          |
   +---------------------+-----------------------------------------------+
   | ``color``           | Vector with custom colors for the different   |
   |                     | features. If empty, we will use our own       |
   |                     | hand-picked pallete if N<=15, and the default |
   |                     | ggplot2 palette otherwise (default ``NULL``). |
   +---------------------+-----------------------------------------------+
   | ``base_size``       | numeric. Base font size (default ``11``).     |
   +---------------------+-----------------------------------------------+
   | ``max_scale_value`` | numeric. Maximum value to include in the y    |
   |                     | axis. By default it is handled automatically  |
   |                     | by ggplot2 (default ``NULL``).                |
   +---------------------+-----------------------------------------------+
   | ``metadata_groups`` | list. Split the plot into groups defined by   |
   |                     | the user: list('G1' = c('sample1', sample2'), |
   |                     | 'G2' = c('sample3', 'sample4')) default       |
   |                     | ``NULL``).                                    |
   +---------------------+-----------------------------------------------+

   .. rubric:: Value
      :name: value

   a ggplot2 plot object.

   .. rubric:: See Also
      :name: see-also

   ``plotBins`` for plotting the most abundant bins of a SQM object;
   ``plotBars`` and ``plotHeatmap`` for plotting barplots or heatmaps
   with arbitrary data.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Bins distribution.
      plotBins(Hadza)
