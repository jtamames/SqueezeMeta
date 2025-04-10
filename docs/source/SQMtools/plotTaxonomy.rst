************
plotTaxonomy
************

.. container::

   ============ ===============
   plotTaxonomy R Documentation
   ============ ===============

   .. rubric:: Barplot of the most abundant taxa in a SQM object
      :name: plotTaxonomy

   .. rubric:: Description
      :name: description

   This function selects the most abundant taxa across all samples in a
   SQM object and represents their abundances in a barplot.
   Alternatively, a custom set of taxa can be represented.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      plotTaxonomy(
        SQM,
        rank = "phylum",
        count = "percent",
        N = 15,
        tax = NULL,
        others = TRUE,
        samples = NULL,
        nocds = "treat_separately",
        ignore_unmapped = FALSE,
        ignore_unclassified = FALSE,
        no_partial_classifications = FALSE,
        rescale = FALSE,
        color = NULL,
        base_size = 11,
        max_scale_value = NULL,
        metadata_groups = NULL
      )

   .. rubric:: Arguments
      :name: arguments

   +--------------------------------+------------------------------------+
   | ``SQM``                        | A SQM, SQMbunch or a SQMlite       |
   |                                | object.                            |
   +--------------------------------+------------------------------------+
   | ``rank``                       | Taxonomic rank to plot (default    |
   |                                | ``phylum``).                       |
   +--------------------------------+------------------------------------+
   | ``count``                      | character. Either ``"percent"``    |
   |                                | for percentages, or ``"abund"``    |
   |                                | for raw abundances (default        |
   |                                | ``"percent"``).                    |
   +--------------------------------+------------------------------------+
   | ``N``                          | integer Plot the ``N`` most        |
   |                                | abundant taxa (default ``15``).    |
   +--------------------------------+------------------------------------+
   | ``tax``                        | character. Custom taxa to plot. If |
   |                                | provided, it will override ``N``   |
   |                                | (default ``NULL``).                |
   +--------------------------------+------------------------------------+
   | ``others``                     | logical. Collapse the abundances   |
   |                                | of least abundant taxa, and        |
   |                                | include the result in the plot     |
   |                                | (default ``TRUE``).                |
   +--------------------------------+------------------------------------+
   | ``samples``                    | character. Character vector with   |
   |                                | the names of the samples to        |
   |                                | include in the plot. Can also be   |
   |                                | used to plot the samples in a      |
   |                                | custom order. If not provided, all |
   |                                | samples will be plotted (default   |
   |                                | ``NULL``).                         |
   +--------------------------------+------------------------------------+
   | ``nocds``                      | character. Either                  |
   |                                | ``"treat_separately"`` to treat    |
   |                                | reads annotated as No CDS          |
   |                                | separately,                        |
   |                                | ``"treat_as_unclassified"`` to     |
   |                                | treat them as Unclassified or      |
   |                                | ``"ignore"`` to ignore them in the |
   |                                | plot (default                      |
   |                                | ``"treat_separately"``).           |
   +--------------------------------+------------------------------------+
   | ``ignore_unmapped``            | logical. Don't include unmapped    |
   |                                | reads in the plot (default         |
   |                                | ``FALSE``).                        |
   +--------------------------------+------------------------------------+
   | ``ignore_unclassified``        | logical. Don't include             |
   |                                | unclassified reads in the plot     |
   |                                | (default ``FALSE``).               |
   +--------------------------------+------------------------------------+
   | ``no_partial_classifications`` | logical. Treat reads not fully     |
   |                                | classified at the requested level  |
   |                                | (e.g. "Unclassified Bacteroidota"  |
   |                                | at the class level or below) as    |
   |                                | fully unclassified. This takes     |
   |                                | effect before                      |
   |                                | ``ignore_unclassified``, so if     |
   |                                | both are ``TRUE`` the plot will    |
   |                                | only contain fully classified      |
   |                                | contigs (default ``FALSE``).       |
   +--------------------------------+------------------------------------+
   | ``rescale``                    | logical. Re-scale results to       |
   |                                | percentages (default ``FALSE``).   |
   +--------------------------------+------------------------------------+
   | ``color``                      | Vector with custom colors for the  |
   |                                | different features. If empty, we   |
   |                                | will use our own hand-picked       |
   |                                | pallete if N<=15, and the default  |
   |                                | ggplot2 palette otherwise (default |
   |                                | ``NULL``).                         |
   +--------------------------------+------------------------------------+
   | ``base_size``                  | numeric. Base font size (default   |
   |                                | ``11``).                           |
   +--------------------------------+------------------------------------+
   | ``max_scale_value``            | numeric. Maximum value to include  |
   |                                | in the y axis. By default it is    |
   |                                | handled automatically by ggplot2   |
   |                                | (default ``NULL``).                |
   +--------------------------------+------------------------------------+
   | ``metadata_groups``            | list. Split the plot into groups   |
   |                                | defined by the user: list('G1' =   |
   |                                | c('sample1', sample2'), 'G2' =     |
   |                                | c('sample3', 'sample4')) default   |
   |                                | ``NULL``).                         |
   +--------------------------------+------------------------------------+

   .. rubric:: Value
      :name: value

   a ggplot2 plot object.

   .. rubric:: See Also
      :name: see-also

   ``plotFunctions`` for plotting the most abundant functions of a SQM
   object; ``plotBars`` and ``plotHeatmap`` for plotting barplots or
   heatmaps with arbitrary data.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      Hadza.amin = subsetFun(Hadza, "Amino acid metabolism")
      # Taxonomic distribution of amino acid metabolism ORFs at the family level.
      plotTaxonomy(Hadza.amin, "family")
