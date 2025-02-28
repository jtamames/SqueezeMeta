***********
plotHeatmap
***********

.. container::

   =========== ===============
   plotHeatmap R Documentation
   =========== ===============

   .. rubric:: Plot a heatmap using ggplot2
      :name: plotHeatmap

   .. rubric:: Description
      :name: description

   Plot a ggplot2 heatmap from a matrix or data frame. The data should
   be in tabular format (e.g. features in rows and samples in columns).

   .. rubric:: Usage
      :name: usage

   .. code:: R

      plotHeatmap(
        data,
        label_x = "Samples",
        label_y = "Features",
        label_fill = "Abundance",
        gradient_col = c("ghostwhite", "dodgerblue4"),
        base_size = 11,
        metadata_groups = NULL
      )

   .. rubric:: Arguments
      :name: arguments

   +---------------------+-----------------------------------------------+
   | ``data``            | numeric matrix or data frame.                 |
   +---------------------+-----------------------------------------------+
   | ``label_x``         | character Label for the x axis (default       |
   |                     | ``"Samples"``).                               |
   +---------------------+-----------------------------------------------+
   | ``label_y``         | character Label for the y axis (default       |
   |                     | ``"Features"``).                              |
   +---------------------+-----------------------------------------------+
   | ``label_fill``      | character Label for color scale (default      |
   |                     | ``"Abundance"``).                             |
   +---------------------+-----------------------------------------------+
   | ``gradient_col``    | A vector of two colors representing the low   |
   |                     | and high ends of the color gradient (default  |
   |                     | ``c("ghostwhite", "dodgerblue4")``).          |
   +---------------------+-----------------------------------------------+
   | ``base_size``       | numeric. Base font size (default ``11``).     |
   +---------------------+-----------------------------------------------+
   | ``metadata_groups`` | list. Split the plot into groups defined by   |
   |                     | the user: list('G1' = c('sample1', sample2'), |
   |                     | 'G2' = c('sample3', 'sample4')) default       |
   |                     | ``NULL``).                                    |
   +---------------------+-----------------------------------------------+

   .. rubric:: Value
      :name: value

   A ggplot2 plot object.

   .. rubric:: See Also
      :name: see-also

   ``plotFunctions`` for plotting the top functional categories of a SQM
   object; ``plotBars`` for plotting a barplot with arbitrary data;
   ``mostAbundant`` for selecting the most abundant rows in a dataframe
   or matrix.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      topPFAM = mostAbundant(Hadza$functions$PFAM$tpm)
      topPFAM = topPFAM[rownames(topPFAM) != "Unclassified",] # Take out the Unclassified ORFs.
      plotHeatmap(topPFAM, label_x = "Samples", label_y = "PFAMs", label_fill = "TPM")
