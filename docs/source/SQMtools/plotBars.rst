********
plotBars
********

.. container::

   ======== ===============
   plotBars R Documentation
   ======== ===============

   .. rubric:: Plot a barplot using ggplot2
      :name: plotBars

   .. rubric:: Description
      :name: description

   Plot a ggplot2 barplot from a matrix or data frame. The data should
   be in tabular format (e.g. features in rows and samples in columns).

   .. rubric:: Usage
      :name: usage

   .. code:: R

      plotBars(
        data,
        label_x = "Samples",
        label_y = "Abundances",
        label_fill = "Features",
        color = NULL,
        base_size = 11,
        max_scale_value = NULL,
        metadata_groups = NULL
      )

   .. rubric:: Arguments
      :name: arguments

   +---------------------+-----------------------------------------------+
   | ``data``            | Numeric matrix or data frame.                 |
   +---------------------+-----------------------------------------------+
   | ``label_x``         | character Label for the x axis (default       |
   |                     | ``"Samples"``).                               |
   +---------------------+-----------------------------------------------+
   | ``label_y``         | character Label for the y axis (default       |
   |                     | ``"Abundances"``).                            |
   +---------------------+-----------------------------------------------+
   | ``label_fill``      | character Label for color categories (default |
   |                     | ``"Features"``).                              |
   +---------------------+-----------------------------------------------+
   | ``color``           | Vector with custom colors for the different   |
   |                     | features. If empty, the default ggplot2       |
   |                     | palette will be used (default ``NULL``).      |
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

   ``plotTaxonomy`` for plotting the most abundant taxa of a SQM object;
   ``plotHeatmap`` for plotting a heatmap with arbitrary data;
   ``mostAbundant`` for selecting the most abundant rows in a dataframe
   or matrix.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      sk = Hadza$taxa$superkingdom$abund
      plotBars(sk, label_y = "Raw reads", label_fill = "Superkingdom")
