************
mostAbundant
************

.. container::

   ============ ===============
   mostAbundant R Documentation
   ============ ===============

   .. rubric:: Get the N most abundant rows (or columns) from a numeric
      table
      :name: mostAbundant

   .. rubric:: Description
      :name: description

   Return a subset of an input matrix or data frame, containing only the
   N most abundant rows (or columns), sorted. Alternatively, a custom
   set of rows can be returned.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      mostAbundant(
        data,
        N = 10,
        items = NULL,
        extra_items = NULL,
        ignore = NULL,
        others = FALSE,
        rescale = FALSE,
        bycol = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +-----------------+---------------------------------------------------+
   | ``data``        | numeric matrix or data frame                      |
   +-----------------+---------------------------------------------------+
   | ``N``           | integer Number of rows to return (default         |
   |                 | ``10``).                                          |
   +-----------------+---------------------------------------------------+
   | ``items``       | character vector. Custom row names to return. If  |
   |                 | provided, it will override ``N`` and              |
   |                 | ``extra_items`` (default ``NULL``).               |
   +-----------------+---------------------------------------------------+
   | ``extra_items`` | character vector. Extra row names to return on    |
   |                 | top of the N most abundant (default ``NULL``)     |
   +-----------------+---------------------------------------------------+
   | ``ignore``      | character. Custom row names to drop before        |
   |                 | abundance calculation.                            |
   +-----------------+---------------------------------------------------+
   | ``others``      | logical. If ``TRUE``, an extra row will be        |
   |                 | returned containing the aggregated abundances of  |
   |                 | the elements not selected with ``N`` or ``items`` |
   |                 | (default ``FALSE``).                              |
   +-----------------+---------------------------------------------------+
   | ``rescale``     | logical. Scale result to percentages column-wise  |
   |                 | (default ``FALSE``).                              |
   +-----------------+---------------------------------------------------+
   | ``bycol``       | logical. Operate on columns instead of rows       |
   |                 | (default ``FALSE``).                              |
   +-----------------+---------------------------------------------------+

   .. rubric:: Value
      :name: value

   A matrix or data frame (same as input) with the selected rows (or
   columns).

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
      # Which are the 20 most abundant KEGG functions in the ORFs related to carbohydrate metabolism?
      topCarb = mostAbundant(Hadza.carb$functions$KEGG$tpm, N=20)
      # Now print them with nice names.
      rownames(topCarb) = paste(rownames(topCarb),
                                Hadza.carb$misc$KEGG_names[rownames(topCarb)], sep="; ")
      topCarb
      # We can pass this to any R function.
      heatmap(topCarb)
      # But for convenience we provide wrappers for plotting ggplot2 heatmaps and barplots.
      plotHeatmap(topCarb, label_y="TPM")
      plotBars(topCarb, label_y="TPM")
