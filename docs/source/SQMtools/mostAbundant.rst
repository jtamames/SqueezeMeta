************
mostAbundant
************

============ ===============
mostAbundant R Documentation
============ ===============

Get the N most abundant rows (or columns) from a numeric table
--------------------------------------------------------------

Description
~~~~~~~~~~~

Return a subset of an input matrix or data frame, containing only the N
most abundant rows (or columns), sorted. Alternatively, a custom set of
rows can be returned.

Usage
~~~~~

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

Arguments
~~~~~~~~~

+-----------------+----------------------------------------------------+
| ``data``        | numeric matrix or data frame                       |
+-----------------+----------------------------------------------------+
| ``N``           | integer Number of rows to return (default ``10``). |
+-----------------+----------------------------------------------------+
| ``items``       | character vector. Custom row names to return. If   |
|                 | provided, it will override ``N`` and               |
|                 | ``extra_items`` (default ``NULL``).                |
+-----------------+----------------------------------------------------+
| ``extra_items`` | character vector. Extra row names to return on top |
|                 | of the N most abundant (default ``NULL``)          |
+-----------------+----------------------------------------------------+
| ``ignore``      | character. Custom row names to drop before         |
|                 | abundance calculation.                             |
+-----------------+----------------------------------------------------+
| ``others``      | logical. If ``TRUE``, an extra row will be         |
|                 | returned containing the aggregated abundances of   |
|                 | the elements not selected with ``N`` or ``items``  |
|                 | (default ``FALSE``).                               |
+-----------------+----------------------------------------------------+
| ``rescale``     | logical. Scale result to percentages column-wise   |
|                 | (default ``FALSE``).                               |
+-----------------+----------------------------------------------------+
| ``bycol``       | logical. Operate on columns instead of rows        |
|                 | (default ``FALSE``).                               |
+-----------------+----------------------------------------------------+

Value
~~~~~

A matrix or data frame (same as input) with the selected rows (or
columns).

Examples
~~~~~~~~

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
