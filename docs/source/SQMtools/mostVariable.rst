************
mostVariable
************

============ ===============
mostVariable R Documentation
============ ===============

Get the N most variable rows (or columns) from a numeric table
--------------------------------------------------------------

Description
~~~~~~~~~~~

Return a subset of an input matrix or data frame, containing only the N
most variable rows (or columns), sorted. Variability is calculated as
the Coefficient of Variation (sd/mean).

Usage
~~~~~

.. code:: R

   mostVariable(data, N = 10, bycol = FALSE)

Arguments
~~~~~~~~~

+-----------+----------------------------------------------------------+
| ``data``  | numeric matrix or data frame                             |
+-----------+----------------------------------------------------------+
| ``N``     | integer Number of rows to return (default ``10``).       |
+-----------+----------------------------------------------------------+
| ``bycol`` | logical. Operate on columns instead of rows (default     |
|           | ``FALSE``).                                              |
+-----------+----------------------------------------------------------+

Value
~~~~~

A matrix or data frame (same as input) with the selected rows or
columns.

Examples
~~~~~~~~

.. code:: R

   data(Hadza)
   Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
   # Which are the 20 most variable KEGG functions in the ORFs related to carbohydrate metabolism?
   topCarb = mostVariable(Hadza.carb$functions$KEGG$tpm, N=20)
   # Now print them with nice names
   rownames(topCarb) = paste(rownames(topCarb),
                             Hadza.carb$misc$KEGG_names[rownames(topCarb)], sep="; ")
   topCarb
   # We can pass this to any R function
   heatmap(topCarb)
   # But for convenience we provide wrappers for plotting ggplot2 heatmaps and barplots
   plotHeatmap(topCarb, label_y="TPM")
   plotBars(topCarb, label_y="TPM")
