***********************
remove_contigs_from_bin
***********************

.. container::

   ======================= ===============
   remove_contigs_from_bin R Documentation
   ======================= ===============

   .. rubric:: Remove contigs from a given bin
      :name: remove_contigs_from_bin

   .. rubric:: Description
      :name: description

   Remove contigs from a given bin

   .. rubric:: Usage
      :name: usage

   .. code:: R

      remove_contigs_from_bin(SQM, bin, contigs)

   .. rubric:: Arguments
      :name: arguments

   +-------------+-------------------------------------------------------+
   | ``SQM``     | A SQM object.                                         |
   +-------------+-------------------------------------------------------+
   | ``bin``     | character. Name of the bin from which the contigs     |
   |             | will be removed.                                      |
   +-------------+-------------------------------------------------------+
   | ``contigs`` | character. Vector with the names of the contigs that  |
   |             | will be removed from the new bin.                     |
   +-------------+-------------------------------------------------------+

   .. rubric:: Value
      :name: value

   SQM object with the new binning information, including recalculated
   bin statistics if possible.

   .. rubric:: See Also
      :name: see-also

   ``find_redundant_contigs``, ``create_bin``
