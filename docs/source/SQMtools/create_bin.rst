**********
create_bin
**********

.. container::

   ========== ===============
   create_bin R Documentation
   ========== ===============

   .. rubric:: Create a bin from a vector of contigs
      :name: create_bin

   .. rubric:: Description
      :name: description

   Create a bin from a vector of contigs

   .. rubric:: Usage
      :name: usage

   .. code:: R

      create_bin(SQM, bin, contigs, delete_overlapping_bins = FALSE)

   .. rubric:: Arguments
      :name: arguments

   +-----------------------------+---------------------------------------+
   | ``SQM``                     | A SQM object.                         |
   +-----------------------------+---------------------------------------+
   | ``bin``                     | character. Name of the bin to be      |
   |                             | created.                              |
   +-----------------------------+---------------------------------------+
   | ``contigs``                 | character. Vector with the names of   |
   |                             | the contigs that will be included in  |
   |                             | the new bin.                          |
   +-----------------------------+---------------------------------------+
   | ``delete_overlapping_bins`` | logical. If ``TRUE``, bins that       |
   |                             | originally contained any of the       |
   |                             | provided contigs will be removed from |
   |                             | the object. Default ``FALSE``.        |
   +-----------------------------+---------------------------------------+

   .. rubric:: Value
      :name: value

   SQM object with the new binning information, including recalculated
   bin statistics if possible.

   .. rubric:: See Also
      :name: see-also

   ``find_redundant_contigs``, ``remove_contigs_from_bin``
