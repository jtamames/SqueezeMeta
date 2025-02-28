**********************
find_redundant_contigs
**********************

.. container::

   ====================== ===============
   find_redundant_contigs R Documentation
   ====================== ===============

   .. rubric:: Find redundant contigs within a bin
      :name: find_redundant_contigs

   .. rubric:: Description
      :name: description

   Find contigs with overlapping marker genes in a given bin, and
   suggest contigs to be removed in order to reduce contamination
   without affecting completeness. Note that this can give a quick idea
   of the contigs that are sources of contamination within a bin, but is
   not a replacement for proper bin refininement with other tools such
   as anvi\\'o.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      find_redundant_contigs(SQM, bin, minimum_overlap_for_removal = 1)

   .. rubric:: Arguments
      :name: arguments

   +---------------------------------+-----------------------------------+
   | ``SQM``                         | A SQM object.                     |
   +---------------------------------+-----------------------------------+
   | ``bin``                         | character. Name of the bin to be  |
   |                                 | created.                          |
   +---------------------------------+-----------------------------------+
   | ``minimum_overlap_for_removal`` | numeric. Fraction of marker genes |
   |                                 | in the contigs present in another |
   |                                 | contig needed to suggest it for   |
   |                                 | removal. If set to ``1``          |
   |                                 | (default), contigs will only      |
   |                                 | suggested for removal if their    |
   |                                 | markers fully overlap with those  |
   |                                 | in another contig (and thus       |
   |                                 | completeness will not change      |
   |                                 | after removing them). Smaller     |
   |                                 | values will result in more        |
   |                                 | contigs being suggested for       |
   |                                 | removal, which will further       |
   |                                 | reduce contamination at the       |
   |                                 | expense of some completeness.     |
   +---------------------------------+-----------------------------------+

   .. rubric:: Value
      :name: value

   A character vector with the contigs deemed to be redundant. A heatmap
   showing how marker genes overlap over different contigs will also be
   produced.

   .. rubric:: See Also
      :name: see-also

   ``create_bin``, ``remove_contigs_from_bin``

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      bin_name = "Hadza2merged.concoct.28.fa.contigs"
      # Get redundant contigs that could be removed from our bin
      candidates_for_removal = find_redundant_contigs(Hadza, bin_name)
      # We can now remove them from the bin
      Hadza.new.1 = remove_contigs_from_bin(Hadza, bin_name, candidates_for_removal)
      # Or we can create a new bin out of them
      #  which will also remove them from the original bin
      Hadza.new.2 = create_bin(Hadza, "new_bin_name", candidates_for_removal)
