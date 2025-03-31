************
seqvec2fasta
************

.. container::

   ============ ===============
   seqvec2fasta R Documentation
   ============ ===============

   .. rubric:: Print a named vector of sequences as a fasta-formatted
      string
      :name: seqvec2fasta

   .. rubric:: Description
      :name: description

   Print a named vector of sequences as a fasta-formatted string

   .. rubric:: Usage
      :name: usage

   .. code:: R

      seqvec2fasta(seqvec, output_name = "")

   .. rubric:: Arguments
      :name: arguments

   +-----------------+---------------------------------------------------+
   | ``seqvec``      | vector. The vector to be written as a fasta       |
   |                 | string.                                           |
   +-----------------+---------------------------------------------------+
   | ``output_name`` | A connection, or a character string naming the    |
   |                 | file to print to. If "" (the default), sequences  |
   |                 | will be printed to the standard output            |
   |                 | connection.                                       |
   +-----------------+---------------------------------------------------+

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      seqvec2fasta(Hadza$orfs$seqs[1:10])
