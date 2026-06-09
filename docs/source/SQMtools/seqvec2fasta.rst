************
seqvec2fasta
************

============ ===============
seqvec2fasta R Documentation
============ ===============

Print a named vector of sequences as a fasta-formatted string
-------------------------------------------------------------

Description
~~~~~~~~~~~

Print a named vector of sequences as a fasta-formatted string

Usage
~~~~~

.. code:: R

   seqvec2fasta(seqvec, output_name = "")

Arguments
~~~~~~~~~

+-----------------+----------------------------------------------------+
| ``seqvec``      | vector. The vector to be written as a fasta        |
|                 | string.                                            |
+-----------------+----------------------------------------------------+
| ``output_name`` | A connection, or a character string naming the     |
|                 | file to print to. If "" (the default), sequences   |
|                 | will be printed to the standard output connection. |
+-----------------+----------------------------------------------------+

Examples
~~~~~~~~

.. code:: R

   data(Hadza)
   seqvec2fasta(Hadza$orfs$seqs[1:10])
