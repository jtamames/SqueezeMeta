# Changelog

A list of changes per version.
When changing something in the develop branch, it should be added here.

## [v1.1.0] 2019-08-02

### `Changed`

 - [#236](https://github.com/BinPro/CONCOCT/pull/236) - Always add suffix to contigs at cutup, even when they are not cut.
 - [#254](https://github.com/BinPro/CONCOCT/pull/254) - Slight cleanup of concoct refine
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - New suffices (.concoct_part_XX) are now used for contig parts 
 - [#261](https://github.com/BinPro/CONCOCT/pull/261) - Epsilon argument removed as it was not working and is not very useful
 - [#262](https://github.com/BinPro/CONCOCT/pull/262) - Rewrote documentation, including installation instructions
 - [#264](https://github.com/BinPro/CONCOCT/pull/264) - `concoct_part_` suffix is enforced in subcontig for coverage script 
 - [#264](https://github.com/BinPro/CONCOCT/pull/264) - Header line is enforced for input for `merge_cutup_clustering.py` and `extract_fasta_bins.py`
 - [#267](https://github.com/BinPro/CONCOCT/pull/267) - Updated documentation

### `Added`

 - [#253](https://github.com/BinPro/CONCOCT/pull/253) - A dockerfile useful to test the conda installation
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - Tests for all fundamental scripts, including a new integration test data repository
 - [#259](https://github.com/BinPro/CONCOCT/pull/259) - This changelog
 - [#262](https://github.com/BinPro/CONCOCT/pull/262) - Added documentation for the core scripts used with concoct
 - [#265](https://github.com/BinPro/CONCOCT/pull/265) - A warning is now printed when concoct runs in single threaded mode

### `Fixed`

 - [#230](https://github.com/BinPro/CONCOCT/pull/230) - Enable at least single threaded installation on Mac OSX
 - [#231](https://github.com/BinPro/CONCOCT/pull/231) - Replace pandas .ix with .loc to fix deprecation warnings
 - [#246](https://github.com/BinPro/CONCOCT/pull/246) - Limit some dependency version numbers for python 2
 - [#254](https://github.com/BinPro/CONCOCT/pull/254) - Concoct refine now works with python 3
 - [#258](https://github.com/BinPro/CONCOCT/pull/258) - Seed tests now working again
 - [#260](https://github.com/BinPro/CONCOCT/pull/260) - Fix the dockerfile build by adding integration test data
