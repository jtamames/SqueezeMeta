SQM_DIR=$PREFIX/SqueezeMeta
PYTHON=$PREFIX/bin/python
mkdir $SQM_DIR
cp -r $SRC_DIR/* $SQM_DIR
cd $PREFIX/bin
ln -s $SQM_DIR/scripts/* .
ln -s $SQM_DIR/utils/*.pl .
ln -s $SQM_DIR/utils/*.py .
ln -s $SQM_DIR/utils/install_utils/* .
ln -s $SQM_DIR/utils/binning_utils/* .
ln -s $SQM_DIR/utils/anvio_utils/* .
ln -s $SQM_DIR/bin/pogenom.pl .
ln -s $SQM_DIR/bin/run-pogenom.py .
ln -s $SQM_DIR/bin/checkm2 .
R CMD INSTALL $SQM_DIR/bin/DAS_Tool/package/DASTool_*.tar.gz
R CMD INSTALL $SQM_DIR/lib/SQMtools
cd $SQM_DIR/bin/CONCOCT-1.1.0/
$PYTHON setup.py install
cpanm Linux::MemInfo -L $PREFIX
$PYTHON -m pip install https://files.pythonhosted.org/packages/39/d0/0bb235071b241f7294bfacc54623a98afddf855ba36807a8dae1717b5694/mOTUlizer-0.2.4.tar.gz # we let mOTUlizer pull also igraph
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/5c/3a/8729cac4fac0b3af6e91c5c5a542bfb0268cbdf53d0fe6fde7417b49de6f/speedict-0.3.12-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/22/46/16da52fd09105dc2f9e29028e4971eeeab6b4a1ebcee17c837b7b09c24b2/superpang-1.3.1.tar.gz
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/c6/bd/b66510b12e1b143fcd6fd45a78a8c50cf8750e95d98ebcb5b496421ab279/gtdbtk-2.4.0.tar.gz
## Check that we won't break anything when manually moving libraries
#if [ -d $PREFIX/lib/perl5/5.32/site_perl/ ]; then exit 1; fi # This would fail on purpose
if [ -d $PREFIX/lib/perl5/5.32/site_perl/Linux ]; then exit 1; fi
if [ -d $PREFIX/lib/perl5/5.32/site_perl/auto ]; then exit 1; fi
if [ -d $PREFIX/lib/perl5/5.32/site_perl/x86_64-linux-thread-multi/ ]; then exit 1; fi
mv $PREFIX/lib/perl5/Linux $PREFIX/lib/perl5/5.32/site_perl/
mv $PREFIX/lib/perl5/auto $PREFIX/lib/perl5/5.32/site_perl/
mv $PREFIX/lib/perl5/x86_64-linux-thread-multi/ $PREFIX/lib/perl5/5.32/site_perl/
#cpanm Term::ANSIColor -L $PREFIX # Apparently this is already installed
## Override PERL5LIB to priorize conda libraries
echo "export PERL5LIB=\$CONDA_PREFIX/lib/perl5/5.32/site_perl/:\$CONDA_PREFIX/lib/perl5/5.32/site_perl/x86_64-linux-thread-multi/:\$CONDA_PREFIX/lib/perl5/5.32/:\$PERL5LIB" > $PREFIX/etc/conda/activate.d/activate-perl.sh
echo "export PERL5LIB=\`echo \$PERL5LIB | sed -e \"s,\$CONDA_PREFIX.*\:,,\"\`" > $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
chmod +x $PREFIX/etc/conda/activate.d/activate-perl.sh
chmod +x $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
# fix libs for samtools and mothur
cd $PREFIX/lib
ln -s libncurses.so.6 libncurses.so.5 
ln -s libtinfo.so.6 libtinfo.so.5
ln -s libreadline.so.8 libreadline.so.6
cp libboost_iostreams.so libboost_iostreams.so.1.85.0
cp libboost_system.so libboost_system.so.1.85.0
cp libboost_filesystem.so libboost_filesystem.so.1.85.0
