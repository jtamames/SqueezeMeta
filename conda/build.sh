SQM_DIR=$PREFIX/SqueezeMeta
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
R CMD INSTALL $SQM_DIR/bin/DAS_Tool/package/DASTool_*.tar.gz
R CMD INSTALL $SQM_DIR/lib/SQMtools
cd $SQM_DIR/bin/CONCOCT-1.1.0/
$PYTHON setup.py install
cpanm Linux::MemInfo -L $PREFIX
#$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/6f/9c/5f7d48cbf6c0400abb1a3414cfee8e557600391840f06bbb497c49082d9d/igraph-0.10.8-pp310-pypy310_pp73-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
$PYTHON -m pip install https://files.pythonhosted.org/packages/39/d0/0bb235071b241f7294bfacc54623a98afddf855ba36807a8dae1717b5694/mOTUlizer-0.2.4.tar.gz # we let mOTUlizer pull also igraph
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/f1/2a/3a5c44c25388238cb001e6968de008d765fdd9af881190388421097c84d9/speedict-0.3.11-cp310-cp310-manylinux_2_28_x86_64.whl --no-deps
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/40/b9/b8fa8aa5244430270e0da993b35db855f9e383ef9e174123cd5daf90ea98/SuperPang-1.1.0.post3.tar.gz
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/5a/4b/2ed5ee6a35a0f9b4ec896e0b796aeb931824d82b9ca236665e5839b9a7a7/gtdbtk-2.3.2.tar.gz
# Check that we won't break anything when manually moving libraries
#if [ -d $PREFIX/lib/site_perl/5.26.2/ ]; then exit 1; fi # This would fail on purpose
if [ -d $PREFIX/lib/site_perl/5.26.2/Linux ]; then exit 1; fi
if [ -d $PREFIX/lib/site_perl/5.26.2/auto ]; then exit 1; fi
if [ -d $PREFIX/lib/site_perl/5.26.2/x86_64-linux-thread-multi/auto/Linux ]; then exit 1; fi
mv $PREFIX/lib/perl5/Linux $PREFIX/lib/site_perl/5.26.2/
mv $PREFIX/lib/perl5/auto $PREFIX/lib/site_perl/5.26.2/
mv $PREFIX/lib/perl5/x86_64-linux-thread-multi/auto/Linux $PREFIX/lib/site_perl/5.26.2/x86_64-linux-thread-multi/auto/
rm -r $PREFIX/lib/perl5/
#cpanm Term::ANSIColor -L $PREFIX # Apparently this is already installed
# Override PERL5LIB to priorize conda libraries
echo "export PERL5LIB=\$CONDA_PREFIX/lib/site_perl/5.26.2/:\$CONDA_PREFIX/lib/site_perl/5.26.2/x86_64-linux-thread-multi/:\$CONDA_PREFIX/lib/5.26.2/:\$CONDA_PREFIX/lib/5.26.2/x86_64-linux-thread-multi/:\$PERL5LIB" > $PREFIX/etc/conda/activate.d/activate-perl.sh
echo "export PERL5LIB=\`echo \$PERL5LIB | sed -e \"s,\$CONDA_PREFIX.*\:,,\"\`" > $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
chmod +x $PREFIX/etc/conda/activate.d/activate-perl.sh
chmod +x $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
# fix libs for samtools and mothur
cd $PREFIX/lib
ln -s libncurses.so.6 libncurses.so.5 
ln -s libtinfo.so.6 libtinfo.so.5
ln -s libreadline.so.8 libreadline.so.6
