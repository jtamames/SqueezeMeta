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
ln -s $SQM_DIR/bin/checkm2 .
R CMD INSTALL $SQM_DIR/bin/DAS_Tool/package/DASTool_*.tar.gz
R CMD INSTALL $SQM_DIR/lib/SQMtools
cd $SQM_DIR/bin/CONCOCT-1.1.0/
$PYTHON setup.py install
cpanm Linux::MemInfo -L $PREFIX
$PYTHON -m pip install https://files.pythonhosted.org/packages/39/d0/0bb235071b241f7294bfacc54623a98afddf855ba36807a8dae1717b5694/mOTUlizer-0.2.4.tar.gz # we let mOTUlizer pull also igraph
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/88/bf/cc506602798a6ee578e9ea54b06f12ac01f084be39b0615569447e49673f/speedict-0.3.12-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/4d/c9/307ff96cc6740664c6efad693ff780f5d8f14f8a5cfcf19039f54b5e7c43/SuperPang-1.3.0.tar.gz
$PYTHON -m pip install --no-deps https://files.pythonhosted.org/packages/c6/bd/b66510b12e1b143fcd6fd45a78a8c50cf8750e95d98ebcb5b496421ab279/gtdbtk-2.4.0.tar.gz
# Check that we won't break anything when manually moving libraries
#if [ -d $PREFIX/lib/site_perl/5.32.1/ ]; then exit 1; fi # This would fail on purpose
if [ -d $PREFIX/lib/site_perl/5.32.1/Linux ]; then exit 1; fi
if [ -d $PREFIX/lib/site_perl/5.32.1/auto ]; then exit 1; fi
if [ -d $PREFIX/lib/site_perl/5.32.1/x86_64-linux-thread-multi/auto/Linux ]; then exit 1; fi
mv $PREFIX/lib/perl5/Linux $PREFIX/lib/site_perl/5.32.1/
mv $PREFIX/lib/perl5/auto $PREFIX/lib/site_perl/5.32.1/
mv $PREFIX/lib/perl5/x86_64-linux-thread-multi/auto/Linux $PREFIX/lib/site_perl/5.32.1/x86_64-linux-thread-multi/auto/
rm -r $PREFIX/lib/perl5/
#cpanm Term::ANSIColor -L $PREFIX # Apparently this is already installed
# Override PERL5LIB to priorize conda libraries
echo "export PERL5LIB=\$CONDA_PREFIX/lib/site_perl/5.32.1/:\$CONDA_PREFIX/lib/site_perl/5.32.1/x86_64-linux-thread-multi/:\$CONDA_PREFIX/lib/5.32.1/:\$CONDA_PREFIX/lib/5.32.1/x86_64-linux-thread-multi/:\$PERL5LIB" > $PREFIX/etc/conda/activate.d/activate-perl.sh
echo "export PERL5LIB=\`echo \$PERL5LIB | sed -e \"s,\$CONDA_PREFIX.*\:,,\"\`" > $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
chmod +x $PREFIX/etc/conda/activate.d/activate-perl.sh
chmod +x $PREFIX/etc/conda/deactivate.d/deactivate-perl.sh
# fix libs for samtools and mothur
cd $PREFIX/lib
ln -s libncurses.so.6 libncurses.so.5 
ln -s libtinfo.so.6 libtinfo.so.5
ln -s libreadline.so.8 libreadline.so.6
