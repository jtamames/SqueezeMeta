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
R CMD INSTALL $SQM_DIR/bin/DAS_Tool/package/DASTool_*.tar.gz
R CMD INSTALL $SQM_DIR/lib/SQMtools
cd $SQM_DIR/bin/CONCOCT-1.1.0/
python3 setup.py install
cpanm Linux::MemInfo -L $PREFIX
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
