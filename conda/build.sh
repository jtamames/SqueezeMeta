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
