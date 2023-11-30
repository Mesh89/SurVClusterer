tar -xjf htslib-1.18.tar.bz2
cd htslib-1.18
wget -O config.guess 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD'
wget -O config.sub 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD'
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
