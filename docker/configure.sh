#!/usr/bin/env bash 
 
CFLAGS="-g -Wall -Werror-implicit-function-declaration -fno-common -Wextra -Wunused" \ 
            CXXFLAGS="-g -Wall"  \ 
            ./configure --prefix=/usr/local \ 
            --with-postgres --with-postgres-includes=/usr/include/postgresql \ 
            --with-gdal \ 
            --with-proj \ 
            --with-cxx --enable-largefile --with-sqlite \ 
            --with-motif --with-glw --with-nls --with-readline \ 
            --with-freetype --with-freetype-includes=/usr/include/freetype2 \ 
            --with-odbc --with-python=/usr/bin/python-config --with-wxwidgets \ 
            --with-geos --with-pthread \ 
            --with-cairo \ 
            --with-readline-libs=/lib/x86_64-linux-gnu --with-openmp \ 
            --with-netcdf --with-x --x-includes=/usr/include --x-libraries=/usr/lib/x86_64-linux-gnu \ 
            --with-lapack --with-blas --with-bzlib --without-pdal 
 
exit 0
