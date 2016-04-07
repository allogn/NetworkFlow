# Installation

Using boost > 1.60.0 (for rtree with iterative queries)

cmake -DCMAKE_INSTALL_PREFIX=/usr/local/boost_1_60_0/bin.v2 -DBoost_NO_BOOST_CMAKE=TRUE -DBoost_NO_SYSTEM_PATHS=TRUE -DBOOST_ROOT:PATHNAME=/usr/local/boost_1_60_0 -DBoost_LIBRARY_DIRS:FILEPATH=/usr/local/boost_1_60_0/bin.v2/libs
