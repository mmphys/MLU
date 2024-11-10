# MLU

# Include MLU's header when compiling this directory.
#AM_CPPFLAGS = -I$(ac_top_srcdir)

# Build the library in the MLU subdirectory first.
SUBDIRS =
if WANT_BUILD_GSL
SUBDIRS += gsl-2.7
endif
SUBDIRS += MLU Analyse Fit

nobase_include_HEADERS = MLU/CorrFile.hpp MLU/DataSet.hpp \
  MLU/DebugInfo.hpp MLU/FitRange.hpp MLU/FitRangeImp.hpp \
  MLU/Fold.hpp MLU/GSLVecMat.hpp MLU/GSLVecMatImp.hpp MLU/HDF5.hpp \
  MLU/JackBoot.hpp MLU/Math.hpp MLU/MLU.hpp MLU/MLUFirst.hpp \
  MLU/Model.hpp MLU/Param.hpp MLU/Physics.hpp MLU/Posix.hpp \
  MLU/Sample.hpp MLU/String.hpp MLU/Utility.hpp
