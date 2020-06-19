/// @file
///
/// Package config file. ONLY include in ".cc" files never in header files!

// std include
#include <string>

#ifndef ASKAP_SCIMATH_H
#define ASKAP_SCIMATH_H

  /// The name of the package
#define ASKAP_PACKAGE_NAME "scimath"

/// askap namespace
namespace askap {
  /// @return version of the package
  std::string getAskapPackageVersion_scimath();
}

  /// The version of the package
#define ASKAP_PACKAGE_VERSION askap::getAskapPackageVersion_scimath()

#endif
