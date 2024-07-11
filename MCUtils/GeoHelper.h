#ifndef GEO_HELPER_H
#define GEO_HELPER_H

// Standard Library Includes
#include <iostream>
#include <fstream>

// // ROOT Includes
#include "TObject.h"
#include "TVector3.h"


namespace analysis {


class GeoHelper {

public:
        GeoHelper() {}
        ~GeoHelper() {}
	

	// return if a position is in the active volume
	bool isActive(double pos_x, double pos_y, double pos_z) {

	  // Get the Geometry
	  art::ServiceHandle<geo::Geometry> geom;

	  // Get active volume boundary.
	  double xmin = -2.0 * geom->DetHalfWidth() - 1e-8;
	  double xmax = 2.0 * geom->DetHalfWidth() + 1e-8;
	  double ymin = -geom->DetHalfHeight() -1e-8;
	  double ymax = geom->DetHalfHeight() + 1e-8;
	  double zmin = 0. -1e-8;
	  double zmax = geom->DetLength() + 1e-8;


	  if ( xmin <= pos_x && pos_x <= xmax
	    && ymin <= pos_y && pos_y <= ymax
	    && zmin <= pos_z && pos_z <= zmax ) {

	    return true;

	  }
	  else {
	    return false;
	  }

	} // end of isActive



private:


};

} // end namespace analysis

#endif
