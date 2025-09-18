#pragma once
#include <core/Shape.hpp>
#include <lib/high-precision/Constants.hpp>

// FIXED: add namespace yade, see https://gitlab.com/yade-dev/trunk/issues/57 and (old site, fixed bug) https://bugs.launchpad.net/yade/+bug/528509

namespace yade {

class Sphere : public Shape {
public:
	Sphere(Real _radius)
	        : radius(_radius)
	{
		createIndex();
	}
	virtual ~Sphere() {};
	Real getVolume() override {return (4 / 3.) * Mathr::PI * pow(radius, 3);}
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Sphere,Shape,"Geometry of spherical particle.",
		((Real,radius,NaN,,"Radius [m]")),
		createIndex(); /*ctor*/
		,
		.def("getVolume",&Sphere::getVolume,"Returns the shape volume.")
	);
	// clang-format on
	REGISTER_CLASS_INDEX(Sphere, Shape);
};

REGISTER_SERIALIZABLE(Sphere);

}; // namespace yade
