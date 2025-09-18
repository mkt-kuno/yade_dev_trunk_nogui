/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <lib/base/Logging.hpp>
#include <lib/multimethods/Indexable.hpp>
#include <lib/serialization/Serializable.hpp>
#include <core/Dispatcher.hpp>

#define BV_FUNCTOR_CACHE

namespace yade { // Cannot have #include directive inside.

class BoundFunctor;

class InternalForceFunctor;

class Shape : public Serializable, public Indexable {
public:
	virtual ~Shape() {}; // vtable
#ifdef BV_FUNCTOR_CACHE
	shared_ptr<BoundFunctor> boundFunctor;
#endif
#ifdef YADE_FEM
	//! cache functor that are called for this type of DeformableElement. Used by FEInternalForceEngine
	shared_ptr<InternalForceFunctor> internalforcefunctor;
#endif
	virtual Real getVolume()
	{ // a number of overriden methods in child classes may execute some update operations before returning a value, and constness can not apply
		LOG_ERROR( // LOG_ERROR shall suffice here ? Maybe the calculation can still continue (with the results being all wrong as with LOG_ERROR)
		        "Shape " << getClassName()
		                 << " calling virtual method Shape::getVolume(). This shall not happen, please submit bug report at "
		                    "https://gitlab.com/yade-dev/trunk/issues");
		return -1.;
	}

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Shape,Serializable,"Geometry of a body",
		((Vector3r,color,Vector3r(1,1,1),,"Color for rendering (normalized RGB)."))
		((bool,wire,false,,"Whether this Shape is rendered using color surfaces, or only wireframe (can still be overridden by global config of the renderer)."))
		((bool,highlight,false,,"Whether this Shape will be highlighted when rendered.")),
		/*ctor*/,
		/*py*/ YADE_PY_TOPINDEXABLE(Shape)
	);
	// clang-format on
	REGISTER_INDEX_COUNTER(Shape);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Shape);

} // namespace yade
