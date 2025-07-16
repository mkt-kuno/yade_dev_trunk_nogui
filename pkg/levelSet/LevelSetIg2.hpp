/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/levelSet/LevelSet.hpp>
#include <pkg/levelSet/LevelSetIGeom.hpp>

namespace yade {

class Ig2_Sphere_LevelSet_ScGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override;
	bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c) override
	{
		c->swapOrder();
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	};
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Sphere_LevelSet_ScGeom,IGeomFunctor,R"""(
	Creates or updates a :yref:`ScGeom` instance representing the intersection of one :yref:`LevelSet`-shaped body with one :yref:`Sphere`-shaped body, where overlap is always chosen to occur inside the level set body (i.e. spheres will always be expelled). 
	:yref:`Contact normal<ScGeom.normal>` $\vec{n}$ is given by the level set normal at the centre of the sphere $\vec{c}$ while :yref:`overlap<ScGeom.penetrationDepth>` is given by $R - \varphi$ with $R$ the radius and $\varphi$ the level set value. 
	And :yref:`contact points<ScGeom.contactPoint>` is defined as $\vec{c}-\varphi \vec{n}$. This functionality does not require the level set to have surface nodes. 
	Approximations for $\phi$ outside the level set may become inaccurate if the spheres are of similar size or larger than the level set body. 
	Accuracy is guaranteed it the largest sphere is around the same size, or smaller than, the smallest grid cell in the level set.)""");
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Sphere, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(Sphere, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_Sphere_LevelSet_ScGeom);

class Ig2_LevelSet_LevelSet_ScGeom : public IGeomFunctor {
public:
	static void geomPtrForLaterRemoval(
	        const State&                   rbp1,
	        const State&                   rbp2,
	        const shared_ptr<Interaction>& c,
	        const Scene*                   scene); // static because also used in e.g., Ig2_Box_LevelSet_ScGeom
	std::pair<std::pair<Vector3r, Vector3r>, std::pair<bool, bool>>
	                                boundOverlap(bool single, const State&, const State&, const shared_ptr<Interaction>&, const Vector3r&);
	std::tuple<Real, int, Vector3r> smallestDistanceLoop(
	        vector<int> // indices of nodes to consider of "S" (with the special case of an empty vector triggering full consideration)
	        ,
	        const shared_ptr<LevelSet>&,
	        const shared_ptr<LevelSet>& // both shapes
	        ,
	        Vector3r,
	        Quaternionr // describing "S" configuration
	        ,
	        Vector3r,
	        Quaternionr // describing "B" configuration
	        ,
	        Vector3r,
	        Vector3r // the min and max of the bounds overlap
	); // return the obtained smallest distance (without sign condition) in get<0> ; the idx of the corresponding node (wrt 1st argument) in get<1> ; its current position in global frame in get<2>. The latter two can be garbage quantity in case get<0> is infinite (which means no node was detected to be in contact)
	std::tuple<Real, int, Vector3r>
	     smallestDistanceFullLoop( // variant of the above that will take care once for all of defining a full list of node indices
                const shared_ptr<LevelSet>&,
                const shared_ptr<LevelSet>&,
                Vector3r,
                Quaternionr,
                Vector3r,
                Quaternionr,
                Vector3r,
                Vector3r);
	bool goSingleOrMulti(
	        bool single,
	        const shared_ptr<Shape>&,
	        const shared_ptr<Shape>&,
	        const State&,
	        const State&,
	        const bool&,
	        const shared_ptr<Interaction>&,
	        const Vector3r&);
	bool go(const shared_ptr<Shape>&,
	        const shared_ptr<Shape>&,
	        const State&,
	        const State&,
	        const Vector3r&,
	        const bool&,
	        const shared_ptr<Interaction>&)
	        override; // reminder: method signature is imposed by InteractionLoop.cpp and also somewhat inherited from template class FunctorWrapper
	bool
	goReverse(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override
	{
		LOG_ERROR(
		        "We ended up calling goReverse.. How is this possible for symmetric IgFunctor ? Anyway, we now have to code something"); /* nothing, such as in TTetraGeom, mixed examples elsewhere*/
		return false;
	};
	Real normalMismatch(const shared_ptr<Interaction>&);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ig2_LevelSet_LevelSet_ScGeom,IGeomFunctor,R"""(Creates or updates a :yref:`ScGeom` instance representing the contact of two convex :yref:`LevelSet`-shaped bodies after executing a master-slave algorithm that combines distance function $\phi$ (:yref:`LevelSet.distField`) with surface nodes $\vec{N}$ (:yref:`LevelSet.surfNodes`) [Duriez2021a]_ [Duriez2021b]_. Denoting $S$, resp. $B$, the smallest, resp. biggest, contacting body, $\vec{N_c}$ the surface node of $S$ with the greatest penetration depth into $B$ (its current position), $u_n$ the corresponding :yref:`overlap<ScGeom.penetrationDepth>`, $\vec{C}$ the :yref:`contact point<ScGeom.contactPoint>` and $\vec{n}$ the contact :yref:`normal<ScGeom.normal>`, we have:

* $u_n = - \phi_B(\vec{N_c})$
* $\vec{n} = \pm \vec{\nabla} \phi_S(\vec{N_c})$  chosen to be oriented from :yref:`1<Interaction.id1>` to :yref:`2<Interaction.id2>`
* $\vec{C} = \vec{N_c} - \dfrac{u_n}{2} \vec{n}$

.. note:: in case the two :yref:`LevelSet grids<LevelSet.lsGrid>` no longer overlap for a previously existing interaction, the above workflow does not apply and $u_n$ is assigned an infinite tensile value that should insure interaction removal in the same DEM iteration (for sure with :yref:`Law2_ScGeom_FrictPhys_CundallStrack`).

While :yref:`Ig2_LevelSet_LevelSet_MultiScGeom` accomodates non-convex cases, :yref:`Ig2_LevelSet_LevelSet_LSnodeGeom` might be a better choice, towards faster simulations, when handling convex bodies.
)"""
	,/*attrs*/,/*ctor*/
	,.def("normalMismatch",&Ig2_LevelSet_LevelSet_ScGeom::normalMismatch,(boost::python::arg("cont")),"Normal orientation mismatch (in radians) of given interaction *cont* if master and slave particles were to be exchanged. Does not support periodic boundary conditions at the moment."));
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(LevelSet, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(LevelSet, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_LevelSet_LevelSet_ScGeom);

class Ig2_LevelSet_LevelSet_LSnodeGeom : public Ig2_LevelSet_LevelSet_ScGeom {
public:
	bool go(const shared_ptr<Shape>&,
	        const shared_ptr<Shape>&,
	        const State&,
	        const State&,
	        const Vector3r&,
	        const bool&,
	        const shared_ptr<Interaction>&) override; // reminder: method signature is imposed by InteractionLoop.cpp
	// clang-format-off
	bool
	goReverse(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override
	{
		LOG_ERROR(
		        "We ended up calling goReverse.. How is this possible for symmetric IgFunctor ? Anyway, we now have to code something"); /* nothing, such as in TTetraGeom, mixed examples elsewhere*/
		return false;
	};
	YADE_CLASS_BASE_DOC(
	        Ig2_LevelSet_LevelSet_LSnodeGeom,
	        Ig2_LevelSet_LevelSet_ScGeom,
	        "Same as :yref:`Ig2_LevelSet_LevelSet_ScGeom` except for the additional information in :yref:`LSnodeGeom` which should enable faster "
	        "simulations, provided that :yref:`Law2_ScGeom_FrictPhys_CundallStrack.neverErase` is defined as True. This functor actually creates an "
	        "interaction geom, enabling the interaction to turn real, as soon as a closest node can be detected (even though there is not yet mechanical "
	        "contact). On the contrary, in order to avoid a continuous increase in the interaction list, it will take the responsibility to ask for "
	        "interaction removal as soon as node treatment is impossible (due to an excessive gap), even when the interaction was previously real (at "
	        "variance with general YADE design), prompting the need for setting :yref:`InteractionLoop.warnRoleIg2` = False to avoid unnecessary warning "
	        "messages in such cases, see :ysrc:`examples/levelSet/lsNodeGeom.py` for an example. Please also note that using this Functor is expected to bias (increase) unbalancedForce measurements as it artificially increases the number of real interactions, with a number of those carrying no force, making for a smaller mean average interaction force, hence a higher unbalanced force")
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(LevelSet, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(LevelSet, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_LevelSet_LevelSet_LSnodeGeom);

class Ig2_LevelSet_LevelSet_MultiScGeom : public Ig2_LevelSet_LevelSet_ScGeom {
public:
	bool go(const shared_ptr<Shape>&,
	        const shared_ptr<Shape>&,
	        const State&,
	        const State&,
	        const Vector3r&,
	        const bool&,
	        const shared_ptr<Interaction>&) override; // reminder: method signature is imposed by InteractionLoop.cpp
	// clang-format-off
	bool
	goReverse(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override
	{
		LOG_ERROR(
		        "We ended up calling goReverse.. How is this possible for symmetric IgFunctor ? Anyway, we now have to code something"); /* nothing, such as in TTetraGeom, mixed examples elsewhere*/
		return false;
	};
	YADE_CLASS_BASE_DOC(
	        Ig2_LevelSet_LevelSet_MultiScGeom,
	        Ig2_LevelSet_LevelSet_ScGeom,
	        "Multiple contact points version of :yref:`Ig2_LevelSet_LevelSet_ScGeom` for a :yref:`MultiScGeom` description of a contact between two "
	        "(non-convex typically) :yref:`LevelSet`-shaped bodies (with a :yref:`ScGeom` interaction at each contacting surface node). Does not support "
	        "periodic boundary conditions at the moment. It is designed to be used in combination with :yref:`MultiFrictPhys` for what concerns the "
	        ":yref:`interaction physics<Interaction.phys>` (which is here also touched by that Ig2 in some contrast with general YADE design, from a "
	        "developer point of view) [Duriez2023]_.");
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(LevelSet, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(LevelSet, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_LevelSet_LevelSet_MultiScGeom);

class Ig2_Box_LevelSet_ScGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override;
	bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c) override
	{
		c->swapOrder();
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	};
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Box_LevelSet_ScGeom,IGeomFunctor,"Creates or updates a :yref:`ScGeom` instance representing the intersection of one :yref:`LevelSet` body with one :yref:`Box` body. Normal is given by the box geometry while overlap and contact points are defined likewise to :yref:`Ig2_LevelSet_LevelSet_ScGeom`. Restricted to the case of Boxes for which local and global axes coincide, and with non zero thickness, and assuming the center of the level set body never enters into the box (ie excluding big overlaps). You may prefer using :yref:`Ig2_Wall_LevelSet_ScGeom`.");
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Box, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(Box, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_Box_LevelSet_ScGeom);

class Ig2_Wall_LevelSet_ScGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override;
	bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c) override
	{
		c->swapOrder();
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	};
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Wall_LevelSet_ScGeom,IGeomFunctor,"Creates or updates a :yref:`ScGeom` instance representing the intersection of one :yref:`LevelSet`-shaped body with one :yref:`Wall`-shaped body, where overlap is chosen to occur on the opposite wall side than the LevelSet body's center. :yref:`Contact normal<ScGeom.normal>` is given by the wall normal (relative orientation of wall wrt global axes is not supported) while :yref:`overlap<ScGeom.penetrationDepth>` and :yref:`contact points<ScGeom.contactPoint>` are defined likewise to :yref:`Ig2_LevelSet_LevelSet_ScGeom`.");
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Wall, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(Wall, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_Wall_LevelSet_ScGeom);

class Ig2_Wall_LevelSet_MultiScGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override; // reminder: method signature is imposed by InteractionLoop.cpp
	bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c) override
	{
		c->swapOrder();
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	};
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Wall_LevelSet_MultiScGeom,IGeomFunctor,"Creates or updates a :yref:`MultiScGeom` instance representing the multiple contact points interaction kinematics of one :yref:`LevelSet` body with one :yref:`Wall` body, extending :yref:`Ig2_Wall_LevelSet_ScGeom` to non-convex LevelSet-shaped bodies. Relative orientation of wall wrt global axes is again not supported. TODO: time cost could / should be improved (wrt Ig2_LevelSet_LevelSet_MultiScGeom; jduriez note see aor8* and aor9*)");
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Wall, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(Wall, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_Wall_LevelSet_MultiScGeom);

} // namespace yade
#endif // YADE_LS_DEM
