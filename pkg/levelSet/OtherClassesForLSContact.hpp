/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/levelSet/LevelSet.hpp>

namespace yade {
class Bo1_LevelSet_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*) override;
	FUNCTOR1D(LevelSet);
	YADE_CLASS_BASE_DOC(Bo1_LevelSet_Aabb, BoundFunctor, "Creates/updates an :yref:`Aabb` of a :yref:`LevelSet`");
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Bo1_LevelSet_Aabb);

class MultiPhys: public IPhys {
	public:
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(MultiPhys,IPhys,"Describes the physical part of an interaction with multiple contact points, e.g., between two :yref:`LevelSet` bodies, as a set of :yref:`IPhys` items. This class is actually not intended to be used directly but to serve as a common ancestor for children classes such as :yref:`MultiFrictPhys`.",
	((vector< shared_ptr<IPhys> >,contacts,,,"The actual list of :yref:`IPhys` (for the mother class, actually obtained types can be different for derived classes such as :yref:`MultiFrictPhys` which will include :yref:`FrictPhys` instances) items corresponding to the different contact points."))
	((vector< int >,nodesIds,,,"The physics counterpart of :yref:`MultiScGeom.nodesIds` (both should be equal by design).")) // let us keep both for verification purpose
	,
	createIndex(); // this class will enter InteractionLoop dispatch, we need a create_index() here, and a REGISTER_*_INDEX below (https://yade-dem.org/doc/prog.html#indexing-dispatch-types)
	);
	// clang-format on
	REGISTER_CLASS_INDEX(MultiPhys, IPhys); // see createIndex() remark
};
REGISTER_SERIALIZABLE(MultiPhys);

class MultiFrictPhys : public MultiPhys {
public:
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(MultiFrictPhys,MultiPhys,"Describes the physical part of an interaction with multiple frictional contact points, e.g., between two :yref:`LevelSet` bodies, through a set of :yref:`FrictPhys` instances in :yref:`contacts<MultiFrictPhys.contacts>`. To combine with :yref:`MultiScGeom` and associated classes.",
	,
	createIndex(); // this class will enter InteractionLoop dispatch, we need a create_index() here, and a REGISTER_*_INDEX below (https://yade-dem.org/doc/prog.html#indexing-dispatch-types)
	);
	// clang-format on
	REGISTER_CLASS_INDEX(MultiFrictPhys, MultiPhys); // see createIndex() remark
};
REGISTER_SERIALIZABLE(MultiFrictPhys);

class Ip2_FrictMat_FrictMat_MultiFrictPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	FUNCTOR2D(FrictMat, FrictMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ip2_FrictMat_FrictMat_MultiFrictPhys,IPhysFunctor,"Handles the :yref:`MultiFrictPhys` physical description of the contact between two :yref:`FrictMats<FrictMat>`. Contact stiffnesses (for every :yref:`contact<MultiFrictPhys.contacts>`) are directly assigned from below attributes, independent of FrictMat properties. Contact friction angle is taken as the minimum of the 2 material friction angles (:yref:`FrictMat.frictionAngle`).",
		((Real,kn,0,,"Chosen value for :yref:`MultiFrictPhys.kn`"))
		((Real,ks,0,,"Chosen value for :yref:`MultiFrictPhys.ks`"))
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_MultiFrictPhys);

} // namespace yade
#endif // YADE_LS_DEM
