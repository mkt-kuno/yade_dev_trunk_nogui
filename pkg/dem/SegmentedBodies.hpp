// 2024 © Karol Brzeziński <karol.brze@gmail.com>
// 2024 © Vasileios Angelidakis <angelidakis@qub.ac.uk>
#pragma once
#include <core/Shape.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/CohesiveFrictionalContactLaw.hpp>

namespace yade {


/*********************************************************************************
*
* S T A T E
*
*********************************************************************************/

class SegmentedState : public State {
public:
        struct Segment {
            int segmentPosition; //segment position in state
            Real thickness;
            Real volume;
            Real area;
        };
        struct SegmentContact { // structure representing geometry of two neighbor ("contacting") segments of the same sphere
            int segmentPosition1; //segment position in state
            int segmentPosition2; //segment position in state
            Real distance;
            Real edgeLen;
        };
        vector<SegmentContact> segmentContacts; //stores information about segment "interactions" within one sphere (segments "contacting" to each other)
        int thetaResolution = -1;// positive only for spheres - this is also the way to know the shape of the particle
        int phiResolution = -1;
        Vector3r facetRefNormal; // If we store facet normal during initialization, we will be able to know it without accessing the shape 
        Segment getSegmentFromContactPoint(Vector3r midPoint_);
        void initializeSegmentContacts();
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(SegmentedState,State,"State information about body coating map used by :yref:`CohFrictMatSeg`. It allows spheres or facets to store information about the coating thickness. The coating thickness for the facets is stored in two fields (positive and negative side of the facets). In the case of spheres, the surface of sphere is divided into segments of equal area (controlled by :yref:`SegmentedStateUpdater::thetaResolution` and :yref:`SegmentedStateUpdater.phiResolution` parameters).",
		((vector<Real>,coatingThickness,,,"Vector storing information about coating thickness in each segment. It is assumed that the linear unit of coating thickness is 1000 times smaller than the unit of radius (e.g., if the simulation basic unit is meter, then coating thickness is provided in milimeters). It is one-dimensional, hence for easier intepretation can be reshaped based on the information about thetaResolution and phiResolution. "))
		((vector<Real>,coatingVolume,,,"Vector storing information about coating volume in each segment (unit is $1000^3$ smaller than volume unit of the simulation, see :yref:`SegmentedState.coatingThickness` for explanation). It is one-dimensional, hence for easier intepretation can be reshaped based on the information about thetaResolution and phiResolution."))
		((Real,segmentArea,0,Attr::readonly,"Area of one segment of the body  (unit is $1000^2$ smaller than surface unit of the simulation, see :yref:`SegmentedState::coatingThickness` for explanation)."))
		((bool,needsInit,true,Attr::readonly,"Tells :yref:`SegmentedStateUpdater` whether thickness was initialized for given body."))
		((bool,allowMatExchange,true,,"If true, material can be exchanged between bodies.This condition needs to be true for material exchange, but also body state needs to be initialized with :yref:`SegmentedStateUpdater`.")),
		/*ctor*/ createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(SegmentedState, State);
};
REGISTER_SERIALIZABLE(SegmentedState);

/*********************************************************************************
*
* S T A T E   U P D A T E R 
*
*********************************************************************************/

/* SegmentedStateUpdater is designed to be setter of thickness, but we inherited it from PeriodicEngine in case one would like to propose more sothetasticated law for exchanging coating material between two bodies (e.g. time-dependent). */

class SegmentedStateUpdater : public PeriodicEngine {
public:

	void action() override; //
	void initialize(Body::id_t bId, Real thickness_ = 0);
	void setThicknessToSegmentRange(Body::id_t bId, int thetaMin_, int thetaMax_, int phiMin_, int phiMax_, Real thickness_ = 0); // general function
	void setSegmentThickness(
	        Body::id_t bId, int segmentPos_, Real thickness_); // directly sets thickness of segment (segmentPos in coatingThickness needs to be known)
	Real computeVolumeExchange(Real thickness1_, Real thickness2_, Real volume1_, Real volume2_, Real exchangeRate_); // returns volume change in segments (subtract from first segment and add to second segment)
	void setThicknessToSphereSeg(Body::id_t bId, int theta_, int phi_, Real thickness_ = 0);
	void setThicknessToFacetSide(Body::id_t bId, int side_, Real thickness_ = 0);
	void setThicknessToSpheres(vector<Body::id_t> bIds_, int thetaMin_, int thetaMax_, int phiMin_, int phiMax_, Real thickness_ = 0);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(SegmentedStateUpdater,PeriodicEngine,"This class provides methods for setting the coating thickness of bodies and performing the coating material transfer (both interparticle and intraparticle).",
		((int,thetaResolution,2,,"In how many regions a coating map of sphere should be divided vertically (theta is polar angle ranging from 0 to 180 degrees)."))
		((int,phiResolution,4,,"In how many regions a coating map of sphere should be divided horizontally (phi is azimuthal angle ranging from 0 to 360 degrees)."))
		// I will try to better comment the below listed coefficients in the code of the functions and we will decide how to document it.
		((Real,interParticleExchangeCoeff,1,,"Rate of coating material exchange rate coefficient between two contacting particles."))
		((Real,intraParticleExchangeCoeff,1,,"Rate of coating material exchange rate coefficient between two segments of the same particle."))
		((int,activateWettability,false,,"If true, particles are wettable (coating material is distributed over particles with time). It works only for spherical particles. For clumps of spheres it is recommended to set false."))
		((Real,prevTime,0,Attr::readonly,"Virtual time of the last run. Compared to the :yref:`PeriodicEngine::virtLast`, it is updated after the action of the engine."))//The difference won't be visible in the console, because no-one access such variable with Puthon during execution. However, putting it into the macro ensures proper save and load (https://gitlab.com/yade-dev/trunk/-/merge_requests/1143#note_2973402989)
		((int,reInitCohesion,1,,"This parameter is introduced to reinitialize the cohesion after two bodies exchange the coating material. That means it sets :yref:`CohFrictPhys::initCohesion` true for the interaction. It is later determined by Ip2, whether the interaction should be cohesive (based on the material properties). Default value is 1 which means that reinitialization works only for interactions that are not yet cohesive (:yref:`CohFrictPhys::cohesionBroken` is true). Set 2, in order to reinitialize cohesion on every run of :yref:`SegmentedStateUpdater`."))
		,//this comma after all the params
		/*ctor*/,
	.def("initialize",&SegmentedStateUpdater::initialize,"Initializes :yref:`SegmentedState::coatingThickness`, so it has a proper length and next, sets the same thickness to all the segments. Not initialized bodies cannot exchange the coating material.")
	.def("setThicknessToSphereSeg",&SegmentedStateUpdater::setThicknessToSphereSeg,"Set thickness to the segment of sphere at position given by phi and theta.")
	.def("setThicknessToFacetSide",&SegmentedStateUpdater::setThicknessToFacetSide,"Set thickness to facet side (0 to negative side, 1 to positive side, according to facet normal).") // maybe -1 would be more appropriate than zero, but I it would be counter-intuitive when compared with segment addressing
	.def("setThicknessToSpheres",&SegmentedStateUpdater::setThicknessToSpheres,"For each sphere from provided list (of Ids) set thickness at some range of segments specified by minimum and maximum theta and phi.")
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SegmentedStateUpdater);


/*********************************************************************************
*
* S E G M E N T E D   M A T   S P R I N K L E R 
*
*********************************************************************************/

/* A periodic engine distributing new material of particles having CohFrictMatSeg. Determines the particles to cover by ray tracing algorithm. */

class SegmentedMatSprinkler : public PeriodicEngine {
private:
	bool needsInit = true;
	void init();
    vector<Vector3r> cornerVersors; // corner versors required to define limiting planes of the volume potentially within the reach of the sprinkler.
    vector<Body::id_t> roiSpheres; // Spheres in the region of interest. Storing preselected spheres that can be reached by the ray.
    struct rayToSphereInfo { // structure to pass information information about versor / ray hitting a sphere
        bool bodyFound;
        Real minB;
        Body::id_t bId;
    };    
public:
	void action() override; //
	rayToSphereInfo oneLineAction(Vector3r lineVersor_); // performs action for one line (searches for sphere and applies to the closest segment)
	/// functions exposed to Pyhon below
	void setVersors(Quaternionr ori_);//
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(SegmentedMatSprinkler,PeriodicEngine,"A periodic engine distributing coating material of particles having :yref:`CohFrictMatSeg`. Determines the particles to cover by ray tracing algorithm.",
		((Vector3r,pos,Vector3r(0,0,0),,"Origin of the sprinkling source."))
		((Vector3r,refDir,Vector3r(0,0,-1),Attr::readonly,"Reference direction of SegmentedMatSprinkler (center of the stream). It should not be directly modified. When the direction of sprinkler is rotated in setVersors function, this direction is taken as reference.")) 
		((Real,feedRate,1,,"Volume of the material distributed per second (unit is $1000^3$ smaller than volume unit of the simulation, see :yref:`SegmentedState.coatingThickness` for explanation).")) 
		((Real,alpha,45,,"Angle (degrees) defining width of the stream in the first direction.")) 
		((Real,beta,45,,"Angle (degrees) defining width of the stream in the second direction."))
		((int,alphaResolution,50,,"How many searching lines scan one cross-section in first direction"))
		((int,betaResolution,50,,"How many searching lines scan one cross-section in first direction"))	
		((vector<Vector3r>,versors,,,"Versors of the searching lines."))
		((Real,prevTime,0,Attr::readonly,"Virtual time of the last run. Compared to the :yref:`PeriodicEngine::virtLast`, it is updated after the action of the engine."))//The difference won't be visible in the console, because no-one access such variable with Puthon during execution. However, putting it into the macro ensures proper save and load (https://gitlab.com/yade-dev/trunk/-/merge_requests/1143#note_2973402989)
		((Real,releasedVolume,0,Attr::readonly,"Volume released by the sprinkler since the begining of its work.")) 
		((Real,lostVolume,0,Attr::readonly,"Released volume that did not reach any target.")) 
		,/*ctor*/ needsInit=true,
	.def("setVersors",&SegmentedMatSprinkler::setVersors,"Sets versors of lines searching for the surfaces of closest spheres. Should always be used after changing any parameters of SegmentedMatSprinkler other than feedRate.") 
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SegmentedMatSprinkler);

/*********************************************************************************
*
* M A T E R I A L
*
*********************************************************************************/
class CohFrictMatSeg : public CohFrictMat {
public:
	shared_ptr<State> newAssocState() const override { return shared_ptr<State>(new SegmentedState); }
	bool              stateTypeOk(State* s) const override { return (bool)dynamic_cast<SegmentedState*>(s); }

	virtual ~CohFrictMatSeg() {};
	/// Serialization
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CohFrictMatSeg,CohFrictMat,"Material extending :yref:`CohFrictMat` to use dual material properties with segmented bodies (transition between properties depends on the coating thickness). Should be used with :yref:`Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys` to benefit from dual material capabilites. The other functors can be the same as for :yref:`CohFrictMat` (e.g. with contact law :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`). Besides primary set of material parameter of :yref:`CohFrictMat`, some are doubled and called secondary (e.g. secondaryYoung, secondaryNormalCohesion etc.). When this material is used with spheres or facets, a special state (:yref:`SegmentedState`) is assigned to them, allowing to store information about coating thickness.",
		//((bool,isCohesive,true,,"Whether this body can form possibly cohesive interactions (if true and depending on other parameters such as :yref:`Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys`).")) // !!! This description would be better (new refference), but now `isCohesive` is inherited from CohFrictMat - I am not sure how overshadowing would work. TODO
		((Real,secondaryNormalCohesion,-1,,"Tensile strength (secondary), homogeneous to a pressure. If negative the normal force is purely elastic."))
		((Real,secondaryShearCohesion,-1,,"Shear strength (secondary), homogeneous to a pressure. If negative the shear force is purely elastic."))
		((Real,secondaryYoung,1e9,,"elastic modulus (secondary) [Pa]. It has different meanings depending on the Ip functor."))
		((Real,secondaryPoisson,.25,,"Poisson's ratio (secondary) or the ratio between shear and normal stiffness [-]. It has different meanings depending on the Ip functor. "))
		((Real,secondaryFrictionAngle,.5,,"Friction angle (secondary)."))
		((Real,minThickness,0.0,,"Value of thickness parameter, below which only primary material parameters are taken into account."))
		((Real,maxThickness,1.0,,"Value of thickness parameter, above which only secondary material parameters are taken into account."))
		,
		createIndex();
		);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(CohFrictMatSeg, CohFrictMat);
};

REGISTER_SERIALIZABLE(CohFrictMatSeg);


/*********************************************************************************
*
* I P 2
*
*********************************************************************************/

class Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys : public IPhysFunctor {
public:
	// new structure for storing blended properties of the material based on the thickness
        struct BlendedProp {
            Real blendedFrictionAngle; //segment position in state
            Real blendedYoung;
            Real blendedPoisson;
            Real blendedNormalCohesion;
            Real blendedShearCohesion;
        };
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	int  cohesionDefinitionIteration;
	// assign cohesion, return false only if typecasting fails
	bool setCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, BlendedProp matProp1, BlendedProp matProp2);
	/* this is a shorthand and preserving it would be very complicated, while similar effect can be obtained via Puthon console. Hence, I propose to get rid of this function.
	// the python call will go through dyn_cast since the passed pointer is null
	void pySetCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, bool resetDisp)
	{
		bool assigned = setCohesion(interaction, cohesive, nullptr);
		if (assigned and resetDisp)
			YADE_CAST<CohFrictPhys*>(interaction->phys.get())->unp = YADE_CAST<ScGeom6D*>(interaction->geom.get())->penetrationDepth;
	}*/
	// method for blending properties of the material based on the thickness
        BlendedProp blendMatProp(Body::id_t bId, int segmentPosition_); // returns 'blended' primary and secondary material properties of the body based on the coating thickness
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys,IPhysFunctor,
		"This is a version of Ip2 to be used with :yref:`CohFrictMatSeg` (material with coating thickness-dependent properties). It generates cohesive-frictional interactions with moments, used in the contact law :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`. For details regarding the computation of cohesion, twisting, bending, etc. please refer to :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys`.",
		((bool,setCohesionNow,false,,"If true, assign cohesion to all existing contacts in current time-step. The flag is turned false automatically, so that assignment is done in the current timestep only."))
		((bool,setCohesionOnNewContacts,false,,"If true, assign cohesion at all new contacts. If false, only existing contacts can be cohesive (also see :yref:`Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::setCohesionNow`), and new contacts are only frictional."))
		// !!! I am not sure whether we should provide machmakers. I think it will be to complex, and confusing. TODO
		/*((shared_ptr<MatchMaker>,normalCohesion,,,"Instance of :yref:`MatchMaker` determining tensile strength"))
		((shared_ptr<MatchMaker>,shearCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the shear strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`)"))
		((shared_ptr<MatchMaker>,rollingCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the rolling strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`). The default is $\\frac{r}{4}R_t$ with $R_t$ the shear strength (inspired by stress in beams with circular cross-section)."))
		((shared_ptr<MatchMaker>,twistingCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the twisting strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`). The default is $\\frac{r}{2}R_s$ with $R_s$ the shear strength (inspired by stress in beams with circular cross-section)."))
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))*/
		,//ctor
		cohesionDefinitionIteration = -1;
		,//py
		//.def("setCohesion",&Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::pySetCohesion,(boost::python::arg("interaction"),boost::python::arg("cohesive"),boost::python::arg("resetDisp")),"Bond or un-bond an interaction with cohesion.\n\n  When ``True``, the resulting state is the same as what's obtained by executing an :yref:`InteractionLoop` with the functor's :yref:`setCohesionOnNewContacts<Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesionOnNewContacts>` or the interaction's :yref:`CohFrictPhys::initCohesion` ``True``. It will use the matchmakers if defined. The only difference is that calling this function explicitly will make the contact cohesive even if not both materials have :yref:`CohFrictMat::isCohesive`=``True``.\n\n When ``False``, the resulting state is the same as after breaking a fragile interaction. If `resetDisp` is ``True``, the current distance is taken as the reference for computing normal displacement and normal force.")
		);
	// clang-format on
	FUNCTOR2D(CohFrictMatSeg, CohFrictMatSeg);
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys);


}; // namespace yade
