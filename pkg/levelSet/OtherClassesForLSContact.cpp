/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <pkg/levelSet/OtherClassesForLSContact.hpp>
#include <pkg/levelSet/ShopLS.hpp>

namespace yade {
YADE_PLUGIN((Bo1_LevelSet_Aabb)(MultiPhys)(MultiFrictPhys)(MultiViscElPhys)(Ip2_FrictMat_FrictMat_MultiFrictPhys)(Ip2_ViscElMat_ViscElMat_MultiViscElPhys));
CREATE_LOGGER(Bo1_LevelSet_Aabb);
CREATE_LOGGER(Ip2_FrictMat_FrictMat_MultiFrictPhys);
CREATE_LOGGER(Ip2_ViscElMat_ViscElMat_MultiViscElPhys);

void Bo1_LevelSet_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*)
{ //TODO: use Eigen Aligned Box.extend to avoid the 2*6 if below ?
	// NB: see BoundDispatcher::processBody() (called by BoundDispatcher::action(), called by InsertionSortCollider::action()) in pkg/common/Dispatching.cpp for the attributes used upon calling
	if (!bv) { bv = shared_ptr<Bound>(new Aabb); }
	Aabb*     aabb    = static_cast<Aabb*>(bv.get()); // no need to bother deleting that raw pointer: e.g. https://stackoverflow.com/q/53908753/9864634
	LevelSet* lsShape = static_cast<LevelSet*>(cm.get());
	Real      inf     = std::numeric_limits<Real>::infinity();
	// We compute the bounds from LevelSet->corners serving as an Aabb in local frame, and considering the transformation from that local frame
	// NB: it is useless to try to make something much similar to Bo1_Box_Aabb::go(), in case se3.position of our level set body would not be in the middle of the lsShape->corners box, contrary to Box bodies and their extents and "halfSize" in Bo1_Box_Aabb::go()

	if (!lsShape->corners.size()) { // kind of "first iteration"
		Vector3i nGP(lsShape->lsGrid->nGP);
		Vector3r gp;
		int      nGPx(nGP[0]), nGPy(nGP[1]), nGPz(nGP[2]);
		Real     xGP, yGP, zGP;
		Real     xMin(inf), xMax(-inf), yMin(inf), yMax(-inf), zMin(inf), zMax(-inf); // extrema values for the "inside" gridpoints
		// identifying the extrema coordinates within the surface using a nGridPoints^3 loop. Would be much better to do 3*nGridpoints^2 loops, if possible, searching for find_if(firstIterator to lsValAlongZAxis,lastIterator to lsValAlongZAxis, within surface) which is easy for the zAxis = lsGrid->gridPoint[xInd][yInd], but how accessing an yAxis ~ lsGrid->gridPoint[xInd][:][zInd] which does not exist in C++ ?
		for (int xInd = 0; xInd < nGPx; xInd++) {
			for (int yInd = 0; yInd < nGPy; yInd++) {
				for (int zInd = 0; zInd < nGPz; zInd++) {
					gp  = lsShape->lsGrid->gridPoint(xInd, yInd, zInd);
					xGP = gp[0];
					yGP = gp[1];
					zGP = gp[2];
					if (lsShape->distField[xInd][yInd][zInd] <= 0) {
						if (xGP < xMin) xMin = xGP;
						if (xGP > xMax) xMax = xGP;
						if (yGP < yMin) yMin = yGP;
						if (yGP > yMax) yMax = yGP;
						if (zGP < zMin) zMin = zGP;
						if (zGP > zMax) zMax = zGP;
					}
				}
			}
		}
		if ((xMin == xMax) or (yMin == yMax) or (zMin == zMax))
			LOG_WARN("One flat LevelSet body, as detected by shape.corners computation, was that expected ? (is the grid too coarse ?)");
		// right now, our *Min, *Max define a downwards-rounded Aabb (smaller than true surface), let's make it upwards-rounded below:
		Real g(lsShape->lsGrid->spacing);
		for (int iInd = 0; iInd < 2; iInd++) {
			for (int jInd = 0; jInd < 2; jInd++) {
				for (int kInd = 0; kInd < 2; kInd++)
					lsShape->corners.push_back(
					        Vector3r(iInd == 0 ? xMin - g : xMax + g, jInd == 0 ? yMin - g : yMax + g, kInd == 0 ? zMin - g : zMax + g));
			}
		}
	}
	if (lsShape->corners.size() != 8) LOG_ERROR("We have a LevelSet-shaped body with some shape.corners computed but not 8 of them !");
	std::array<Vector3r, 8> cornersCurrent; // current positions of corners (in global frame)
	for (int corner = 0; corner < 8; corner++)
		cornersCurrent[corner] = ShopLS::rigidMapping(lsShape->corners[corner], Vector3r::Zero(), se3.position, se3.orientation);
	// NB: corners should average to the origin. This is kind of tested through LevelSet.center in levelSetBody() Py function

	Real xCorner, yCorner, zCorner;                                           // the ones of cornersCurrent
	Real xMin(inf), xMax(-inf), yMin(inf), yMax(-inf), zMin(inf), zMax(-inf); // extrema values for the cornersCurrent
	for (int corner = 0; corner < 8; corner++) {
		xCorner = cornersCurrent[corner][0];
		yCorner = cornersCurrent[corner][1];
		zCorner = cornersCurrent[corner][2];
		if (xCorner < xMin) xMin = xCorner;
		if (xCorner > xMax) xMax = xCorner;
		if (yCorner < yMin) yMin = yCorner;
		if (yCorner > yMax) yMax = yCorner;
		if (zCorner < zMin) zMin = zCorner;
		if (zCorner > zMax) zMax = zCorner;
	}
	aabb->min = Vector3r(xMin, yMin, zMin);
	aabb->max = Vector3r(xMax, yMax, zMax);
}


void Ip2_FrictMat_FrictMat_MultiFrictPhys::go(const shared_ptr<Material>& mat1, const shared_ptr<Material>& mat2, const shared_ptr<Interaction>& interaction)
{
	//NB: we will execute that Ip2::go at every iteration, because new items in interaction->phys->contacts have to be touched here (e.g., assigned kn)
	shared_ptr<MultiFrictPhys> multiFrictPhysPtr(new MultiFrictPhys);
	if (interaction->phys->getClassName() == "MultiPhys"){ // this interaction has just been created (and already populated with a phys, of MultiPhys type), by the Ig2*MultiScGeom::go
		LOG_DEBUG("Interaction " << interaction->id1 << " - " << interaction->id2 << " seen for the first time in the Ip2" << std::endl);
		// we will backport interaction->phys->contacts and ->nodesIds;
		shared_ptr<MultiPhys> iPhysPtrAsMultiPhys(YADE_PTR_CAST<MultiPhys>(interaction->phys));
		multiFrictPhysPtr->contacts = iPhysPtrAsMultiPhys->contacts;
		multiFrictPhysPtr->nodesIds = iPhysPtrAsMultiPhys->nodesIds;
	}
	else if(interaction->phys->getClassName() == "MultiFrictPhys") { // it has existed for a while: i->phys is already of MultiFrictPhys type
		LOG_DEBUG("Interaction " << interaction->id1 << " - " << interaction->id2 << " had already been Ip2-handled in the past, its phys is of " << interaction->phys->getClassName() << " at least" << std::endl);
		multiFrictPhysPtr = YADE_PTR_CAST<MultiFrictPhys>(interaction->phys);
	}
	else{
		LOG_ERROR("Interaction " << interaction->id1 << " - " << interaction->id2 << " has a phys of " << interaction->phys->getClassName() << ", that is clearly not expected ! (only MultiPhys or MultiFrictPhys are)" << std::endl);
	}
	LOG_DEBUG(
		"YADE_PTR_DYN_CAST<MultiFrictPhys>(interaction->phys) = "
		<< YADE_PTR_DYN_CAST<MultiFrictPhys>(interaction->phys)
		<< " (0 would be a problem)"); // the dynamic cast will never be executed here unless LOG_DEBUG level is actually effective
	// Computing once for all (outside of below loop) the to-be-assigned contact friction from the Materials frictionAngle(s):
	const shared_ptr<FrictMat>& fmat1 = YADE_PTR_CAST<FrictMat>(mat1);
	const shared_ptr<FrictMat>& fmat2 = YADE_PTR_CAST<FrictMat>(mat2);
	Real frictAngle( std::min(fmat1->frictionAngle, fmat2->frictionAngle) );
	// Looping over contacting nodes to assign mechanical properties into individual contacts items:
	LOG_DEBUG("About to looping over contacting nodes");
	for (unsigned int idx = 0; idx < multiFrictPhysPtr->contacts.size(); idx++) {
		if (multiFrictPhysPtr->contacts[idx]->getClassName() == "IPhys") { // this contacts[idx] item has just been created in Ig2::go
			LOG_DEBUG("Contact at node " << multiFrictPhysPtr->nodesIds[idx] << " is detected by the Ip2 as being just created" << std::endl);
			shared_ptr<FrictPhys> frictPhysPtr(new FrictPhys);
			// Direct assignment of values of stiffnesses:
			frictPhysPtr->kn                     = kn;
			frictPhysPtr->ks                     = ks;
			// Assignment of (tangent of) friction angle:
			frictPhysPtr->tangensOfFrictionAngle = std::tan(frictAngle);
			multiFrictPhysPtr->contacts[idx] = frictPhysPtr;
		}
		else { // this contacts[idx] item shall have existed for a while, there is nothing to do
			LOG_DEBUG("Contact at node " << multiFrictPhysPtr->nodesIds[idx] << " is detected by the Ip2 as existing for a while" << std::endl);
		}
	}
	LOG_DEBUG("Nodes loop done");
	interaction->phys = multiFrictPhysPtr;
};


void Ip2_ViscElMat_ViscElMat_MultiViscElPhys::go(const shared_ptr<Material>& mat1, const shared_ptr<Material>& mat2, const shared_ptr<Interaction>& interaction)
{
	//NB: we will execute that Ip2::go at every iteration, because new items in interaction->phys->contacts have to be touched here (e.g., assigned kn)
	shared_ptr<MultiViscElPhys> multiViscElPhysPtr(new MultiViscElPhys);
	if (interaction->phys->getClassName() == "MultiPhys"){ // this interaction has just been created (and already populated with a phys, of MultiPhys type), by the Ig2*MultiScGeom::go
		LOG_DEBUG("Interaction " << interaction->id1 << " - " << interaction->id2 << " seen for the first time in the Ip2" << std::endl);
		// we will backport interaction->phys->contacts and ->nodesIds;
		shared_ptr<MultiPhys> iPhysPtrAsMultiPhys(YADE_PTR_CAST<MultiPhys>(interaction->phys));
		multiViscElPhysPtr->contacts = iPhysPtrAsMultiPhys->contacts;
		multiViscElPhysPtr->nodesIds = iPhysPtrAsMultiPhys->nodesIds;
	}
	else if(interaction->phys->getClassName() == "MultiViscElPhys") { // it has existed for a while: i->phys is already of MultiViscElPhys type
		LOG_DEBUG("Interaction " << interaction->id1 << " - " << interaction->id2 << " had already been Ip2-handled in the past, its phys is of " << interaction->phys->getClassName() << " at least" << std::endl);
		multiViscElPhysPtr = YADE_PTR_CAST<MultiViscElPhys>(interaction->phys);
	}
	else{
		LOG_ERROR("Interaction " << interaction->id1 << " - " << interaction->id2 << " has a phys of " << interaction->phys->getClassName() << ", that is clearly not expected ! (only MultiPhys or MultiViscElPhys are)" << std::endl);
	}
	LOG_DEBUG(
		"YADE_PTR_DYN_CAST<MultiViscElPhys>(interaction->phys) = "
		<< YADE_PTR_DYN_CAST<MultiViscElPhys>(interaction->phys)
		<< " (0 would be a problem)"); // the dynamic cast will never be executed here unless LOG_DEBUG level is actually effective
	// Computing once for all (outside of below loop) the to-be-assigned contact friction from the Materials frictionAngle(s):
	
	const shared_ptr<ViscElMat>& fmat1 = YADE_PTR_CAST<ViscElMat>(mat1);
	const shared_ptr<ViscElMat>& fmat2 = YADE_PTR_CAST<ViscElMat>(mat2);
	
	Real frictAngle( std::min(fmat1->frictionAngle, fmat2->frictionAngle) );

	Real                              kna             = fmat1->kn;
	Real                              knb             = fmat2->kn;

	Real                              ksa             = fmat1->ks;
	Real                              ksb             = fmat2->ks;
	
	Real                              ena             = fmat1->en;
	Real                              enb             = fmat2->en;
	
	Real                              esa             = fmat1->et; // ViscElMat notes the tangential coefficient of restitution as "et" - here we use "es" where "s" stands for "shear"
	Real                              esb             = fmat2->et; 
	
	Real mass1 = Body::byId(interaction->getId1())->state->mass;
	Real mass2 = Body::byId(interaction->getId2())->state->mass;
	if (mass1 == 0.0 and mass2 > 0.0) {
		mass1 = mass2;
	} else if (mass2 == 0.0 and mass1 > 0.0) {
		mass2 = mass1;
	}
	
	// See [Pournin2001, just below equation (19)]
	const Real massR = mass1 * mass2 / (mass1 + mass2);

	Real Kn = (kn) ? (*kn)(mat1->id, mat2->id) : 2. * kna * knb / (kna + knb);
	Real Ks = (ks) ? (*ks)(mat1->id, mat2->id) : 2. * ksa * ksb / (ksa + ksb);
	
	Real En = (en) ? (*en)(mat1->id, mat2->id) : 2. * ena * enb / (ena + enb);
	Real Es = (et) ? (*et)(mat1->id, mat2->id) : 2. * esa * esb / (esa + esb);
	
	Real cn = 2.0 * find_cn_from_en(En, massR, contactParameterCalculation(Kn, Kn), interaction);
	Real cs = 2.0 * find_cn_from_en(Es, massR, contactParameterCalculation(Ks, Ks), interaction);
	
	// Looping over contacting nodes to assign mechanical properties into individual contacts items:
	LOG_DEBUG("About to looping over contacting nodes");
	for (unsigned int idx = 0; idx < multiViscElPhysPtr->contacts.size(); idx++) {
		if (multiViscElPhysPtr->contacts[idx]->getClassName() == "IPhys") { // this contacts[idx] item has just been created in Ig2::go
			LOG_DEBUG("Contact at node " << multiViscElPhysPtr->nodesIds[idx] << " is detected by the Ip2 as being just created" << std::endl);
			shared_ptr<ViscElPhys> viscElPhysPtr(new ViscElPhys);
			// Direct assignment of values of stiffnesses:
			viscElPhysPtr->kn                     = Kn;
			viscElPhysPtr->ks                     = Ks;
			// Direct assignment of values of viscous damping coefficients:
			viscElPhysPtr->cn                     = cn;
			viscElPhysPtr->cs                     = cs;
			// Assignment of (tangent of) friction angle:
			viscElPhysPtr->tangensOfFrictionAngle = std::tan(frictAngle);
			multiViscElPhysPtr->contacts[idx] = viscElPhysPtr;
		}
		else { // this contacts[idx] item shall have existed for a while, there is nothing to do
			LOG_DEBUG("Contact at node " << multiViscElPhysPtr->nodesIds[idx] << " is detected by the Ip2 as existing for a while" << std::endl);
		}
	}
	LOG_DEBUG("Nodes loop done");
	interaction->phys = multiViscElPhysPtr;
};

} // namespace yade
#endif //YADE_LS_DEM
