/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <pkg/common/Sphere.hpp>
#include <pkg/levelSet/LevelSetIGeom.hpp>
#include <pkg/levelSet/LevelSetIg2.hpp>
#include <pkg/levelSet/OtherClassesForLSContact.hpp>
#include <pkg/levelSet/ShopLS.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade {
YADE_PLUGIN((Ig2_Sphere_LevelSet_ScGeom)(Ig2_Box_LevelSet_ScGeom)(Ig2_Wall_LevelSet_ScGeom)(Ig2_Wall_LevelSet_MultiScGeom)(Ig2_LevelSet_LevelSet_ScGeom)(Ig2_LevelSet_LevelSet_LSnodeGeom)(
        Ig2_LevelSet_LevelSet_MultiScGeom));
CREATE_LOGGER(Ig2_Sphere_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_Box_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_LevelSet_LevelSet_LSnodeGeom);
CREATE_LOGGER(Ig2_LevelSet_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_LevelSet_LevelSet_MultiScGeom);
CREATE_LOGGER(Ig2_Wall_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_Wall_LevelSet_MultiScGeom);

bool Ig2_Sphere_LevelSet_ScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	shared_ptr<Sphere>   sphereSh = YADE_PTR_CAST<Sphere>(shape1);
	shared_ptr<LevelSet> lsSh     = YADE_PTR_CAST<LevelSet>(shape2);
	Vector3r             lsPos(state2.pos + shift2), spherePos(state1.pos);

	// Orientation of the level set (sphere orientation doesn't matter).
	// ori = rotation from reference configuration (local axes) to current one (global axes)
	// ori.conjugate() from the current configuration (global axes) to the reference one (local axes)
	Quaternionr rotLS(state2.ori), rotConjLS(state2.ori.conjugate());

	// Position of the sphere relative to the level set, then rotate to local coordinates of the LS grid
	Vector3r spherePosLocal = rotConjLS * (lsPos - spherePos);

	// Evaluate the level set and the normal, sphere centre may be outside so set true.
	Real     lev         = lsSh->distance(spherePosLocal, true);
	Vector3r normalLocal = lsSh->normal(spherePosLocal, true);

	// Compute the maximum overlap
	Real rad(sphereSh->radius); // The sphere radius
	Real maxOverlap = rad - lev;

	// Bring normal back to local coordinates, adjust to point from body 1 to body 2
	Vector3r normal = (rotLS * normalLocal).normalized();

	// Determine contact point. Middle of overlapping volumes, as usual.
	Vector3r contactPoint = spherePos + (rad - maxOverlap / 2.) * normal;

	if (maxOverlap < 0 && !c->isReal() && !force)
		return false; // We won't create the interaction in this case (but it is not our job here to delete it in case it already exists).

	shared_ptr<ScGeom> geomPtr;
	bool               isNew = !c->geom;

	if (isNew) {
		geomPtr = shared_ptr<ScGeom>(new ScGeom());
		c->geom = geomPtr;
	} else
		geomPtr = YADE_PTR_CAST<ScGeom>(c->geom);
	geomPtr->doIg2Work(
	        contactPoint,
	        maxOverlap, // Doesn't work for maxOverlap > R, like how pretty much any contact law doesn't.
	        rad,        // Inconsequential value since the contact point to centre distance is used for the torque anyway.
	        rad,        // Inconsequential value since the contact point to centre distance is used for the torque anyway.
	        state1,
	        state2,
	        scene,
	        c,
	        normal,
	        shift2,
	        isNew,
	        false);
	return true;
}

bool Ig2_Box_LevelSet_ScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	// 1.1 Preliminary declarations
	shared_ptr<Box>      boxSh = YADE_PTR_CAST<Box>(shape1);
	shared_ptr<LevelSet> lsSh  = YADE_PTR_CAST<LevelSet>(shape2);
	Vector3r             lsCenter(state2.pos + shift2), boxCenter(state1.pos);
	Vector3r             axisContact(Vector3r::Zero()); // I will need to .dot this vector with Vector3r => Vector3r (and not Vector3i) as well

	// 1.2 Checking the "contact direction", adopting the two following hypothesis:
	// - H1: lsCenter is outside of the box
	// - H2: projecting lsCenter along the contact direction towards the box will make lsCenter fall onto the  box, and not alongside. This is not always fulfilled, but will always be during triaxial box setups, for instance (as long as the box is closed).
	for (int axis = 0; axis < 3; axis++) {
		if (math::abs(lsCenter[axis] - boxCenter[axis]) > boxSh->extents[axis]) // TODO: account for a change in orientation of the box ?
			axisContact[axis] = 1; // TODO: once I'm sure this will happen only once, I could stop the loop here...
	}
	if (axisContact.norm() != 1) {
		LOG_ERROR(
		        "Problem while determining contact direction for "
		        << c->id1 << "-" << c->id2 << " : we got " << axisContact
		        << ". (0 0 0) means the LevelSet'd body has its center inside the box, which is not supported. Indeed, center = " << lsCenter
		        << " while boxCenter = " << boxCenter << " and extents = " << boxSh->extents);
	}
	LOG_DEBUG("axisContact = " << axisContact);
	Real boxC(axisContact.dot(boxCenter)), boxE(boxSh->extents.dot(axisContact)), lsC(lsCenter.dot(axisContact));

	// 2.1. Preliminary declarations for the surface nodes loop
	// clang-format off
	const int nNodes(lsSh->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");
	vector<Vector3r> lsNodes; // current positions of the boundary nodes (some of of those, at least) for the level set body. See for loop below
	lsNodes.reserve(nNodes); // nNodes will be a maximum size, reserve() is appropriate, not resize() (see also https://github.com/isocpp/CppCoreGuidelines/issues/493)
	vector<Real> distList; // will include all distance values from the node(s)On1 to shape2. It is a std::vector because we do not know beforehand the number of elements in this "list
	distList.reserve(nNodes); // nNodes might be a maximum size, reserve() is appropriate, not resize() (see also https://github.com/isocpp/CppCoreGuidelines/issues/493)
	vector<int> indicesNodes; // distList will include distance values corresponding to nodes e.g. 2, 5, 7 only out of 10 nodes among surfNodes. This indicesNodes vector will store these 2,5,7 indices
	indicesNodes.reserve(nNodes);
	// NB: it might actually be somewhat faster to not use these vectors and just compare new node with previous node, as done in Ig2_LevelSet_LevelSet*
	// I do not think it is critical for present Ig2_Box_LevelSet_ScGeom
	Real distToNode; // one distance value, for one node
	Real xNode;      // current position (on the axisContact of interest, can be something else than x-axis..) of the level set boundary node
	// clang-format on

	// 2.2. Actual loop over surface nodes
	for (int node = 0; node < nNodes; node++) {
		lsNodes[node] = ShopLS::rigidMapping(lsSh->surfNodes[node], Vector3r::Zero(), lsCenter, state2.ori);
		xNode         = lsNodes[node].dot(axisContact);
		if (xNode < boxC - boxE || xNode > boxC + boxE) continue;
		distToNode = (lsC > boxC ? boxC + boxE - xNode : xNode - (boxC - boxE));
		if (distToNode < 0) LOG_ERROR("Unexpected case ! We obtained " << distToNode << " while waiting for a positive quantity");
		distList.push_back(distToNode);
		indicesNodes.push_back(node);
	}

	// 2.3. Finishing the work when there is no contact
	if (!distList.size()) { // all boundary nodes are outside the bounds' overlap,
		if (!c->isReal()) return false;
		else {
			Ig2_LevelSet_LevelSet_ScGeom::geomPtrForLaterRemoval(state1, state2, c, scene);
			return true;
		}
	}
	Real maxOverlap;
	maxOverlap = *std::max_element(distList.begin(), distList.end());
	if (maxOverlap < 0 && !c->isReal() && !force) // inspired by Ig2_Sphere_Sphere_ScGeom:
		return false;

	// 2.4. Finishing the work when there is a contact
	int indexCont = std::min_element(distList.begin(), distList.end())
	        - distList.begin();                                       //this works: it seems min_element is returning here a random access iterator
	Vector3r           normal((lsC > boxC ? 1. : -1.) * axisContact); // normal from the box body to the level set body, ie from 1 to 2, as expected.
	Real               rad((lsNodes[indicesNodes[indexCont]] - lsCenter).norm());
	shared_ptr<ScGeom> geomPtr;
	bool               isNew = !c->geom;
	if (isNew) {
		geomPtr = shared_ptr<ScGeom>(new ScGeom());
		c->geom = geomPtr;
	} else
		geomPtr = YADE_PTR_CAST<ScGeom>(c->geom);
	geomPtr->doIg2Work(
	        lsNodes[indicesNodes[indexCont]] + maxOverlap / 2. * normal, // middle of overlapping volumes, as usual
	        maxOverlap,                                                  // does not work for very big/huge overlap
	        rad,
	        rad,
	        state1,
	        state2,
	        scene,
	        c,
	        normal,
	        shift2,
	        isNew,
	        false);
	return true;
}

bool Ig2_Wall_LevelSet_ScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
#ifdef USE_TIMING_DELTAS
	timingDeltas->start();
#endif
	shared_ptr<Wall>     wallSh = YADE_PTR_CAST<Wall>(shape1);
	shared_ptr<LevelSet> lsSh   = YADE_PTR_CAST<LevelSet>(shape2);
	Real                 lsPos(state2.pos[wallSh->axis] + shift2[wallSh->axis]), wallPos(state1.pos[wallSh->axis]);

	const int nNodes(lsSh->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");
	if (wallSh->sense != 0) LOG_ERROR("Ig2_Wall_LevelSet_* only supports Wall.sense = 0");

	Real distToNode, // one wall-level set distance value (< 0 when contact), for one node
	        nodePos, // current position along the Wall->axis of one given boundary node
	        maxOverlap(-std::numeric_limits<Real>::infinity());
	Vector3r currNode,   // current position of one given boundary node
	        contactNode(math::NaN,math::NaN,math::NaN); // the boundary node which is the most inside the wall
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("Until nodes loop");
#endif
	for (int node = 0; node < nNodes; node++) {
		currNode = ShopLS::rigidMapping(lsSh->surfNodes[node], Vector3r::Zero(), state2.pos + shift2, state2.ori);
		nodePos  = currNode[wallSh->axis];
		if (wallPos >= lsPos && wallPos <= nodePos) // first possibility for the wall to intersect the LevelSet body
			distToNode = wallPos - nodePos;
		else if (wallPos >= nodePos && wallPos <= lsPos) // second possibility for intersection
			distToNode = nodePos - wallPos;
		else
			continue; // go directly to next node
		if (distToNode > 0) LOG_ERROR("Unexpected case ! We obtained " << distToNode << " while waiting for a negative quantity");
		if (-distToNode > maxOverlap) {
			maxOverlap  = -distToNode;
			contactNode = currNode;
		}
	}
	// After that loop we may have (-)infinite maxOverlap and other garbage data in contactNode in case no node was contacting, which would be a problem in case of a neverErase
	if (!std::isfinite(maxOverlap)) {
		LOG_DEBUG("Interaction " << c->id1 << "-" << c->id2 << ": passing here since no node was found to be contacting");
		maxOverlap  = -1;                                      // any arbitrary negative value
		contactNode = state2.pos + shift2 + Vector3r(1, 1, 1); // some quite arbitrary value, that still aims to avoid a 0 "rad" below
	}
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("After completing nodes loop");
#endif
	if (maxOverlap < 0 && !c->isReal() && !force)
		return false; // we won't create the interaction in this case (but it is not our job here to delete it in case it already exists)
	Vector3r wallNormal(Vector3r::Zero()), normal(Vector3r::Zero());
	wallNormal[wallSh->axis] = 1;
	normal                   = (wallPos - lsPos > 0 ? -1 : 1) * wallNormal; // Points from wall to particle center
	Real rad((contactNode - state2.pos - shift2).norm());                   // Distance from surface to center of level-set body

	shared_ptr<ScGeom> geomPtr;
	bool               isNew = !c->geom;
	if (isNew) {
		geomPtr = shared_ptr<ScGeom>(new ScGeom());
		c->geom = geomPtr;
	} else
		geomPtr = YADE_PTR_CAST<ScGeom>(c->geom);
	geomPtr->doIg2Work(
	        contactNode + maxOverlap / 2. * normal, // middle of overlapping volumes, as usual, with garbage but finite data in case of tensile separation
	        maxOverlap,                             // does not work for very big/huge overlap
	        rad, // considering the 2* feature of radius* (see comments in ScGeom::doIg2work), this is what makes most sense ?
	        rad, // we keep in particular radius1/2 slightly greater than (contactPoint-center).norm. And we just use the same radii for the two particles, as in Ig2_Box_Sphere_ScGeom
	        state1,
	        state2,
	        scene,
	        c,
	        normal,
	        shift2,
	        isNew,
	        false);
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("Terminating ::go");
#endif
	return true;
}

bool Ig2_Wall_LevelSet_MultiScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	bool                       newIgeom(!bool(c->geom));
	shared_ptr<MultiScGeom>    geomMulti(new MultiScGeom);
	shared_ptr<MultiPhys> physMulti(new MultiPhys); // As a reminder, Ig2_*MultiScGeom have to touch the physics
	if (!newIgeom) {
		geomMulti = YADE_PTR_CAST<MultiScGeom>(c->geom);
		physMulti = YADE_PTR_CAST<MultiPhys>(c->phys);
	} // nothing to do in else
	shared_ptr<Wall>     wallSh = YADE_PTR_CAST<Wall>(shape1);
	shared_ptr<LevelSet> lsSh   = YADE_PTR_CAST<LevelSet>(shape2);
	Real                 lsPos(state2.pos[wallSh->axis] + shift2[wallSh->axis]), wallPos(state1.pos[wallSh->axis]);

	const int nNodes(lsSh->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");
	if (wallSh->sense != 0) LOG_ERROR("Ig2_Wall_LevelSet_* only supports Wall.sense = 0");

	// Before even looking at nodes below, we can compute the normal (which would be a waste if no nodes are touching, but let it be this way):
	Vector3r wallNormal(Vector3r::Zero()), normal(Vector3r::Zero());
	wallNormal[wallSh->axis] = 1;
	normal                   = (wallPos - lsPos > 0 ? -1 : 1) * wallNormal; // points from wall to particle center
	// We will now look at nodes:
	Real distToNode{math::NaN},   // one wall-level set distance value (< 0 when contact), for one node
	        nodePos{math::NaN};   // current position along the Wall->axis of one given boundary node
	Vector3r currNode(math::NaN,math::NaN,math::NaN); // current position of one given boundary node
	for (int node = 0; node < nNodes; node++) {
		currNode = ShopLS::rigidMapping(lsSh->surfNodes[node], Vector3r::Zero(), state2.pos + shift2, state2.ori);
		nodePos  = currNode[wallSh->axis];
		LOG_DEBUG(
		        "Looking at node " << node << " located (in current configuration) at " << currNode << " while unitialized distToNode is "
		                           << distToNode)
		if (wallPos >= lsPos
		    && wallPos
		            <= nodePos) // first possibility for the wall to intersect the LevelSet body in a LScenter - Wall - LSsurface alignment when walking along axis
			distToNode = wallPos - nodePos;
		else if (
		        wallPos >= nodePos
		        && wallPos <= lsPos) // second possibility for intersection, in a reverse (LSsurface - Wall - LScenter) alignment along axis
			distToNode = nodePos - wallPos;
		else { // that node does not currently involve contact with the wall
			// we still have to wonder whether that was the case before (to take it out of lists if yes)
			ShopLS::handleNonTouchingNodeForMulti(geomMulti, physMulti, node);
			continue; // and then move on to the next node
		}
		if (distToNode > 0) LOG_ERROR("Unexpected case ! We obtained " << distToNode << " while waiting for a negative quantity");
		if (force) LOG_WARN("force is true, let us think about that");
		// from now on, we know this node is contacting
		Real rad((currNode - state2.pos - shift2).norm());
		ShopLS::handleTouchingNodeForMulti(
		        geomMulti,
		        physMulti,
		        node,
		        currNode - distToNode / 2. * normal, // middle of overlapping volumes, as usual
		        -distToNode,
		        rad,
		        rad,
		        state1,
		        state2,
		        scene,
		        c,
		        normal,
		        shift2);
	}
	LOG_DEBUG(
	        "Nodes loop done, we finally got " << geomMulti->contacts.size() << " in geomMulti.contacts and " << physMulti->contacts.size()
	                                           << " in physMulti for that interaction which is " << (c->isReal() ? "real" : "not real"));
	if (geomMulti->contacts.size() == 0
	    && !c->isReal() // this typically happens just after Aabb overlap for the first time but there is no node-surface contact yet
	    && !force)
		return false; // we won't create the interaction in this case
	if (newIgeom) {
		c->geom = geomMulti;
		c->phys = physMulti;
	}
	LOG_DEBUG("Returning true");
	return true;
}

void Ig2_LevelSet_LevelSet_ScGeom::geomPtrForLaterRemoval(const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c, const Scene* scene)
{
	// to use when we can not really compute anything, e.g. bodies lsGrid do not overlap anymore, but still need to have some geom data (while returning true as per general InteractionLoop workflow because it is an existing interaction. Otherwise we would need to update InteractionLoop itself to avoid LOG_WARN messages). Data mostly include an infinite tensile stretch to insure subsequent interaction removal (by Law2)
	// this function is only applied onto a real interaction c, i.e. with an existing geom
	shared_ptr<ScGeom> geomPtr(YADE_PTR_CAST<ScGeom>(c->geom));
	LOG_DEBUG(
	        "YADE_PTR_DYN_CAST<ScGeom>(c->geom) = "
	        << YADE_PTR_DYN_CAST<ScGeom>(geomPtr)
	        << " (0 would be a problem maybe because of an unexpected MultiScGeom)"); // the dynamic cast will never be executed here unless LOG_DEBUG (level) is actually desired
	geomPtr->doIg2Work(
	        Vector3r::Zero() /* inconsequential bullsh..*/,
	        -std::numeric_limits<Real>::infinity() /* arbitrary big tensile value to trigger interaction removal by Law2*/,
	        1, /* inconsequential bullsh..*/
	        1, /* inconsequential bullsh..*/
	        rbp1,
	        rbp2,
	        scene,
	        c,
	        Vector3r::UnitX() /* inconsequential bullsh..*/,
	        Vector3r::Zero(), /* inconsequential bullsh..*/
	        false,            /* inconsequential bullsh..*/
	        false);
}

bool Ig2_LevelSet_LevelSet_LSnodeGeom::go(
        const shared_ptr<Shape>& shape1,
        const shared_ptr<Shape>& shape2,
        const State&             state1,
        const State&             state2,
        const Vector3r&          shift2,
        const bool&, // "force" (always false for the record), unused here
        const shared_ptr<Interaction>& c)
{
#ifdef USE_TIMING_DELTAS
	timingDeltas->start();
#endif
	// 1. We first use boundOverlap() very much just like Ig2_LevelSet_LevelSet_ScGeom::goSingleOrMulti
	auto returnOfBoundOverlap(boundOverlap(true, state1, state2, c, shift2));
	if (returnOfBoundOverlap.second.first)                                                                    // an empty bound overlap was detected
		return (returnOfBoundOverlap.second.second);                                                      // we shall stop Ig2*::go work
	Vector3r minBoOverlap(returnOfBoundOverlap.first.first), maxBoOverlap(returnOfBoundOverlap.first.second); // current configuration obviously
	// stop of ~ copied-pasted 1.

	// 2. We go now for the master-slave contact algorithm with surface ie boundary nodes:
	// 2.1 Preliminary declarations, starting some copy-paste from Ig2_LevelSet_LevelSet_ScGeom::goSingleOrMulti:

	// ****** following block initially very much copied-pasted ******** !
	shared_ptr<LevelSet> sh1 = YADE_PTR_CAST<LevelSet>(shape1);
	shared_ptr<LevelSet> sh2 = YADE_PTR_CAST<LevelSet>(shape2);
	bool                 id1isBigger(
                sh1->getVolume() > sh2->getVolume()); // identifying the smallest particle (where the master-slave contact algorithm will locate boundary nodes)

	shared_ptr<LevelSet> shS(id1isBigger ? sh2 : sh1);                          // the shape of the small body, useful to have it defined here
	shared_ptr<LevelSet> shB(id1isBigger ? sh1 : sh2);                          // the shape of the big body, useful to have it defined here
	Vector3r             ctrS(id1isBigger ? (state2.pos + shift2) : state1.pos) // center of the small body
	        ,
	        ctrB(id1isBigger ? state1.pos : (state2.pos + shift2)) // center of the big body
	        ;
	Quaternionr rotS(
	        id1isBigger ? state2.ori
	                    : state1.ori) // ori = rotation from reference configuration (local axes) to current one, useful to have it defined here
	        ,
	        rotB(id1isBigger ? state1.ori : state2.ori) // orientation of the big body
	        ;
	// ****** end of copy-paste ? ********

	// 2.2 We now handle surface nodes searching for the closest one with various workflows depending on surfNodeIdx available data
	std::tuple<Real, int, Vector3r> closestNodeData; // its (distance value; idx; position in global configuration)
	shared_ptr<LSnodeGeom>          geom(new LSnodeGeom);
	bool                            newIgeom(!bool(c->geom));
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("Until the history-dependent stuff with nodes loop");
#endif
	if (newIgeom) {
		LOG_INFO("A brand new LSnodeGeom at " << c->id1 << "-" << c->id2 << " contact");
		geom->usedWorkflow = 0;
		closestNodeData    = smallestDistanceFullLoop(
                        shS,
                        shB // both shapes
                        ,
                        ctrS // center of the small body
                        ,
                        rotS // orientation of the small body
                        ,
                        ctrB // center of the big body
                        ,
                        rotB // orientation of the big body
                        ,
                        minBoOverlap,
                        maxBoOverlap);
	} else {
		geom = YADE_PTR_CAST<LSnodeGeom>(c->geom); // we have an existing geom we can pick
		LOG_INFO("An existing LSnodeGeom..");
		if (geom->surfNodeIdx
		    < 0) { // but we have not yet identified a closest-node, maybe because they re still all outside the big body bound / lsGrid, let us then keep the whole loop
			geom->usedWorkflow = 1;
			LOG_INFO(".. with no surfNodeIdx defined at " << c->id1 << "-" << c->id2 << " contact");
			// exact same as above:
			closestNodeData = smallestDistanceFullLoop(
			        shS,
			        shB // both shapes
			        ,
			        ctrS // center of the small body
			        ,
			        rotS // orientation of the small body
			        ,
			        ctrB // center of the big body
			        ,
			        rotB // orientation of the big body
			        ,
			        minBoOverlap,
			        maxBoOverlap);
		} else { // using surfNodeIdx data
			geom->usedWorkflow = 2;
			LOG_INFO(".. with surfNodeIdx defined at " << c->id1 << "-" << c->id2 << " contact");
			// we will move accross neighbors as long as we find decreasing phi values, which is defined by this boolean:
			bool move(true), lastThoroughSearch(false) // that other boolean is for triggering when necessary a final full search
			        ;
			if (shS->n_neighborsNodes < 2) // all default -1, 0 and 1 would not be very good..
				LOG_ERROR(
				        "Impossible to use an optimized surfNodeIdx-based node loop with just "
				        << shS->n_neighborsNodes << " neighboring nodes (+ the guy) defined in LevelSet.n_neighborsNodes");
			while (move) {
				LOG_DEBUG("Searching for decreasing phi in neighbors of " << geom->surfNodeIdx);
				closestNodeData = smallestDistanceLoop(
				        shS->neighborsNodes[geom->surfNodeIdx],
				        shS,
				        shB // both shapes
				        ,
				        ctrS,
				        rotS // center and orientation of the small body
				        ,
				        ctrB,
				        rotB // center and orientation of the big body
				        ,
				        minBoOverlap,
				        maxBoOverlap);
				if (std::get<1>(closestNodeData) == geom->surfNodeIdx) {
					// we found back the same node: we are at phi-minimum ! Let us stop the search
					move = false;
					LOG_DEBUG("Closest distance still found at the same node, stoping to move");
				} else if (std::get<1>(closestNodeData) < 0) {
					// this is possible when previous contacting node and all its neighbors are now impossible to test in the distance field
					lastThoroughSearch = true;  // in this case we decide to redo (below) a full search to be sure
					move               = false; // we naturally need to stop working with neighbors
				} else {
					LOG_DEBUG("Closest distance now found at " << std::get<1>(closestNodeData));
					geom->surfNodeIdx = std::get<1>(closestNodeData);
				}
			}
			if (lastThoroughSearch) {
				closestNodeData = smallestDistanceFullLoop(
				        shS,
				        shB // both shapes
				        ,
				        ctrS // center of the small body
				        ,
				        rotS // orientation of the small body
				        ,
				        ctrB // center of the big body
				        ,
				        rotB // orientation of the big body
				        ,
				        minBoOverlap,
				        maxBoOverlap);
				geom->usedWorkflow = 3;
			}
		}
	}
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("After completing the history-dependent stuff");
#endif
	// 2.3 We finish the work

	// we now know about closestNodeData (and may have already filled geom->surfNodeIdx, for the record), with the possibility still of garbage data (e.g., nodeIdxOfInterest < 0) in case no node could actually be tested against the other distance field:
	// NB: after several tries, we will in this case indirectly trigger interaction removal by returning false (be it real or not)
	int nodeIdxOfInterest(std::get<1>(closestNodeData));
	if (nodeIdxOfInterest < 0) return false;
	Real     smallestPhi(std::get<0>(closestNodeData));
	Vector3r contactNode(std::get<2>(closestNodeData)); // in global configuration
	LOG_INFO("Got closest node as " << nodeIdxOfInterest << " with phi = " << smallestPhi);
	if (!std::isfinite(smallestPhi))
		LOG_FATAL(
		        "Unexpected case at interaction " << c->id1 << "-" << c->id2 << " with nodeIdx = " << nodeIdxOfInterest
		                                          << " but corresponding phi = " << smallestPhi);
	// NB: apart the above case and unlike, e.g., Ig2_Sphere_Sphere_ScGeom, we will here always create the geom (i.e., never return false) in order to use surfNodeIdx info even if id1, id2 are not yet in actual geometric contact
	Vector3r normal(
	        nodeIdxOfInterest >= 0 ? // that test for a special case is now no longer justified ?
	                rotS
	                        * shS->normal(
	                                shS->surfNodes
	                                        [nodeIdxOfInterest]) // shS->surfNodes[..] (and normal() to itself) refers to the initial shape, current normal is obtained with rotS *. It is for now the outward normal to the small particle
	                               : Vector3r(1, 0, 0)           // assigning a completely random (but finite) normal for those distant interactions
	);
	if (id1isBigger) normal *= -1;         // if necessary, we make the normal from 1 to 2, as expected
	geom->surfNodeIdx = nodeIdxOfInterest; // might be redundant with a line above for usedWorkflow = 2, let it be

	// Defining geometric quantities with special cases for the possibility of a distant interaction with nothing meaningfull computed (that no longer exists ? and the special cases could be removed)
	Vector3r middle12((state1.pos + state2.pos) / 2);
	Vector3r ctctPt(
	        std::isfinite(smallestPhi)
	                ? (contactNode + (id1isBigger ? -1 : 1) * smallestPhi / 2. * normal) // middle of overlapping volumes (or gap), as usual
	                : middle12 // in the case of no overlap we do not want to deal with smallestPhi in the quantity since it is infinite.. Ideal would be to have the middle point of the gap, but it is not obvious to know it in case of non-spherical particles.. We thus quite arbitrarily take the middle of the centers (even though it is not even sure it is in the gap)
	);
	Real rad1(
	        std::isfinite(smallestPhi)
	                ? (contactNode - state1.pos + smallestPhi * (id1isBigger ? normal : Vector3r::Zero()))
	                          .norm() // radius1 value: taking center to surface distance for id1, expressed from contactNode and possible offset (depending on which particle contactNode came from)
	                : sh1->maxRad // another finite arbitrary value
	        ),
	        rad2(std::isfinite(smallestPhi) ? (contactNode - state2.pos - shift2 + smallestPhi * (id1isBigger ? Vector3r::Zero() : normal))
	                                                  .norm() // radius2: same thing for id2. Reminder: contactNode does belong to id2 if id1 is bigger)
	                                        : sh2->maxRad     // another finite arbitrary value
	        );
	geom->doIg2Work(
	        ctctPt,
	        std::isfinite(smallestPhi) ? -smallestPhi // does not work for very big/huge overlap, eg when one sphere's center gets into another.
	                                   : -1,          // any negative value
	        rad1,
	        rad2,
	        state1,
	        state2,
	        scene,
	        c,
	        normal,
	        shift2,
	        newIgeom,
	        false);
	if (newIgeom) c->geom = geom;
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("Terminating ::go");
#endif
	return true;
}

bool Ig2_LevelSet_LevelSet_ScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
#ifdef USE_TIMING_DELTAS
	timingDeltas->start();
#endif
	bool ret(goSingleOrMulti(true, shape1, shape2, state1, state2, force, c, shift2));
#ifdef USE_TIMING_DELTAS
	timingDeltas->checkpoint("Completing ::go");
#endif
	return ret;
}

bool Ig2_LevelSet_LevelSet_MultiScGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	return goSingleOrMulti(false, shape1, shape2, state1, state2, force, c, shift2);
}

std::pair<std::pair<Vector3r, Vector3r>, std::pair<bool, bool>>
Ig2_LevelSet_LevelSet_ScGeom::boundOverlap(bool single, const State& state1, const State& state2, const shared_ptr<Interaction>& c, const Vector3r& shift2)
{
	// Shall determine the Aabb zone where bodies bounds overlap (in .first) , and also prepare for a premature return of the Ig2::go if there is none, when .second.first will be true, with .second.second the return value to be sent by Ig2::go
	// TODO: possible use of Eigen AlignedBox ?
	std::array<Real, 6>     overlap; // the xA,xB, yA,yB, zA,zB defining the Aabb where bounds overlap: for x in [xA,xB] ; y in [yA,yB] ; ...
	const shared_ptr<Bound> bound1 = Body::byId(c->id1, scene)->bound;
	const shared_ptr<Bound> bound2 = Body::byId(c->id2, scene)->bound;
	bool                    emptyBoOverlap(false) // whether there is no bounds overlap
	        ,
	        returnValue(true) // will be passed to Ig2*::go as the return value to send if(emptyBoOverlap)
	        ;
	for (int axis = 0; axis < 3; axis++) {
		overlap[2 * axis]     = math::max(bound1->min[axis], bound2->min[axis] + shift2[axis]);
		overlap[2 * axis + 1] = math::min(bound1->max[axis], bound2->max[axis] + shift2[axis]);
		if (overlap[2 * axis + 1]
		    < overlap[2 * axis]) {     // an empty bound overlap here (= is possible and OK: it happens when bodies bounds are just tangent)
			emptyBoOverlap = true; // we already know there is no bound overlap if there is no overlap on one axis
			if (!c->isReal()) returnValue = false;
			else {
				if (single) geomPtrForLaterRemoval(state1, state2, c, scene);
				else {
					shared_ptr<MultiScGeom> multiGeom(YADE_PTR_CAST<MultiScGeom>(c->geom));
					multiGeom->contacts
					        .clear(); // and Law2_MultiScGeom_MultiFrictPhys_CundallStrack::go will also return false and eventually trigger interaction removal
				}
				returnValue = true;
			}
			break; // no need to check other axes
		}
	}
	return std::pair<std::pair<Vector3r, Vector3r>, std::pair<bool, bool>>(
	        std::make_pair(
	                Vector3r(overlap[0], overlap[2], overlap[4]),
	                Vector3r(
	                        overlap[1],
	                        overlap[3],
	                        overlap[5])) // (some of) these could be garbage values if emptyBoOverlap but in this case we will never use them anyway
	        ,
	        std::make_pair(emptyBoOverlap, returnValue));
}

Real Ig2_LevelSet_LevelSet_ScGeom::normalMismatch(const shared_ptr<Interaction>& c)
{
	std::pair<Vector3r, Vector3r> normals; // .first, resp. .second, will be computed considering nodes are on 1, resp. on 2.
	const shared_ptr<Body>&       b1 = Body::byId(c->id1, scene), b2 = Body::byId(c->id2, scene);
	const shared_ptr<LevelSet>    sh1 = YADE_PTR_CAST<LevelSet>(b1->shape);                         // let s go for the static cast
	const shared_ptr<LevelSet>    sh2 = YADE_PTR_CAST<LevelSet>(b2->shape);                         // and leave the responsability to the user
	auto returnOfBoundOverlap(boundOverlap(true, *(b1->state), *(b2->state), c, Vector3r::Zero())); // TODO: shift2 here if PBC are to be supported
	if (returnOfBoundOverlap.second.first)                                                          // an empty bound overlap was detected
		LOG_ERROR("Given interaction shows no Bound overlap, how is this possible ?)");
	Vector3r minBoOverlap(returnOfBoundOverlap.first.first), maxBoOverlap(returnOfBoundOverlap.first.second); // current configuration obviously
	auto     closestNodeOn1Data
	        = smallestDistanceFullLoop(sh1, sh2, b1->state->pos, b1->state->ori, b2->state->pos, b2->state->ori, minBoOverlap, maxBoOverlap);
	auto closestNodeOn2Data
	        = smallestDistanceFullLoop(sh2, sh1, b2->state->pos, b2->state->ori, b1->state->pos, b1->state->ori, minBoOverlap, maxBoOverlap);
	normals.first  = b1->state->ori * sh1->normal(sh1->surfNodes[std::get<1>(closestNodeOn1Data)]);
	normals.second = b2->state->ori * sh2->normal(sh2->surfNodes[std::get<1>(closestNodeOn2Data)]);
	LOG_DEBUG("Normal computed on 1 is " << normals.first << " while on 2 it is " << -normals.second);
	const shared_ptr<ScGeom> geom = YADE_PTR_CAST<ScGeom>(c->geom);
	if ((normals.first.cross(geom->normal).norm() > 1.e-7) and (normals.second.cross(geom->normal).norm() > 1.e-7))
		LOG_ERROR("Could not reproduce cont.geom.normal, how is that possible ?");
	return math::asin(normals.first.cross(normals.second).norm()); // a - normals.second would be logical but useless because of .norm()
}

std::tuple<Real, int, Vector3r> Ig2_LevelSet_LevelSet_ScGeom::smallestDistanceFullLoop(
        const shared_ptr<LevelSet>& shS,
        const shared_ptr<LevelSet>& shB,
        Vector3r                    ctrS,
        Quaternionr                 rotS // describing "S" configuration
        ,
        Vector3r    ctrB,
        Quaternionr rotB // describing "B" configuration
        ,
        Vector3r minBoOverlap,
        Vector3r maxBoOverlap // the min and max of the bounds overlap
)
{
	vector<int> allIndices;
	allIndices.resize(shS->surfNodes.size());
	std::iota(allIndices.begin(), allIndices.end(), 0); // negative performance impact of that line ? or not ?
	return smallestDistanceLoop(allIndices, shS, shB, ctrS, rotS, ctrB, rotB, minBoOverlap, maxBoOverlap);
}

std::tuple<Real, int, Vector3r> Ig2_LevelSet_LevelSet_ScGeom::smallestDistanceLoop(
        vector<int>                 idxOfNodesToLoopOn,
        const shared_ptr<LevelSet>& shS,
        const shared_ptr<LevelSet>& shB,
        Vector3r                    ctrS,
        Quaternionr                 rotS // describing "S" configuration
        ,
        Vector3r    ctrB,
        Quaternionr rotB // describing "B" configuration
        ,
        Vector3r minBoOverlap,
        Vector3r maxBoOverlap // the min and max of the bounds overlap
)
{
	Vector3r minGridB(shB->lsGrid->min), maxGridB(shB->lsGrid->max()); // the min and max of the lsGrid of "B"
	Real     smallestPhi(std::numeric_limits<Real>::infinity())        // will be the first return value
	        ,
	        distToNode // current distance value that will often change below
	        ;
	Vector3r nodeOfS // node at hand, in global configuration
	        ,
	        nodeOfSinB // mapped into "B" configuration
	        ,
	        smallestDistanceNode // nodeOfS value showing smallestPhi
	        ;
	int smallestDistanceNodeIdx(-1) // second return item
	        ;
	LOG_DEBUG("Will loop over " << idxOfNodesToLoopOn.size() << " surface nodes");
	for (const int& nodeIdx : idxOfNodesToLoopOn) {
		if (nodeIdx < 0 || nodeIdx >= int(shS->surfNodes.size()))
			LOG_ERROR("Major problem with nodeIdx " << nodeIdx << " for looping over " << shS->surfNodes.size() << " elements.");
		// 1. We first look at bounds overlap and grid considerations:
		nodeOfS = ShopLS::rigidMapping(
		        shS->surfNodes[nodeIdx],
		        Vector3r::Zero(),
		        ctrS,
		        rotS); // current position of this boundary node. Note that with a correct ls body creation and use, does not matter whether we take LevelSet->getCenter() or 0
		if (!Shop::isInBB(nodeOfS, minBoOverlap, maxBoOverlap)) continue;
		nodeOfSinB = ShopLS::rigidMapping(nodeOfS, ctrB, Vector3r::Zero(), rotB.conjugate()); // mapping this node into the big Body local frame
		if (!Shop::isInBB(
		            nodeOfSinB,
		            minGridB,
		            maxGridB)) // possible when bodies (and their shape.corners) rotate, leading their bounds to possibly "inflate" (think of a sphere)
			continue;
		// 2. Actual distance considerations:
		distToNode = shB->distance(nodeOfSinB);
		LOG_DEBUG("For node " << nodeIdx << " just computed distance = " << distToNode);
		if (distToNode < smallestPhi) {
			// storing just-computed stuff, while the rest (like computing normal) will be done once for all outside that method
			smallestDistanceNodeIdx = nodeIdx;
			smallestPhi             = distToNode;
			smallestDistanceNode    = nodeOfS;
		}
	} // at this stage it is possible we return an infinite phi together with a negative node idx (if all were outside of bounding box/grid considerations) and a garbage nodeOfS
	// this will be taken care of in calling functions
	return std::make_tuple(smallestPhi, smallestDistanceNodeIdx, smallestDistanceNode);
}

bool Ig2_LevelSet_LevelSet_ScGeom::goSingleOrMulti(
        bool                           single,
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const bool&                    force,
        const shared_ptr<Interaction>& c,
        const Vector3r&                shift2)
{
	// 1. We first use boundOverlap() that looked at the Aabb zone where bodies bounds overlap
	auto returnOfBoundOverlap(boundOverlap(single, state1, state2, c, shift2));
	if (returnOfBoundOverlap.second.first)                                                                    // an empty bound overlap was detected
		return (returnOfBoundOverlap.second.second);                                                      // we shall stop Ig2*::go work
	Vector3r minBoOverlap(returnOfBoundOverlap.first.first), maxBoOverlap(returnOfBoundOverlap.first.second); // current configuration obviously

	// 2. We go now for the master-slave contact algorithm with surface ie boundary nodes:
	// 2.1 Preliminary declarations:
	bool                       newIgeom(!bool(c->geom));
	shared_ptr<ScGeom>         geomSingle(new ScGeom);     // useful only if single but let s also make it here in full scope
	shared_ptr<MultiScGeom>    geomMulti(new MultiScGeom); // useful only if !single but has to be in full scope
	shared_ptr<MultiPhys> physMulti(
	        new MultiPhys); // useful only if !single but has to be in full scope. As a reminder, Ig2_*MultiScGeom have to touch the physics
	if (!newIgeom) {
		if (single) geomSingle = YADE_PTR_CAST<ScGeom>(c->geom);
		else {
			geomMulti = YADE_PTR_CAST<MultiScGeom>(c->geom);
			physMulti = YADE_PTR_CAST<MultiPhys>(c->phys);
		}
	} // nothing to do in else
	shared_ptr<LevelSet> sh1 = YADE_PTR_CAST<LevelSet>(shape1);
	shared_ptr<LevelSet> sh2 = YADE_PTR_CAST<LevelSet>(shape2);
	bool                 id1isBigger(
                sh1->getVolume() > sh2->getVolume()); // identifying the smallest particle (where the master-slave contact algorithm will locate boundary nodes)
	shared_ptr<LevelSet> shB(id1isBigger ? sh1 : sh2); // the "big" shape
	shared_ptr<LevelSet> shS(id1isBigger ? sh2 : sh1); // the "small" shape
	// centr*end, and rot* below define the bodies' changes in configuration since the beginning of the simulation:
	Vector3r    centrSend(id1isBigger ? (state2.pos + shift2) : state1.pos), centrBend(id1isBigger ? state1.pos : (state2.pos + shift2));
	Quaternionr rotS(id1isBigger ? state2.ori : state1.ori), // ori = rotation from reference configuration (local axes) to current one
	        rotB(id1isBigger ? state1.ori : state2.ori)      // will be .conjugate() below to be from the current configuration to the reference one
	        ;
	const int nNodes(shS->surfNodes.size());
	if (!nNodes) LOG_ERROR("We will loop over no boundary nodes for contact detection. Will probably crash");
	Real     distToNode = 0;                                     // one distance value, for one node
	Real     smallestPhi(std::numeric_limits<Real>::infinity()); // useful only if single but has to be in full scope
	Vector3r minLSgrid   = shB->lsGrid->min;
	Vector3r maxLSgrid   = shB->lsGrid->max();
	Vector3r nodeOfS     = Vector3r::Zero(); // some node of smaller Body, in current configuration
	Vector3r nodeOfSinB  = Vector3r::Zero(); // mapped into initial configuration of larger Body
	Vector3r normal      = Vector3r::Zero();
	Vector3r contactNode = Vector3r::Zero();
	int      surfNodeIdx(-1);

	// 2.2 Actual loop over surface nodes:
	LOG_DEBUG("Will loop over " << nNodes << " surface nodes")
	for (int node = 0; node < nNodes; node++) {
		nodeOfS = ShopLS::rigidMapping(
		        shS->surfNodes[node],
		        Vector3r::Zero(),
		        centrSend,
		        rotS); // current position of this boundary node. Note that with a correct ls body creation and use, does not matter whether we take LevelSet->getCenter() or 0
		if (!Shop::isInBB(nodeOfS, minBoOverlap, maxBoOverlap)) {
			if (!single) ShopLS::handleNonTouchingNodeForMulti(geomMulti, physMulti, node);
			continue;
		}
		nodeOfSinB = ShopLS::rigidMapping(nodeOfS, centrBend, Vector3r::Zero(), rotB.conjugate()); // mapping this node into the big Body local frame
		if (!Shop::isInBB(
		            nodeOfSinB,
		            minLSgrid,
		            maxLSgrid)) { // possible when bodies (and their shape.corners) rotate, leading their bounds to possibly "inflate" (think of a sphere)
			if (!single) ShopLS::handleNonTouchingNodeForMulti(geomMulti, physMulti, node);
			continue;
		}
		distToNode = shB->distance(nodeOfSinB);
		if (single) {
			if (distToNode < 0 and distToNode < smallestPhi) { // when single, only greatest penetration depth is relevant
				surfNodeIdx = node; // storing the idx of that node. The rest (like computing normal) will be done once for all outside the loop
				smallestPhi = distToNode; // also storing just corresponding overlap
			}
		} else { // !single
			if (distToNode < 0) {
				normal = rotS * shS->normal(shS->surfNodes[node]); // getting the normal as above (single)
				if (id1isBigger) normal *= -1;                     // as above (single)
				contactNode = nodeOfS;                             // as above (single)
				ShopLS::handleTouchingNodeForMulti(
				        geomMulti,
				        physMulti,
				        node,
				        contactNode
				                + (id1isBigger ? -1 : 1) * distToNode / 2.
				                        * normal // middle of overlapping volumes, as usual TODO: is this correct or do we apply twice the id1isBigger -1
				        ,
				        -distToNode // does not work for very big/huge overlap, eg when one sphere's center gets into another
				        ,
				        (contactNode - state1.pos).norm() //TODO PBC
				        ,
				        (contactNode - state2.pos).norm() //TODO PBC
				        ,
				        state1,
				        state2,
				        scene,
				        c,
				        normal,
				        shift2);
			} else { // distToNode >= 0 and !single: we do not keep this as a contact point (even if distToNode = 0)
				ShopLS::handleNonTouchingNodeForMulti(geomMulti, physMulti, node);
				// and we automatically move on to the next node since we are at the end of the loop commands
			} // nothing to do in else (when the node was also not contacting before)
		}
	}
	if (single && surfNodeIdx >= 0) /* 2nd condition checks whether we did detect contact */ {
		normal = rotS
		        * shS->normal(
		                shS->surfNodes
		                        [surfNodeIdx]); // shS->surfNodes[..] (and normal() to itself) refers to the initial shape, current normal is obtained with rotS *. It is for now the outward normal to the small particle
		if (id1isBigger) normal *= -1;          // if necessary, we make the normal from 1 to 2, as expected
		contactNode = ShopLS::rigidMapping(
		        shS->surfNodes[surfNodeIdx],
		        Vector3r::Zero(),
		        centrSend,
		        rotS); // current position of this boundary node, this was already computed during the loop ?
	} else
		LOG_DEBUG(
		        "Having here " << geomMulti->contacts.size() << " guys in ig->contacts, vs " << physMulti->contacts.size() << " in ip->contacts\n And"
		                       << geomMulti->nodesIds.size() << " guys in ig->nodesIds, vs " << physMulti->nodesIds.size() << " in ip->nodesIds");

	// 2.3 Finishing the work:
	if (!c->isReal() && !force && ((single && smallestPhi > 0) or (!single && geomMulti->contacts.size() == 0)))
		return false; // we won't create the geom in this case, but it is not our job here to delete the interaction. We just return false on non-real ones, unless force (which will never happen since it is hard-coded to false in InteractionLoop). See e.g. Ig2_Sphere_Sphere_ScGeom
	if (single) {
		geomSingle->doIg2Work(
		        contactNode + (id1isBigger ? -1 : 1) * smallestPhi / 2. * normal, // middle of overlapping volumes, as usual
		        -smallestPhi, // does not work for very big/huge overlap, eg when one sphere's center gets into another.
		        (contactNode - state1.pos + smallestPhi * (id1isBigger ? normal : Vector3r::Zero()))
		                .norm(), // radius1 value: taking center to surface distance for id1, expressed from contactNode and possible offset (depending on which particle contactNode came from)
		        (contactNode - state2.pos - shift2 + smallestPhi * (id1isBigger ? Vector3r::Zero() : normal))
		                .norm(), // radius2: same thing for id2. Reminder: contactNode does belong to id2 if id1 is bigger
		        state1,
		        state2,
		        scene,
		        c,
		        normal,
		        shift2,
		        newIgeom,
		        false); // NB: the avoidGranularRatcheting=1 expression of relative velocity at contact does not make sense with radius1 and radius2 in all cases of LevelSet-sthg interaction (because the shapes are non-spherical): our present radius1+radius2 may have nothing in common with branch vector
	}
	if (newIgeom) {
		if (single) c->geom = geomSingle;
		else {
			c->geom = geomMulti; // ternary expression with above "= geomSingle" not possible because geomSingle and geomMulti are of different types
			c->phys = physMulti;
		}
	}
	if (!single) LOG_DEBUG("Will return true");
	return true;
}
} // namespace yade
#endif //YADE_LS_DEM
