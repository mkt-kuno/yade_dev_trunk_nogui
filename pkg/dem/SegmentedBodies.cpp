// 2024 © Karol Brzeziński <karol.brze@gmail.com>
// 2024 © Vasileios Angelidakis <angelidakis@qub.ac.uk>
#include <lib/high-precision/Constants.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/SegmentedBodies.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((CohFrictMatSeg)(Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys)(SegmentedStateUpdater)(SegmentedState)(SegmentedMatSprinkler)); //
CREATE_LOGGER(Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys);

/*********************************************************************************
* S T A T E 
*
*********************************************************************************/

SegmentedState::Segment SegmentedState::getSegmentFromContactPoint(Vector3r midPoint_)
{
        if (needsInit) throw runtime_error("Cannot get segment from contact point for body with uninitialized SegmentedState.");
	// function finds a proper
	int     segmentPosition;
	Segment segment;
	// Quaternion
	Quaternionr q            = ori;
	Vector3r    branchVector = midPoint_ - pos;
	if (thetaResolution > 0 && phiResolution > 0) { // sphere
		Real dH     = 2.0 / thetaResolution;
		Real dPhi = 2.0 * Mathr::PI / phiResolution;
		// branch vector in body coordinate system
		Vector3r v = q.inverse() * branchVector; // Maybe .toRotationMatrix() should be put after inverse(), we will see during validation
		// z coordinate of normalized branch vector
		Real h = v[2] / (v.norm()) + 1; // add 1 so it ranges from 0 to 2 instead of -1 to 1
		// computing angle
		Real phi = atan2(v[1], v[0]);
		if (phi < 0) {
			phi += 2 * Mathr::PI; // result was -pi to pi, I want it from zero to 2*pi
		}
		// determine to which segment it belongs
		int phiNo     = static_cast<int>(phi / dPhi);
		int thetaNo       = static_cast<int>(h / dH);
		segmentPosition = thetaNo * phiResolution + phiNo;
		//std::cout << "Q is: " << q<< " thetaResolution is: " << thetaResolution  << " h is: " << h << " dH is: " << dH << " thetaNo is: " << thetaNo   << std::endl;
	} else { // if Facet shape
		Vector3r fn = q * facetRefNormal;
		if (fn.dot(branchVector) < 0) { //
			segmentPosition = 0;
		} else {
			segmentPosition = 1;
		}
	}
	// get the area and thickness
	segment.segmentPosition = segmentPosition;
	segment.thickness       = coatingThickness[segmentPosition];
	segment.volume          = coatingVolume[segmentPosition];
	segment.area            = segmentArea;

	return segment;
}

void SegmentedState::initializeSegmentContacts()
{
        // note that this function do not have access to the shape of body (sphere radius)
        // here just angles are computed and normalized distances, the are multiplied by radius in SegmentedStateUpdater::initialize()
	if (thetaResolution > 0 && phiResolution > 0) { // only if sphere
		Real dH     = 2.0 / thetaResolution;
		Real dPhi = 2.0 * Mathr::PI / phiResolution;
		for (int thetaNo = 0; thetaNo < thetaResolution; thetaNo++){
		        // vertical components
		        Real vertEdgeAngle1 =  acos(1. - thetaNo * dH);// vertical edge means edge between segments with the same thetaNo
		        Real vertEdgeAngle2 =  acos(1. - (thetaNo + 1) * dH);
		        Real vertEdgeAngleCentral =  acos(1. - (thetaNo + 0.5) * dH);
		        Real vertEdgeAngleCentralNext =  acos(1. - (thetaNo + 1.5) * dH);// center of next segment
		        Real normalizedVertEdgeLen = (vertEdgeAngle2-vertEdgeAngle1); // should be enough since radius here is assumed 1
		        Real normalizedVertSegmentDist = (vertEdgeAngleCentralNext-vertEdgeAngleCentral);
		        // horizontal components also depend only on theta (loop over phi is to assign proper segment numbers)
		        // one just need to choose proper slice/cross-section of sphere and compute section of a circle
		        Real sectionR2 = sin(vertEdgeAngle2); // will be zero for thetaNo = thetaResolution,
		        Real sectionCentralR = sin(vertEdgeAngleCentral);
		        Real normalizedHorizontalEdgeLen2 = dPhi * sectionR2;
		        Real normalizedHorizontalSegmentDist = dPhi * sectionCentralR;
		        for (int phiNo = 0; phiNo < phiResolution; phiNo++){
		                // add "contacts" between segments. 
		                // To avoid duplicates, for segment denoted by (thetaNo, phiNo) add only contacts with segment denoted by
		                // (thetaNo + 1, phiNo), and (thetaNo, phiNo + 1).
		                // The only exception is phiNo = 0, because we need to make full circle, so we add also contact with phiNo = phiResolution.
                                // vertical contact
                                if (thetaNo < thetaResolution - 1){
                                        SegmentContact verticalContact;
                                        verticalContact.segmentPosition1 = thetaNo * phiResolution + phiNo;
                                        verticalContact.segmentPosition2 = (thetaNo + 1) * phiResolution + phiNo;
                                        verticalContact.distance = normalizedVertSegmentDist; // distance measured vertically but edge is horizontal
                                        verticalContact.edgeLen = normalizedHorizontalEdgeLen2;  
                                        // push back created contact
                                        segmentContacts.push_back(verticalContact);                                
                                }
                                // horizontal contact 
                                SegmentContact horizontalContact;
                                int nextPhiNo = phiNo + 1;
                                if (phiNo == phiResolution - 1) nextPhiNo = 0;
                                horizontalContact.segmentPosition1 = thetaNo * phiResolution + phiNo;
                                horizontalContact.segmentPosition2 = thetaNo * phiResolution + nextPhiNo;                                
                                horizontalContact.distance = normalizedHorizontalSegmentDist; // distance measured horizontally but edge is vertical
                                horizontalContact.edgeLen = normalizedVertEdgeLen; 
                                // push back created contact
                                segmentContacts.push_back(horizontalContact);                             
		        }		
		}
	}
        return;
}

/*********************************************************************************
*
* S T A T E   U P D A T E R 
*
*********************************************************************************/

void SegmentedStateUpdater::initialize(Body::id_t bId, Real thickness_)
{
	int  coatingThicknessVectorLen = 0;
	Real segmentArea               = 0;
	// body
	const auto b = Body::byId(bId, scene);
	if (!b) throw runtime_error("There is no body with given Id in the simulation.");
	// material
	shared_ptr<CohFrictMatSeg> mat = YADE_PTR_DYN_CAST<CohFrictMatSeg>(b->material);
	// try to cast shape into Sphere of Facet to check the shape of the body
	shared_ptr<Facet>  f = YADE_PTR_DYN_CAST<Facet>(b->shape);
	shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);

	SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get());
	if (!segState) throw runtime_error("Failed to obtain SegmentedState for body with provided id.");

	if (!s && !f) {
		throw runtime_error("In the current version of the code SegmentedState cannot be initialized for shapes different than Sphere or Facet.");
	} else {
		if (!s) {
			coatingThicknessVectorLen = 2;
			segmentArea               = f->area * 1e6; // * 1e6 converts to squared millimiters
			segState->facetRefNormal  = f->normal;
		} else {
			coatingThicknessVectorLen = thetaResolution * phiResolution;
			segmentArea               = 6 * Mathr::PI * pow(s->radius, 2) / coatingThicknessVectorLen * 1e6; // * 1e6 converts to squared millimiters
			segState->thetaResolution   = thetaResolution;
			segState->phiResolution = phiResolution;
			// initialize segment contacts here
			segState->initializeSegmentContacts();		    
			// Geometry of 'contacts' created by the above function was normalized by radius.
			// Hence, we will need multiply the distances by radius.
			// TODO Consider, whether the distances should be converted to mm, to be consistent with units assumed for thickness. For now - yes.
			for (int i = 0; i < int(segState->segmentContacts.size()); i++){
			        segState->segmentContacts[i].distance *= s->radius * 1000;
                                segState->segmentContacts[i].edgeLen *= s->radius * 1000;
			}

		}
	};


	for (long counter = 0; counter < coatingThicknessVectorLen; counter++) {
		segState->coatingThickness.push_back(thickness_); 
		segState->coatingVolume.push_back(thickness_ * segmentArea); 
	};

	segState->needsInit   = false;
	segState->segmentArea = segmentArea;

	return;
}


void SegmentedStateUpdater::setThicknessToSegmentRange(Body::id_t bId, int thetaMin_, int thetaMax_, int phiMin_, int phiMax_, Real thickness_)
{
	// body
	const auto b = Body::byId(bId, scene);
	if (!b) throw runtime_error("There is no body with given Id in the simulation.");
	// material
	shared_ptr<CohFrictMatSeg> mat = YADE_PTR_DYN_CAST<CohFrictMatSeg>(b->material);
	// checking if provided ranges are valid based on the resolution (from material)
	if (phiMin_ < 0) throw runtime_error("Theta coordinate of segment cannot be negative."); // maybe instead of this check we could use just unsigned type?
	if (thetaMin_ < 0) throw runtime_error("Phi coordinate of segment cannot be negative.");
	if (phiResolution - 1 < phiMax_)
		throw runtime_error("Provided phiMax coordinate range is bigger than maximum coordinate specified by phiResolution.");
	if (thetaResolution - 1 < thetaMax_) throw runtime_error("Provided thetaMax coordinate is bigger than maximum coordinate specified by thetaResolution.");

	SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get());
	if (!segState) throw runtime_error("Failed to obtain SegmentedState for body with provided id.");

	if (segState->needsInit)
		SegmentedStateUpdater::initialize(
		        bId, 0); // If thickness is to be set to one segment of uninitialized body, first initialize thickness with zeros
        Real area = segState->segmentArea;// the same for all the segments
	for (int theta = thetaMin_; theta < thetaMax_ + 1; theta++) {
		for (int phi = phiMin_; phi < phiMax_ + 1; phi++) {
			int segmentPosition                         = theta * phiResolution + phi;
			segState->coatingThickness[segmentPosition] = thickness_;
			segState->coatingVolume[segmentPosition] = thickness_ * area;
		}
	}

	return;
}

void SegmentedStateUpdater::setSegmentThickness(Body::id_t bId, int segmentPos_, Real thickness_)
{
	SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get());
	if (segState->needsInit) SegmentedStateUpdater::initialize(bId, 0);
	int maxSegmentPos = segState->coatingThickness.size() - 1;
	if ((maxSegmentPos < segmentPos_) || (segmentPos_ < 0)) {
		throw runtime_error("Segment position specified incorrectly.");
	} else {
		segState->coatingThickness[segmentPos_] = thickness_;
		segState->coatingVolume[segmentPos_] = thickness_ * segState->segmentArea;
	}

	return;
}


void SegmentedStateUpdater::setThicknessToSphereSeg(Body::id_t bId, int theta_, int phi_, Real thickness_)
{
	// body
	const auto b = Body::byId(bId, scene);
	if (!b) throw runtime_error("There is no body with given Id in the simulation.");
	// try to cast shape into Sphere
	shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
	if (!s) throw runtime_error("This function can only be used with Spheres.");

	setThicknessToSegmentRange(bId, theta_, theta_, phi_, phi_, thickness_);

	return;
}

void SegmentedStateUpdater::setThicknessToFacetSide(Body::id_t bId, int side_, Real thickness_)
{
	if (side_ != 0 && side_ != 1) runtime_error("Facet side should be specified with integer 0 or 1.");
	// body
	const auto b = Body::byId(bId, scene);
	if (!b) throw runtime_error("There is no body with given Id in the simulation.");
	// try to cast shape into Facet
	shared_ptr<Facet> f = YADE_PTR_DYN_CAST<Facet>(b->shape);
	if (!f) throw runtime_error("This function can only be used with Facets.");

	SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get());
	if (!segState) throw runtime_error("Failed to obtain SegmentedState for body with provided id.");
	if (segState->needsInit)
		SegmentedStateUpdater::initialize(
		        bId, 0); // If thickness is to be set to one segment of uninitialized body, first initialize thickness with zeros

	segState->coatingThickness[side_] = thickness_;
	segState->coatingVolume[side_] = thickness_ * segState->segmentArea;

	return;
}

void SegmentedStateUpdater::setThicknessToSpheres(vector<Body::id_t> bIds_, int thetaMin_, int thetaMax_, int phiMin_, int phiMax_, Real thickness_)
{
	long size = bIds_.size();
	for (long counter = 0; counter < size; counter++) {
		Body::id_t bId = bIds_[counter];
		// body
		const auto b = Body::byId(bId, scene);
		if (!b) throw runtime_error("There is no body with given Id in the simulation.");
		// try to cast shape into Sphere
		shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
		if (!s) throw runtime_error("This function can only be used with Spheres.");

		setThicknessToSegmentRange(bId, thetaMin_, thetaMax_, phiMin_, phiMax_, thickness_);
	};

	return;
}

// Computes volume of the material flowing from body 1 to body 2. If volume is negative it means that material flows from body 2 to body 1.
// The exchangeRate_ is to be computed in other parts of the code based on the geometry and coefficients.
Real SegmentedStateUpdater::computeVolumeExchange(Real thickness1_, Real thickness2_, Real volume1_, Real volume2_, Real exchangeRate_)
{
        Real volumeChange = (thickness1_ - thickness2_) * exchangeRate_ * (scene->time - prevTime );
        // If there is enough volume to exchange. If not - limit the volume
        if (volumeChange > 0){ // flow from body 1 to body 2
                volumeChange = min(volumeChange, volume1_);
        } else {// flow from body 2 to body 1
                volumeChange = min(volumeChange, volume2_);        
        }
        return volumeChange;
}

void SegmentedStateUpdater::action()
{
        if (nDone == 0) prevTime = scene->time; 
	//Hence, we could not be sure that the same data is not accessed by different threads.
	//First loop is over bodies (internal wettability) and can be parallelized. The second one is over interactions.
	// LOOP OVER BODIES
	if (activateWettability){
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif                    
                for (unsigned int id = 0; id < scene->bodies->size(); id++) {
                        const auto b = Body::byId(id, scene);
                        if (!b) continue;
                        SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(b->state.get());
                        if (!segState) continue;
                        if (!(segState->thetaResolution > 0 and segState->phiResolution > 0)) continue; // if not sphere
                        // First, exchange volumes between segments
                        int numContacts = segState->segmentContacts.size();
	                for (int j = 0; j < numContacts; j++){
                                SegmentedState::SegmentContact contact = segState->segmentContacts[j];
                                int pos1 = contact.segmentPosition1;
                                int pos2 = contact.segmentPosition2;
                                Real thickness1 = segState->coatingThickness[pos1];
                                Real thickness2 = segState->coatingThickness[pos2];
                                Real volume1 = segState->coatingVolume[pos1];
                                Real volume2 = segState->coatingVolume[pos2];

                                Real exchangeRate = intraParticleExchangeCoeff * contact.edgeLen / contact.distance;
                                Real volumeChange = computeVolumeExchange(thickness1, thickness2, volume1, volume2, exchangeRate);
                                
                                segState->coatingVolume[pos1] -= volumeChange;
                                segState->coatingVolume[pos2] += volumeChange;                
                        }    
                        // Next, apply new thickness based on the actual volumes
                        int numSegments = segState->coatingThickness.size();
                        for (int j = 0; j < numSegments; j++){
                            segState->coatingThickness[j] = segState->coatingVolume[j] / segState->segmentArea;
                        }   
	        }
	}
	// LOOP OVER INTERACTIONS
	// For now, we think it should not be parallelized, because one big facet (one segment) can have multiple interactions.
	for (auto it = scene->interactions->begin(); it != scene->interactions->end(); it++) {
		const shared_ptr<Interaction>& interaction = *it;
		ScGeom6D*                      geom        = YADE_CAST<ScGeom6D*>(interaction->geom.get());
		if (geom) {
			/// We know that it is redundant to the part in Ip2 (see comments in Ip2)
			const auto body1 = Body::byId(interaction->id1, scene);
			const auto body2 = Body::byId(interaction->id2, scene);

			SegmentedState* segState1 = YADE_DYN_CAST<SegmentedState*>(body1->state.get());
			SegmentedState* segState2 = YADE_DYN_CAST<SegmentedState*>(body2->state.get());
			if (segState1 && segState2) {
			        if (segState1->allowMatExchange && segState2->allowMatExchange && !segState1->needsInit && !segState2->needsInit) {
				        Vector3r midPoint = (segState1->pos + segState2->pos) / 2;// use mid point between two bodies instead of geom->contactPoint, the latter is usless because it lies on the facet, and wee need to determine on which side the contacting body is.

				        SegmentedState::Segment segment1;
				        SegmentedState::Segment segment2;

				        segment1 = segState1->getSegmentFromContactPoint(midPoint);
				        segment2 = segState2->getSegmentFromContactPoint(midPoint);

				        // compute new thickness and assign it to segments
				        Real volumeChange = computeVolumeExchange(segment1.thickness, segment2.thickness, segment1.volume, segment2.volume, interParticleExchangeCoeff);

				        segState1->coatingVolume[segment1.segmentPosition] = segment1.volume - volumeChange;
				        segState2->coatingVolume[segment2.segmentPosition] = segment2.volume + volumeChange;
				        
				        segState1->coatingThickness[segment1.segmentPosition] = segState1->coatingVolume[segment1.segmentPosition] / segState1->segmentArea;
				        segState2->coatingThickness[segment2.segmentPosition] = segState2->coatingVolume[segment2.segmentPosition] / segState2->segmentArea;
				        
				        // reinitialize cohesion
				        if (reInitCohesion == 1 || reInitCohesion == 2){
				                CohFrictPhys* contactPhysics = YADE_DYN_CAST<CohFrictPhys*>(interaction->phys.get());
				                if (reInitCohesion == 1 && contactPhysics->cohesionBroken) contactPhysics->initCohesion = true;
				                if (reInitCohesion == 2) contactPhysics->initCohesion = true;
				        }
        
				}
			}
		}
	}
	prevTime = scene->time;
}


/*********************************************************************************
*
* S E G M E N T E D   M A T   S P R I N K L E R 
*
*********************************************************************************/

void SegmentedMatSprinkler::init()
{
	needsInit = false;
	versors.clear();	
    // Let us imagine a matrix of pixels perpendicular to versor refDir.
    // The size of this matrix is defined by alpha and beta, while resulution by alphaResolution and betaResolution
    // We will compute versors poining towards the centers of the "pixels".
    Real alphaRad = alpha * Mathr::PI / 180;
    Real betaRad = beta * Mathr::PI / 180;
    Real xLength = 2 * sin(alphaRad/2);
    Real yLength = 2 * sin(betaRad/2);
    Real dx = xLength / alphaResolution;
    Real dy = yLength / betaResolution;
    Real xCurrent = -xLength/2 + dx/2; // adding dx/2 because versor points toward the center of the 'matrix field'
        
	for (int i = 0; i < alphaResolution; i++) {
	        Real yCurrent = -yLength/2 + dy/2;
                for (int j = 0; j < betaResolution; j++){

                        Vector3r vectXY = Vector3r(xCurrent,yCurrent,0);
                        Vector3r versor = refDir + vectXY;
                        versor = versor.normalized();
                        versors.push_back(versor);
                        yCurrent += dy;
                }
                xCurrent += dx;
	}
	
	// Code for storing the corner vectors of the sprinkled field.
	// Maybe we could select form the versors computed above, but I think it will be easier and cleaner to compute it separetely.
	// Also the vectors below are slightly different use the actual corners (not center of corner fields). Although, I am not sure if it is strictly necessary.
	// Finally, this is computed only during initialization of the engine. Hence, the cost is minimal.
	
	cornerVersors.clear();
	for (int j = -1; j < 2; j +=2) { // We want to set y -1 and then x -1 and 1, next set y 1, and x from 1 and -1, because we need to have consecutive corners.
	    for (int i = -1; i < 2; i +=2){
	        Vector3r cornerVersor = refDir + Vector3r(-1 * i * j * xLength / 2, j * yLength, 0);// If seems overcomplicated, see the comment above.
	        cornerVersor = cornerVersor.normalized();
	        cornerVersors.push_back(cornerVersor);
	    }	
	}
	
	
}

void SegmentedMatSprinkler::setVersors(Quaternionr ori_)
{
    init();
	int versorsSize = versors.size();
	for (int i = 0; i < versorsSize; i++) {
                versors[i] = ori_ * versors[i]; // This is rotation. However, versors are initialized earlier, so it results in setting the orientation with reference to refOri.
	}
	// also rotate corner versors
	for (int i = 0; i < 4; i++) {
                cornerVersors[i] = ori_ * cornerVersors[i];
	}
}



/* How parallelization is performed here:

SegmentedMatSprinkler::oneLineAction should only search within a set preselected spheres 

The spheres should be preselected based on the signed distance of the sphere to limiting plane.
The limiting planes are defined by sprinkler origin and cornerVersors. For example, cornerVersor[0] x cornerVersor[1] gives the normal vector pointed outside my volume.
Hence, for each sphere, I assume it is inside, but then I will iterate over four planes. If the condition is not met for any of the planes, the sphere is outside the area of interest.
Condition is wheter the signed distance d <= R (where R is particle's radius). If positive after four iterations, I add the sphere to a vector to be checked (let's call it roiSpheres).
Since I preasummed that sphere is in roi, I only check if condition is negative.

Next, create a vector of body_ids of length equal to the vector of versors. Also, another versor storing information if sphere was found on the way of versor/ray.

Run searching in parallel (each versor separately). If sphere is reached by the ray, this information is stored in the above vectors.

After parallel run, iterating again over versors to put the droplets on sphere segments.

*/





SegmentedMatSprinkler::rayToSphereInfo SegmentedMatSprinkler::oneLineAction(Vector3r lineVersor_)
{
        Vector3r d = lineVersor_.normalized(); // just making sure that versor is actually a versor
        const long size = roiSpheres.size();
        Body::id_t bId = 0; // id of the closest body (sphere) // it is initialized as zero, bot won't be used until bodyFound == true
        bool bodyFound = false;
        Real minB = 1e15; //min. distance to the closest body surface
     
        for(long i=0; i<size; i++){
            Body::id_t id = roiSpheres[i];
            const auto b = Body::byId(id, scene);
            shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
            if (!s) continue;
            SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(b->state.get()); 
            // Finding spheres that are "pierced" by the line and the closest contact point
            // vector connecting line origin and sphere center
            Vector3r AP = segState->pos - pos;
            //Distance from sphere center to line
            Real D = AP.cross(d).norm();// it should be divided by /d.norm() but d.norm() = 1, because d is a versor
            if (D < s->radius){
                    bodyFound = true;
                    // distance from projection of the center on the line (P_prim) to the contact point C
                    Real P_prim_C = sqrt(pow(s->radius,2)-pow(D,2));
                    // distance from origin to point C on the surface
                    Real B = AP.dot(d)-P_prim_C;
                    if (B < minB){
                            minB = B;
                            bId = b->id;
                    }
            } else {
                    continue;
            }
                  
        }  
        
        rayToSphereInfo singleInfo;
        singleInfo.bodyFound = bodyFound;
        singleInfo.minB = minB;
        singleInfo.bId = bId;
        return singleInfo;
}

void SegmentedMatSprinkler::action()
{       
    int counter = 0;
	if (needsInit) init();
	Real deltaTime = 0;
	if (nDone > 1 ) deltaTime = scene->time - prevTime; // delta time between runs of the sprinkler, but if it is turned on for the first time delta time should remain zero
    // compute required volume of 'droplet' carried by one line (ray)
    Real droplet = feedRate * deltaTime / (alphaResolution * betaResolution);
    
    // presecelct bodies 
    roiSpheres.clear();
    const long sizeBodies = scene->bodies->size();
 
    for(long id=0; id<sizeBodies; id++){
        bool sphereInRoi = true;
        const auto b = Body::byId(id, scene);
        shared_ptr<CohFrictMatSeg> mat = YADE_PTR_DYN_CAST<CohFrictMatSeg>(b->material);
        shared_ptr<Sphere> s = YADE_PTR_DYN_CAST<Sphere>(b->shape);
        if (!mat || !s){
            sphereInRoi = false;
            continue;
        };
        for (int i = 0; i < 4; i++){
            // select corner versors
            Vector3r v1 = cornerVersors[i];
            Vector3r v2 = Vector3r::Zero();
            if (i < 3){
                v2 = cornerVersors[i+1];
            } else {
                v2 = cornerVersors[0];
            }
            // define plane and check the signed distance
            Vector3r vn = v1.cross(v2);
            Real planeD = -(vn[0] * pos[0] + vn[1] * pos[1] + vn[2] * pos[2]);// D = -(A * Px + B * Py + C * Pz)  where P is origin of sprinkler
            Vector3r sphCenter = b->state->pos;
            Real radius = s->radius;
            Real signDist = (vn[0] * sphCenter[0] + vn[1] * sphCenter[1] + vn[2] * sphCenter[2] + planeD) / vn.norm(); 
            if (signDist > radius) sphereInRoi = false; // normal vector is pointed outwards the region of interest. We check if a tleast part of sphere is within ROI.
            
        }
        if (sphereInRoi) roiSpheres.push_back(b->id);
    }
    // prepare extra vectors and run this parallel
    long size = versors.size();
    vector<Body::id_t> bodyIdForVersors(size,0);// id of the closest body (sphere) to each versor// it is initialized as zero, bot won't be used until bodyFound == true
    vector<bool> bodyFoundForVersors(size,false);
    vector<Real> versorSphereDistances(size,-1);

#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif    
    for(long i=0; i<size; i++){  
            rayToSphereInfo singleInfo = oneLineAction(versors[i]);
            bodyIdForVersors[i] = singleInfo.bId;
            bodyFoundForVersors[i] = singleInfo.bodyFound;
            versorSphereDistances[i] = singleInfo.minB;

    }

    
    // another loop for putting the material on spheres
    for(long i=0; i<size; i++){  
        Body::id_t bId = bodyIdForVersors[i];
        bool bodyFound = bodyFoundForVersors[i];
        Real minB = versorSphereDistances[i];      
        if (bodyFound){
	        // state
	        SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get()); 
	        Vector3r surfPoint = pos + minB * versors[i]; // point on the surface
	        SegmentedState::Segment segment = segState->getSegmentFromContactPoint(surfPoint);
	        //update volume and thickness
	        Real newVolume = segment.volume + droplet;
	        segState->coatingVolume[segment.segmentPosition] = newVolume;
	        segState->coatingThickness[segment.segmentPosition] = newVolume / segment.area; 
	        counter++; // counter checking number of bodies hit 
	    } else {
	        lostVolume += droplet; // lost volume if body not found
	    }
	    releasedVolume += droplet; // always release volume     
    }    
    if (counter == 0) std::cout<< "The SegmentedMatSprinkler is active but does not reach any sphere." << std::endl;
    prevTime = scene->time;

}

/*********************************************************************************
*
* I P 2
*
*********************************************************************************/
Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::BlendedProp Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::blendMatProp(Body::id_t bId, int segmentPosition_)
{
        BlendedProp properties;
        Real thickness = 0; //
        Real minThickness;
        Real maxThickness;
	// body
	const auto b = Body::byId(bId, scene);
	// material
	shared_ptr<CohFrictMatSeg> mat = YADE_PTR_DYN_CAST<CohFrictMatSeg>(b->material);
	minThickness = mat->minThickness;
	maxThickness = mat->maxThickness;
	// state
	SegmentedState* segState = YADE_DYN_CAST<SegmentedState*>(Body::byId(bId, scene)->state.get()); 
	if (!segState->needsInit) thickness = segState->coatingThickness[segmentPosition_];
	if (thickness < minThickness || segState->needsInit){ // if state is not initialized or thickness is below min threshold set all properties primary
                properties.blendedFrictionAngle = mat->frictionAngle; 
                properties.blendedYoung = mat->young;
                properties.blendedPoisson = mat->poisson;
                properties.blendedNormalCohesion = mat->normalCohesion;
                properties.blendedShearCohesion = mat->shearCohesion;
	} else if (thickness > maxThickness){//thickness is above max threshold set all properties secondary
                properties.blendedFrictionAngle = mat->secondaryFrictionAngle; 
                properties.blendedYoung = mat->secondaryYoung;
                properties.blendedPoisson = mat->secondaryPoisson;
                properties.blendedNormalCohesion = mat->secondaryNormalCohesion;
                properties.blendedShearCohesion = mat->secondaryShearCohesion;	        
	} else { // thickness between min and max thresholds
	        Real proportion = (thickness-minThickness) / (maxThickness-minThickness);
	        
                properties.blendedFrictionAngle = mat->frictionAngle + proportion * (mat->secondaryFrictionAngle - mat->frictionAngle); 
                properties.blendedYoung = mat->young + proportion * (mat->secondaryYoung - mat->young);
                properties.blendedPoisson = mat->poisson + proportion * (mat->secondaryPoisson - mat->poisson);
                properties.blendedNormalCohesion = mat->normalCohesion + proportion * (mat->secondaryNormalCohesion - mat->normalCohesion);
                properties.blendedShearCohesion = mat->shearCohesion + proportion * (mat->secondaryShearCohesion - mat->shearCohesion);	   	
	}	
	return properties;
}



bool Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::setCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, BlendedProp matProp1, BlendedProp matProp2)
{
        CohFrictPhys* contactPhysics = YADE_DYN_CAST<CohFrictPhys*>(interaction->phys.get());
	if (not contactPhysics) {
		contactPhysics = YADE_DYN_CAST<CohFrictPhys*>(interaction->phys.get());
		if (not contactPhysics) {
			LOG_WARN("Invalid type of interaction, cohesion not set");
			return false;
		}
	}
	// if breaks
	if ((not cohesive) and not contactPhysics->cohesionBroken) {
		contactPhysics->SetBreakingState(true);
		return true;
	}
	// else bond
	if (not scene) scene = Omega::instance().getScene().get();
	auto       b1   = Body::byId(interaction->getId1(), scene);
	auto       b2   = Body::byId(interaction->getId2(), scene);
	auto       mat1 = static_cast<CohFrictMatSeg*>(b1->material.get());
	auto       mat2 = static_cast<CohFrictMatSeg*>(b2->material.get());
	const auto geom = YADE_CAST<ScGeom6D*>(interaction->geom.get());
	// determine adhesion terms with matchmakers or default formulae
	// I guess that properties blending may not work with matchmakers
	// TODO  - to discuss We should think whether matchmakers should be allowed at all, because they may be kind of useless in the segmented workflow
	/*const auto normalCohPreCalculated   = (normalCohesion) ? (*normalCohesion)(b1->id, b2->id) : math::min(matProp1.blendedNormalCohesion, matProp2.blendedNormalCohesion);
	const auto shearCohPreCalculated    = (shearCohesion) ? (*shearCohesion)(b1->id, b2->id) : math::min(matProp1.blendedShearCohesion, matProp2.blendedShearCohesion);
	const auto rollingCohPreCalculated  = (rollingCohesion) ? (*rollingCohesion)(b1->id, b2->id) : normalCohPreCalculated;
	const auto twistingCohPreCalculated = (twistingCohesion) ? (*twistingCohesion)(b1->id, b2->id) : shearCohPreCalculated;*/
        const auto normalCohPreCalculated   =  math::min(matProp1.blendedNormalCohesion, matProp2.blendedNormalCohesion);
	const auto shearCohPreCalculated    =  math::min(matProp1.blendedShearCohesion, matProp2.blendedShearCohesion);
	const auto rollingCohPreCalculated  =  normalCohPreCalculated;
	const auto twistingCohPreCalculated =  shearCohPreCalculated;
	// assign adhesions
	contactPhysics->cohesionBroken = false;
	contactPhysics->normalAdhesion = normalCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
	contactPhysics->shearAdhesion  = shearCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
	if (contactPhysics->momentRotationLaw) {
		// the max stress in pure bending is 4*M/πr^3 = 4*M/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/4)*normalAdhesion
		contactPhysics->rollingAdhesion = 0.25 * rollingCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
		// the max shear stress in pure twisting is 2*Mt/πr^3 = 2*Mt/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/2)*shearAdhesion
		contactPhysics->twistingAdhesion = 0.5 * twistingCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
	}
	geom->initRotations(*(b1->state), *(b2->state));
	contactPhysics->fragile      = (mat1->fragile || mat2->fragile);
	contactPhysics->initCohesion = false;
	return true;
}

void Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys::go(
        const shared_ptr<Material>& b1 // CohFrictMatSeg
        ,
        const shared_ptr<Material>& b2 // CohFrictMatSeg
        ,
        const shared_ptr<Interaction>& interaction)
{
	CohFrictMatSeg* mat1 = static_cast<CohFrictMatSeg*>(b1.get());
	CohFrictMatSeg* mat2 = static_cast<CohFrictMatSeg*>(b2.get());
	ScGeom6D*       geom = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	//Create cohesive interractions only once
	if (setCohesionNow && cohesionDefinitionIteration == -1) cohesionDefinitionIteration = scene->iter;
	if (setCohesionNow && cohesionDefinitionIteration != -1 && cohesionDefinitionIteration != scene->iter) {
		cohesionDefinitionIteration = -1;
		setCohesionNow              = 0;
	}

	if (geom) {
	        // The conditions may be a bit tricky, but it is mostly to preserve option to setting cohesion on existing physics without modifying CohFrictPhys (otherwise we could store information about blendedProperties / cohesion in phys)
	        bool setNewPhys = !interaction->phys;
	        if (setNewPhys)    interaction->phys            = shared_ptr<CohFrictPhys>(new CohFrictPhys());
	        CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
	        bool setCohesionOnExistingPhys = (setCohesionNow or contactPhysics->initCohesion) and (mat1->isCohesive && mat2->isCohesive);
	        
	        if (setNewPhys or setCohesionOnExistingPhys){
	                // This is common part of separable conditions below. 
		        
			const auto body1 = Body::byId(interaction->id1, scene);
			const auto body2 = Body::byId(interaction->id2, scene);

			SegmentedState* segState1 = YADE_DYN_CAST<SegmentedState*>(body1->state.get());
			SegmentedState* segState2 = YADE_DYN_CAST<SegmentedState*>(body2->state.get());
			if (!segState1 || !segState2) throw runtime_error("Ip2_CohFrictMatSeg_CohFrictMatSeg_CohFrictPhys requires that both bodies state have SegmentedState.");

	                SegmentedState::Segment segment1;
	                SegmentedState::Segment segment2;
			Vector3r midPoint = (segState1->pos + segState2->pos) / 2;// use mid point between two bodies instead of geom->contactPoint, the latter is usless because it lies on the facet, and wee need to determine on which side the contacting body is.
		        // Segment can be obtained for bodies with initialized state. If body state is not initialized, set segmentPosition to zero. 
		        ///Blending will not occur for noninitialized states anyway (in blendMatProp() function), 
		        //and will return properties as if the thickness was below minimal threshold.
		        if(segState1->needsInit){
		                segment1.segmentPosition = 0;
		        } else {
		                segment1 = segState1->getSegmentFromContactPoint(midPoint);
		        }
		        if(segState2->needsInit){
		                segment2.segmentPosition = 0;
		        } else {
		                segment2 = segState2->getSegmentFromContactPoint(midPoint);
		        }                       
                        
			// blending material properties
                        BlendedProp matProperties1 = blendMatProp(interaction->id1, segment1.segmentPosition);
                        BlendedProp matProperties2 = blendMatProp(interaction->id2, segment2.segmentPosition);
                        
		        if (setNewPhys) {

			        Real          Ea             = matProperties1.blendedYoung;
			        Real          Eb             = matProperties2.blendedYoung;
			        Real          Va             = matProperties1.blendedPoisson;
			        Real          Vb             = matProperties2.blendedPoisson;
			        Real          Da             = geom->radius1;
			        Real          Db             = geom->radius2;
			        Real          fa             = matProperties1.blendedFrictionAngle;
			        Real          fb             = matProperties2.blendedFrictionAngle;
			        Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			        // no matchmakers for now
			        //Real          frictionAngle  = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(mat1->id, mat2->id, fa, fb);
			        Real          frictionAngle  = math::min(fa, fb);
			        // harmonic average of alphas parameters
			        Real AlphaKr, AlphaKtw;
			        if (mat1->alphaKr && mat2->alphaKr) AlphaKr = 2.0 * mat1->alphaKr * mat2->alphaKr / (mat1->alphaKr + mat2->alphaKr);
			        else
				        AlphaKr = 0;
			        if (mat1->alphaKtw && mat2->alphaKtw) AlphaKtw = 2.0 * mat1->alphaKtw * mat2->alphaKtw / (mat1->alphaKtw + mat2->alphaKtw);
			        else
				        AlphaKtw = 0;

			        Real Ks;
			        if (Va && Vb)
				        Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb
				                / (Ea * Da * Va + Eb * Db * Vb); //harmonic average of two stiffnesses with ks=V*kn for each sphere
			        else
				        Ks = 0;

			        contactPhysics->kn                     = Kn;
			        contactPhysics->ks                     = Ks;
			        contactPhysics->kr                     = Da * Db * Ks * AlphaKr;
			        contactPhysics->ktw                    = Da * Db * Ks * AlphaKtw;
			        contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
			        contactPhysics->momentRotationLaw      = (mat1->momentRotationLaw && mat2->momentRotationLaw);
			        if (contactPhysics->momentRotationLaw) {
				        contactPhysics->maxRollPl  = min(mat1->etaRoll * Da, mat2->etaRoll * Db);
				        contactPhysics->maxTwistPl = min(mat1->etaTwist * Da, mat2->etaTwist * Db);
			        } else {
				        contactPhysics->maxRollPl = contactPhysics->maxTwistPl = 0;
			        }
			        if ((setCohesionOnNewContacts || setCohesionNow) && mat1->isCohesive && mat2->isCohesive) {
				        setCohesion(interaction, true, matProperties1, matProperties2); // 
			        }
		        } else { // !isNew, but if setCohesionNow, all contacts are initialized like if they were newly created
			        if (setCohesionOnExistingPhys)
				        setCohesion(interaction, true, matProperties1, matProperties2);// 
		        }
		}
	}
};


} // namespace yade
