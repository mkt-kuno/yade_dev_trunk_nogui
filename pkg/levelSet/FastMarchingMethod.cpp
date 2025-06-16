/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM

#include <pkg/levelSet/FastMarchingMethod.hpp>
#include <pkg/levelSet/ShopLS.hpp>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(FastMarchingMethod);
YADE_PLUGIN((FastMarchingMethod));

void FastMarchingMethod::iniStates()
{
	// initializes all states to far
	int nGPx(grid->nGP[0]), nGPy(grid->nGP[1]), nGPz(grid->nGP[2]); // how many grid pts we re working with
	// reserving memory for gpStates and ingredients:
	gpStates.reserve(nGPx);
	vector<vector<gpState>> stateInPlane;
	stateInPlane.reserve(nGPy);
	vector<gpState> stateAlongZaxis;
	stateAlongZaxis.reserve(nGPz);
	for (int xInd = 0; xInd < nGPx; xInd++) {
		stateInPlane.clear(); // https://stackoverflow.com/questions/19189699/clearing-a-vector-or-defining-a-new-vector-which-one-is-faster
		for (int yInd = 0; yInd < nGPy; yInd++) {
			stateAlongZaxis.clear();
			for (int zInd = 0; zInd < nGPz; zInd++) {
				stateAlongZaxis.push_back(farState);
			}
			stateInPlane.push_back(stateAlongZaxis);
		}
		gpStates.push_back(stateInPlane);
	}
}

void FastMarchingMethod::iniFront(bool exterior)
{
	// Identifying the boundary conditions: the very initial front band on a given side (exterior?), that will serve as first known gp
	int  nGPx(grid->nGP[0]), nGPy(grid->nGP[1]), nGPz(grid->nGP[2]);
	Real phiVal;
	for (int xInd = 0; xInd < nGPx; xInd++) {
		for (int yInd = 0; yInd < nGPy; yInd++) {
			for (int zInd = 0; zInd < nGPz; zInd++) {
				phiVal = phiField[xInd][yInd][zInd];
				if (math::isfinite(phiVal)) {
					if ((phiField[xInd][yInd][zInd] >= 0 && exterior) || (phiField[xInd][yInd][zInd] <= 0 && !exterior)) {
						knownTmp.push_back(Vector3i(xInd, yInd, zInd));
						gpStates[xInd][yInd][zInd] = knownState;
					}
				}
			}
		}
	} // closing gp loop
}

void FastMarchingMethod::confirm(int xInd, int yInd, int zInd, Real phiVal, bool exterior, bool checkState = true)
{
	if ((checkState) && (gpStates[xInd][yInd][zInd] != trialState)) LOG_ERROR("How comes ?? Current status is " << gpStates[xInd][yInd][zInd]);
	knownTmp.push_back(Vector3i(xInd, yInd, zInd)); // saving the gp among known ones
	                                                // it was already removed from trials in loopTrials()
	phiField[xInd][yInd][zInd] = phiVal;
	gpStates[xInd][yInd][zInd] = knownState;
	trializeFromKnown(xInd, yInd, zInd, exterior);
}

void FastMarchingMethod::printNeighbValues(int i, int j, int k) const
{
	string phiVal;
	for (int neigh = 0; neigh < 6; neigh++) {
		phiVal += " "
		        + std::to_string(double(
		                phiAtNgbr(neigh, i, j, k))); // no yade::math::to_string that would accommodate the precision-variable Real type at the moment
	}
	LOG_WARN(
	        "While looking at " << i << " " << j << " " << k << ", we have:\n- " << int(gpStates[i - 1][j][k]) << " and " << int(gpStates[i + 1][j][k])
	                            << " along x axis\n-   " << int(gpStates[i][j - 1][k]) << " and " << int(gpStates[i][j + 1][k]) << " along y axis\n-   "
	                            << int(gpStates[i][j][k - 1]) << " and " << int(gpStates[i][j][k + 1])
	                            << " along z axis\nWith phi = " << phiVal); // just "C regular casts" int() are fine here
}

void FastMarchingMethod::trializeFromKnown(int xInd, int yInd, int zInd, bool exterior)
{
	// we mark the 6 (or less, in edge cases) neighbors of just-confirmed (xInd,yInd,zInd) gp as trial
	int nGPx(grid->nGP[0]), nGPy(grid->nGP[1]), nGPz(grid->nGP[2]);
	if (xInd > 0)                                     // looking at the x- neighbor, if possible
		trialize(xInd - 1, yInd, zInd, exterior); // test whether that gp actually needs to be trialized will be performed therein
	if (xInd < nGPx - 1)                              // the x+ neighbor, if possible
		trialize(xInd + 1, yInd, zInd, exterior);
	if (yInd > 0) trialize(xInd, yInd - 1, zInd, exterior);        // y- neighbor if possible
	if (yInd < nGPy - 1) trialize(xInd, yInd + 1, zInd, exterior); // y+
	if (zInd > 0) trialize(xInd, yInd, zInd - 1, exterior);        // z-
	if (zInd < nGPz - 1) trialize(xInd, yInd, zInd + 1, exterior); // z+
}

void FastMarchingMethod::trialize(int xInd, int yInd, int zInd, bool exterior)
{
	if ((gpStates[xInd][yInd][zInd] != knownState)                                                         // do not touch known points !
	    && ((exterior && phiField[xInd][yInd][zInd] > 0) || (!exterior && phiField[xInd][yInd][zInd] < 0)) // let s stick on the given side
	) {
		updatePhi(
		        xInd,
		        yInd,
		        zInd,
		        exterior); // maybe that guy already had a finite phi value (was already in trials). But there has been recent changes in the neighbourhood in terms of known gp, and we thus recompute phi.
		if (gpStates[xInd][yInd][zInd] != trialState) { // putting the gp in trials after its distance has been computed (more logical this way)
			gpStates[xInd][yInd][zInd] = trialState;
			trials.push_back(Vector3i(xInd, yInd, zInd));
			if(heapSort)
				std::push_heap(trials.begin(),trials.end(),furthestAway(*this)); // will keep a heap structure for trials where *trials.begin() = trials.front() is the closest to the surface since furthestAway(a,b) is False iff a is the closest
		}
		// in case it was already in trial state, maybe the above updatePhi has changed its value and deteriorated the heap nature of trials. Calling e.g.
		// std::make_heap(trials.begin(),trials.end(),furthestAway(*this));
		// in a else here would make time costs totally explode though, cancelling by far the whole point of heap-sorting (482 s vs 1 s wo the line, vs 131 s before heap-sorting; for ~ 88^3 gp)
		// manual checks with the parent in the heap (at index k/2 where k is the index of current gp) could still be attempted though
	}
}

Real FastMarchingMethod::phiAtNgbr(int idx, int i, int j, int k) const
{
	// returns phi value at some neighbor of i j k, defined by idx
	switch (idx) {
		case 0: return phiField[i - 1][j][k];
		case 1: return phiField[i + 1][j][k];
		case 2: return phiField[i][j - 1][k];
		case 3: return phiField[i][j + 1][k];
		case 4: return phiField[i][j][k - 1];
		case 5: return phiField[i][j][k + 1];
		default: LOG_ERROR(idx << " can not be a neighbor (has to be between 0 and 5)"); return std::numeric_limits<Real>::infinity();
	}
}

Real FastMarchingMethod::phiWhenKnown(int i, int j, int k, bool exterior) const
{
	// for some gp i,j,k returns either phiValue if known or +/- infinity if not (so that we won't use elsewhere the phiValue of that guy)
	Real ret;
	if (gpStates[i][j][k] == knownState) ret = phiField[i][j][k];
	else
		ret = exterior ? std::numeric_limits<Real>::infinity() : -std::numeric_limits<Real>::infinity();
	return ret;
}
std::pair<vector<Real>, Real> FastMarchingMethod::surroundings(int i, int j, int k, bool exterior) const
{
	// will return a pair describing the surroundings, through:
	// first, the minimum known values
	// second, the associated discriminant
	vector<Real> knownSurrVal;
	knownSurrVal.reserve(3);
	Vector3i nGP(grid->nGP);
	Real     neigh[3]; // neighbor values, possibly infinite if not known
	if (i == 0) neigh[0] = phiWhenKnown(i + 1, j, k, exterior);
	else if (i == nGP[0] - 1)
		neigh[0] = phiWhenKnown(i - 1, j, k, exterior);
	else
		neigh[0] = exterior ? math::min(phiWhenKnown(i - 1, j, k, true), phiWhenKnown(i + 1, j, k, true))
		                    : math::max(phiWhenKnown(i - 1, j, k, false), phiWhenKnown(i + 1, j, k, false));
	if (j == 0) neigh[1] = phiWhenKnown(i, j + 1, k, exterior);
	else if (j == nGP[1] - 1)
		neigh[1] = phiWhenKnown(i, j - 1, k, exterior);
	else
		neigh[1] = exterior ? math::min(phiWhenKnown(i, j - 1, k, exterior), phiWhenKnown(i, j + 1, k, exterior))
		                    : math::max(phiWhenKnown(i, j - 1, k, exterior), phiWhenKnown(i, j + 1, k, exterior));
	if (k == 0) neigh[2] = phiWhenKnown(i, j, k + 1, exterior);
	else if (k == nGP[2] - 1)
		neigh[2] = phiWhenKnown(i, j, k - 1, exterior);
	else
		neigh[2] = exterior ? math::min(phiWhenKnown(i, j, k - 1, exterior), phiWhenKnown(i, j, k + 1, exterior))
		                    : math::max(phiWhenKnown(i, j, k - 1, exterior), phiWhenKnown(i, j, k + 1, exterior));
	for (int cpt = 0; cpt < 3; cpt++) {
		if (math::isfinite(neigh[cpt])) knownSurrVal.push_back(neigh[cpt]);
	}
	Real deltaPr = 0;
	switch (knownSurrVal.size()) {
		case 1:
			deltaPr = -1; // no need for a discriminant in 1D propagation anyway
			break;
		case 2:
			deltaPr = eikDiscr(
			        speed == 1 ? grid->spacing : grid->spacing * ShopLS::grad_fioRose(grid->gridPoint(i, j, k)).norm(),
			        knownSurrVal[0],
			        knownSurrVal[1]);
			break;
		case 3:
			deltaPr = eikDiscr(
			        speed == 1 ? grid->spacing : grid->spacing * ShopLS::grad_fioRose(grid->gridPoint(i, j, k)).norm(),
			        knownSurrVal[0],
			        knownSurrVal[1],
			        knownSurrVal[2]);
			break;
		default: LOG_ERROR("Unexpected case with " << knownSurrVal.size() << " known neighbors");
	}
	return std::make_pair(knownSurrVal, deltaPr);
}

Real FastMarchingMethod::eikDiscr(Real spac, Real m0, Real m1) const
{
	// reduced discriminant for the 2D version of the Eikonal equation
	return 2 * spac * spac - pow(m0 - m1, 2);
}

Real FastMarchingMethod::eikDiscr(Real spac, Real m0, Real m1, Real m2) const
{
	// reduced discriminant for the 3D version of the Eikonal equation
	return 3 * spac * spac - (pow(m0 - m1, 2) + pow(m0 - m2, 2) + pow(m1 - m2, 2));
}

void FastMarchingMethod::updatePhi(int i, int j, int k, bool exterior)
{
	// compute a phi-value at i,j,k from surrounding known values
	const auto   sdgs = surroundings(i, j, k, exterior);
	vector<Real> knownPhi(sdgs.first);
	Real         deltaPr(sdgs.second), spac(grid->spacing);
	if (speed != 1) spac *= ShopLS::grad_fioRose(grid->gridPoint(i, j, k)).norm();
	int nKnown(knownPhi.size());
	switch (nKnown) {
		case 0: LOG_ERROR("Gridpoint " << i << " " << j << " " << k << " goes through updatePhi wo any known gp"); break;
		case 1: // 1D propagation along some axis, from phi=knownPhi[0]
			LOG_INFO("1D propagation at " << i << " " << j << " " << k);
			phiField[i][j][k] = exterior ? knownPhi[0] + spac : knownPhi[0] - spac;
			break;
		case 2: { // braces necessary to enclose scope of m0,m1
			Real m0(knownPhi[0]), m1(knownPhi[1]);
			if (deltaPr >= 0) { // 2D propagation
				LOG_INFO("2D propagation at " << i << " " << j << " " << k);
				phiField[i][j][k] = phiFromEik(m0, m1, deltaPr, exterior);
			} else { // switching to 1D propagation
				LOG_INFO("1D propagation (downgraded from 2 known axes) at " << i << " " << j << " " << k);
				phiField[i][j][k] = exterior ? math::min(m0, m1) + spac : math::max(m0, m1) - spac;
			}
		} break; // break can not come before braces
		case 3: {
			Real m0(knownPhi[0]), m1(knownPhi[1]), m2(knownPhi[2]);
			if (deltaPr >= 0) { // 3D propagation
				LOG_INFO("3D propagation at " << i << " " << j << " " << k);
				phiField[i][j][k] = phiFromEik(m0, m1, m2, deltaPr, exterior);
			} else { // downgrading to 2D propagation
				Real         twoDdiscr[3] = { eikDiscr(spac, m0, m1), eikDiscr(spac, m0, m2), eikDiscr(spac, m1, m2) };
				vector<Real> possiblePhi;
				possiblePhi.reserve(3);
				for (int iBis = 0; iBis < 3; iBis++) {
					if (twoDdiscr[iBis] >= 0)
						possiblePhi.push_back(
						        phiFromEik(knownPhi[iBis / 2], knownPhi[iBis > 1 ? 2 : iBis + 1], twoDdiscr[iBis], exterior));
				}
				if (possiblePhi.size()) { // there was at least one positive discriminant
					LOG_INFO("2D propagation downgraded from 3D at " << i << " " << j << " " << k);
					phiField[i][j][k] = exterior
					        ? *std::min_element(possiblePhi.begin(), possiblePhi.end())
					        : *std::max_element(
					                possiblePhi.begin(),
					                possiblePhi.end()); // something simpler and faster than going through min_element ?
				} else { // no positive discriminant (and a empty possiblePhi) downgrading even further to 1D propagation
					LOG_INFO("1D propagation downgraded twice from 3D at " << i << " " << j << " " << k);
					phiField[i][j][k] = exterior ? math::min({m0, m1, m2}) + spac
					                             : math::max({m0, m1, m2}) - spac;
				}
			}
		} break;
		default: LOG_ERROR("Unexpected case of " << nKnown << " known neighbors around gp " << i << " " << j << " " << k);
	}
	if ((phiField[i][j][k] < 0 && exterior) || (phiField[i][j][k] > 0 && !exterior))
		LOG_ERROR(
		        "We finally assigned phi = " << phiField[i][j][k] << " to " << i << " " << j << " " << k << " supposed to be in the "
		                                     << (exterior ? "exterior" : "interior") << ".");
}

Real FastMarchingMethod::phiFromEik(Real m0, Real m1, Real disc, bool exterior) const
{ // 2D version
	return exterior ? (m0 + m1 + sqrt(disc)) / 2 : (m0 + m1 - sqrt(disc)) / 2;
}

Real FastMarchingMethod::phiFromEik(Real m0, Real m1, Real m2, Real disc, bool exterior) const
{ // 3D version
	return exterior ? (m0 + m1 + m2 + sqrt(disc)) / 3 : (m0 + m1 + m2 - sqrt(disc)) / 3;
}

vector<vector<vector<Real>>> FastMarchingMethod::phi()
{
	// computes and returns phiField. It remains to be checked what happens for repeated executions..
	iniStates(); // gpStates has now a correct format and is full of farState
	if (phiIni.size() == 0) LOG_FATAL("Empty (none given ?) phiIni in FastMarchingMethod");
	phiField = phiIni;
	for (int side = 0; side < 2; side++) {
		iniFront(side); // known now has the initial front, with correct values and states
		LOG_INFO("After iniFront on the " << (side ? "exterior" : "interior") << ", we now have: " << knownTmp.size() << " known gp to propagate from");
		for (unsigned int gpKnown = 0; gpKnown < knownTmp.size(); gpKnown++) {
			Vector3i idxOfKnown(knownTmp[gpKnown]);
			trializeFromKnown(
			        knownTmp[gpKnown][0],
			        knownTmp[gpKnown][1],
			        knownTmp[gpKnown][2],
			        side); // this time we define some narrow band with trial gridpoints.
		}
		LOG_INFO("And just after, we have " << trials.size() << " in a neighbour narrowband");
		loopTrials(side);
		known.insert(known.end(), knownTmp.begin(), knownTmp.end());
		knownTmp.clear(); // we need to start from a fresh known on side #2, but there is no real need to touch to the states, on the other hand
		LOG_INFO((side ? "Exterior" : "Interior") << " just completed, we now have " << known.size() << " known gp in total");
	}
	return phiField;
}

void FastMarchingMethod::loopTrials(bool exterior)
{
	// Some variables that will often change below:
	Vector3i closest; // the gp in trials that is the closest-to-surface one (as detected below)
	Real closestPhi; // its distance value
	Vector3i trialGP; // a random gp in trials (if !heapSort)
	vector<Vector3i>::iterator closestIt; // tentative iterator to the closest-to-surface gp in trials
	
	if(heapSort)
		std::make_heap(trials.begin(), trials.end(), furthestAway(*this));
	while (trials.size() != 0) {
		// taking out of trials the closest-to-surface gridpoint, one at a time
		LOG_DEBUG(
		        "A new iteration on the " << (exterior ? "exterior" : "interior") << ", knownTmp currently has " << knownTmp.size()
		                                  << " elements, and we have " << trials.size() << " trial points");

		if(heapSort) {
			closest = trials.front(); // by heap property
			std::pop_heap(trials.begin(),trials.end(),furthestAway(*this));
			trials.pop_back(); // would fit into confirm, but it would require passing the iterator to that function
			closestPhi = phiField[closest[0]][closest[1]][closest[2]];
		}
		else { // brute-force (and much slower) minimum search of the closest value/gp
			//  starting to look at the 1st one as a random possible initialization for that minimu:
			closestIt = trials.begin();
			closestPhi = phiField[trials[0][0]][trials[0][1]][trials[0][2]];
		    Real currPhiValue(closestPhi); // value at phi for any gp, that will often change below
			for (vector<Vector3i>::iterator trialIt = trials.begin()++; trialIt != trials.end(); trialIt++) { // starting from begin()++ is intended
				trialGP      = *trialIt;
				currPhiValue = phiField[trialGP[0]][trialGP[1]][trialGP[2]];
				if ((exterior && currPhiValue == std::numeric_limits<Real>::infinity())
			    	|| (!exterior && currPhiValue == -std::numeric_limits<Real>::infinity())) {
					LOG_ERROR("Skipping GP " << trialGP << " in the loop because it still carries an +/- infinite value");
					continue;
				}
				if ((exterior && currPhiValue < closestPhi)
			    	|| (!exterior && currPhiValue > closestPhi)) { // use a unique math::abs test ? What about short-circuit effects though ?
					closestIt  = trialIt;
					closestPhi = currPhiValue;
				}
			} // minimum now found
			closest = *closestIt;
			// removing it from trials, following https://stackoverflow.com/a/4442529/9864634 and a swap/pop_back workflow was verified to be as fast (or slightly faster) than a erase-remove idiom (which is intended to avoid relocation of elements in trials if invoking erase() not at its end)
			std::swap(*closestIt, trials.back()); // NB: std::swap(closest,trials.back()); is different and worse..
			trials.pop_back(); // would fit into confirm, but it would require passing the iterator to that function
		}
		// LOG_DEBUG("About to confirm "<<closest)
		confirm(closest[0], closest[1], closest[2], closestPhi, exterior); // this will propagate the information
		LOG_DEBUG("In that iteration, " << closest << " with phi = " << closestPhi << " was found to be the closest to the interface\n\n");
	}
}
} // namespace yade
#endif // YADE_LS_DEM
