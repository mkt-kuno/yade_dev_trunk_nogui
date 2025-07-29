/*************************************************************************
*  2022 DLH van der Haven, dannyvdhaven@gmail.com                        *
*  This program is free software, see file LICENSE for details.          *
*  Code here is based on the Gl1_PotentialParticle.*pp by CWBoon 2015    *
*  for the implementation of potential particles in YADE.  				 *
*************************************************************************/
#pragma once
#ifdef YADE_LS_DEM
#ifdef YADE_OPENGL

#include <lib/opengl/GLUtils.hpp>
#include <lib/opengl/OpenGLWrapper.hpp>
#include <pkg/common/GLDrawFunctors.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/levelSet/LevelSet.hpp>
#include <pkg/levelSet/ShopLS.hpp> // Check if this one is really needed.
#include <vector>

namespace yade { // Cannot have #include directive inside.

class Gl1_LevelSet : public GlShapeFunctor {
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;

	// clang-format off
		YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_LevelSet,GlShapeFunctor,"Renders :yref:`LevelSet` object",
			((bool,recompute,false,,"Whether to recompute the triangulation every time it is rendered."))
			((bool,wire,false,,"Only show wireframe."))
			((bool,surfNodes,false,,"Render :yref:`surface nodes<LevelSet.surfNodes>`."))
			((bool,surface,true,,"Render particle surface as a zero contour of the :yref:`distance field<LevelSet.distField>` using the Marching Cubes algorithm."))
		);
	// clang-format on
	RENDERS(LevelSet);
};
REGISTER_SERIALIZABLE(Gl1_LevelSet);
} // namespace yade

#endif // YADE_OPENGL
#endif //YADE_LS_DEM
