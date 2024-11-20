MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
    rho*uSource // multiplied by rho to correct the dimmensions
    + fvm::Sp(rho*uCoeff, U) // multiplied by rho to correct the dimmensions
    + fvc::Sp(rho*uSourceDrag, uParticle) // multiplied by rho to correct the dimmensions
    + fvOptions(rho, U) // commented out for sake of simplicity
);

UEqn.relax();

fvOptions.constrain(UEqn);



  if (pimple.momentumPredictor())
  {
      solve
      (
          UEqn
       ==
          fvc::reconstruct
          (
              (
                  mixture.surfaceTensionForce()
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
              ) * mesh.magSf()
          )
      );

      fvOptions.correct(U);
  }
