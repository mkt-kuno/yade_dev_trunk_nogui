    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     == 
        rho*uSource // multiplied by rho to correct the dimmensions
        + fvm::Sp(rho*uCoeff, U) // multiplied by rho to correct the dimmensions
        + fvc::Sp(rho*uSourceDrag, uParticle) // multiplied by rho to correct the dimmensions //TODO: Return to check the sign of this term
        // fvOptions(rho, U) // commented out for sake of simplicity
    );

    UEqn.relax();

    // fvOptions.constrain(UEqn);


// - Note: All the momentum predictor stuff hass been commented out
//just to simplify the code a bit in the development phase.


// volScalarField rAU(1.0/UEqn.A());
// surfaceScalarField rAUcf("Dp", fvc::interpolate(rAU));

// surfaceScalarField phiForces
// (
//    fvc::flux(rAU*uSource) + rAUf*(g & mesh.Sf())
// );

//     if (pimple.momentumPredictor())
//     {
//         solve
//         (
//             UEqn
//          ==
//             fvc::reconstruct
//             (
//                 (
//                     mixture.surfaceTensionForce()
//                   - ghf*fvc::snGrad(rho)
//                   - fvc::snGrad(p_rgh)
//                 ) * mesh.magSf()
//                 +
//                 phiForces/rAUf - fvc::snGrad(p)*mesh.magSf()
//             )
//         );

//         fvOptions.correct(U);
    // }
