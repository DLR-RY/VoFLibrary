/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an 
	unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

 

Class
    Foam::functionCurvedSurface

Description

Author
    Henning Scheufler, DLR

SourceFiles
    functionCurvedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef functionCurvedSurface_H
#define functionCurvedSurface_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"
#include "cartesianCS.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class functionCurvedSurface Declaration
\*---------------------------------------------------------------------------*/

class functionCurvedSurface
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- coord system varialbe
        const vector centre_;
        const vector normal_;
        const vector tangent_;

        //- principal curvatures
        const scalar k1_;

        const scalar k2_;

        cartesianCS coordSys_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionCurvedSurface(const functionCurvedSurface&);

        //- Disallow default bitwise assignment
        void operator=(const functionCurvedSurface&);


public:

    //- Runtime type information
    TypeName("functionCurvedSurface");


    // Constructors

        //- Construct from components
        functionCurvedSurface
        (
            const vector centre,
            const vector normal,
            const vector tangent,
            const scalar k1,
            const scalar k2
        );

        //- Construct from dictionary (used by implicitFunctions)
        functionCurvedSurface
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~functionCurvedSurface();

    virtual scalar value(const vector p) const
    {
        point localp = coordSys_.localPosition(p);
        scalar x = localp.x();
        scalar y = localp.y();
        scalar angle = 0;
        if(x != 0)
        {
            angle = (atan(y/x));
//            Info << "angle " << angle << endl;
//            Info << "x " << x << endl;
//            Info << "y " << y << endl;
        }


        scalar kn =  k1_*sqr(cos(angle)) + k2_*sqr(sin(angle));
        vector curvCentre( vector(0,0,1)*1/kn);
//        Info << "curvCentre " << curvCentre << endl;
//        Info << "normal_ " << normal_ << endl;
//        Info << "trangent_ " << tangent_ << endl;
//        Info << "mag(localp-curvCentre)- radius " << mag(localp-curvCentre)- radius << endl;
        return -mag(localp-curvCentre) + 1/kn;
    }

    virtual vector grad(const vector p) const
    {  
        point localp = coordSys_.localPosition(p);
        scalar x = localp.x();
        scalar y = localp.y();
        scalar angle = 0;
        if(x != 0)
        {
            angle = (atan(y/x));
        }


        scalar kn =  k1_*sqr(cos(angle)) + k2_*sqr(sin(angle));
        vector curvCentre( vector(0,0,1)*1/kn);
        return coordSys_.globalPosition(curvCentre) - p;
    }

    virtual scalar distanceToSurfaces(const vector p) const
    {

        return mag(value(p));
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
