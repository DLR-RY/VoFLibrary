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
    Foam::implicitFunction::impCurvedSurface

Description
    creates a curved surface with two main curvature axis

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    impCurvedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef impCurvedSurface_H
#define impCurvedSurface_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"
#include "cartesianCS.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impCurvedSurface Declaration
\*---------------------------------------------------------------------------*/

class impCurvedSurface
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
        impCurvedSurface(const impCurvedSurface&);

        //- Disallow default bitwise assignment
        void operator=(const impCurvedSurface&);


public:

    //- Runtime type information
    TypeName("curvedSurface");


    // Constructors

        //- Construct from components
        impCurvedSurface
        (
            const vector centre,
            const vector normal,
            const vector tangent,
            const scalar k1,
            const scalar k2
        );

        //- Construct from dictionary (used by implicitFunctions)
        impCurvedSurface
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impCurvedSurface();

    virtual scalar value(const vector p) 
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
        return -mag(localp-curvCentre) + 1/kn;
    }

    virtual vector grad(const vector p) 
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

    virtual scalar distanceToSurfaces(const vector p) 
    {
        return mag(value(p));
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
