/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::cutCell

Description
    Base Class for cutting a cell, celli, of an fvMesh, mesh_, at its intersection
    with a surface


Author
    Johan Roenby, DHI, all rights reserved.
    Modifications: Henning Scheufler, DLR, all rights reserved.


SourceFiles
    cutCell.C

\*---------------------------------------------------------------------------*/

#ifndef cutCell_H
#define cutCell_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutCell Declaration
\*---------------------------------------------------------------------------*/

class cutCell
{
    private:
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;



    // protected Member Functions
    protected:

            //- calculates volume and centre of the cutted cell
            void calcCellData
            (
                    const DynamicList<point>& cutFaceCentres,
                    const DynamicList<vector>& cutFaceAreas,
                    vector& subCellCentre,
                    scalar& subCellVolume
            );

            //- calculates area and centre of the cutting face
            void calcGeomDataCutFace
            (
                    const DynamicList<DynamicList<point> >& faceEdges,
                    const vector& subCellCentre,
                    vector& faceArea,
                    vector& faceCentre  
            );

            //- calculates the point of the cutting face
            void calcIsoFacePointsFromEdges
            (
                    const vector&  faceArea,
                    const vector&  faceCentre,
                    const DynamicList<DynamicList<point> > & faceEdges,
                    DynamicList<point>& facePoints
            );


    public:

    // Constructors

        //- Construct from fvMesh
        cutCell(const fvMesh&);


    // Member functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
